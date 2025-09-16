#ifndef SHARED_MEMORY_MANAGER_H
#define SHARED_MEMORY_MANAGER_H

#include "shared_memory_data.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <semaphore.h>
#include <stdexcept>
#include <memory>
#include <chrono>
#include <cstring>
#include <iostream>

class SharedMemoryManager {
public:
    enum class AccessMode {
        PRODUCER,  // Write access (Solar Calculations)
        CONSUMER   // Read access (Dyson Sphere, Fuzzy Controller, etc.)
    };

private:
    AccessMode mode_;
    int shm_fd_;
    SolarDataSharedMemory* data_ptr_;
    sem_t* write_semaphore_;
    sem_t* read_semaphore_;
    bool is_creator_;
    
    // Semaphore names
    static constexpr const char* WRITE_SEM_NAME = "/solar_write_sem";
    static constexpr const char* READ_SEM_NAME = "/solar_read_sem";

public:
    explicit SharedMemoryManager(AccessMode mode) 
        : mode_(mode), shm_fd_(-1), data_ptr_(nullptr), 
          write_semaphore_(SEM_FAILED), read_semaphore_(SEM_FAILED), is_creator_(false) {
        
        initialize_shared_memory();
        initialize_semaphores();
    }
    
    ~SharedMemoryManager() {
        cleanup();
    }

    // Delete copy constructor and assignment operator
    SharedMemoryManager(const SharedMemoryManager&) = delete;
    SharedMemoryManager& operator=(const SharedMemoryManager&) = delete;
    
    // Move constructor and assignment
    SharedMemoryManager(SharedMemoryManager&& other) noexcept
        : mode_(other.mode_), shm_fd_(other.shm_fd_), data_ptr_(other.data_ptr_),
          write_semaphore_(other.write_semaphore_), read_semaphore_(other.read_semaphore_),
          is_creator_(other.is_creator_) {
        
        other.shm_fd_ = -1;
        other.data_ptr_ = nullptr;
        other.write_semaphore_ = SEM_FAILED;
        other.read_semaphore_ = SEM_FAILED;
        other.is_creator_ = false;
    }
    
    SharedMemoryManager& operator=(SharedMemoryManager&& other) noexcept {
        if (this != &other) {
            cleanup();
            
            mode_ = other.mode_;
            shm_fd_ = other.shm_fd_;
            data_ptr_ = other.data_ptr_;
            write_semaphore_ = other.write_semaphore_;
            read_semaphore_ = other.read_semaphore_;
            is_creator_ = other.is_creator_;
            
            other.shm_fd_ = -1;
            other.data_ptr_ = nullptr;
            other.write_semaphore_ = SEM_FAILED;
            other.read_semaphore_ = SEM_FAILED;
            other.is_creator_ = false;
        }
        return *this;
    }

private:
    void initialize_shared_memory() {
        if (mode_ == AccessMode::PRODUCER) {
            // Create shared memory segment
            shm_fd_ = shm_open(SharedMemoryConstants::SOLAR_DATA_SHM_NAME, 
                              O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
            is_creator_ = true;
            
            if (shm_fd_ == -1) {
                throw std::runtime_error("Failed to create shared memory segment");
            }
            
            // Set size of shared memory
            if (ftruncate(shm_fd_, SharedMemoryConstants::SOLAR_DATA_SHM_SIZE) == -1) {
                close(shm_fd_);
                shm_unlink(SharedMemoryConstants::SOLAR_DATA_SHM_NAME);
                throw std::runtime_error("Failed to set shared memory size");
            }
        } else {
            // Open existing shared memory segment
            shm_fd_ = shm_open(SharedMemoryConstants::SOLAR_DATA_SHM_NAME, O_RDONLY, 0);
            
            if (shm_fd_ == -1) {
                throw std::runtime_error("Failed to open shared memory segment");
            }
        }
        
        // Map shared memory
        int prot = (mode_ == AccessMode::PRODUCER) ? (PROT_READ | PROT_WRITE) : PROT_READ;
        data_ptr_ = static_cast<SolarDataSharedMemory*>(
            mmap(nullptr, SharedMemoryConstants::SOLAR_DATA_SHM_SIZE, prot, MAP_SHARED, shm_fd_, 0));
        
        if (data_ptr_ == MAP_FAILED) {
            close(shm_fd_);
            if (is_creator_) {
                shm_unlink(SharedMemoryConstants::SOLAR_DATA_SHM_NAME);
            }
            throw std::runtime_error("Failed to map shared memory");
        }
        
        // Initialize data structure if we're the creator
        if (is_creator_) {
            new (data_ptr_) SolarDataSharedMemory();
        }
    }
    
    void initialize_semaphores() {
        if (mode_ == AccessMode::PRODUCER) {
            // Create semaphores
            write_semaphore_ = sem_open(WRITE_SEM_NAME, O_CREAT, S_IRUSR | S_IWUSR, 1);
            read_semaphore_ = sem_open(READ_SEM_NAME, O_CREAT, S_IRUSR | S_IWUSR, 0);
        } else {
            // Open existing semaphores
            write_semaphore_ = sem_open(WRITE_SEM_NAME, 0);
            read_semaphore_ = sem_open(READ_SEM_NAME, 0);
        }
        
        if (write_semaphore_ == SEM_FAILED || read_semaphore_ == SEM_FAILED) {
            cleanup();
            throw std::runtime_error("Failed to initialize semaphores");
        }
    }
    
    void cleanup() {
        if (data_ptr_ != nullptr && data_ptr_ != MAP_FAILED) {
            munmap(data_ptr_, SharedMemoryConstants::SOLAR_DATA_SHM_SIZE);
            data_ptr_ = nullptr;
        }
        
        if (shm_fd_ != -1) {
            close(shm_fd_);
            shm_fd_ = -1;
        }
        
        if (write_semaphore_ != SEM_FAILED) {
            sem_close(write_semaphore_);
            if (is_creator_) {
                sem_unlink(WRITE_SEM_NAME);
            }
            write_semaphore_ = SEM_FAILED;
        }
        
        if (read_semaphore_ != SEM_FAILED) {
            sem_close(read_semaphore_);
            if (is_creator_) {
                sem_unlink(READ_SEM_NAME);
            }
            read_semaphore_ = SEM_FAILED;
        }
        
        if (is_creator_ && shm_fd_ != -1) {
            shm_unlink(SharedMemoryConstants::SOLAR_DATA_SHM_NAME);
        }
    }

public:
    // Producer methods (Solar Calculations)
    bool write_data_begin() {
        if (mode_ != AccessMode::PRODUCER) {
            return false;
        }
        
        struct timespec timeout;
        clock_gettime(CLOCK_REALTIME, &timeout);
        timeout.tv_sec += 1; // 1 second timeout
        
        return sem_timedwait(write_semaphore_, &timeout) == 0;
    }
    
    void write_data_end() {
        if (mode_ == AccessMode::PRODUCER) {
            // Update timestamp and sequence number
            auto now = std::chrono::high_resolution_clock::now();
            auto timestamp = std::chrono::duration_cast<std::chrono::nanoseconds>(
                now.time_since_epoch()).count();
            
            data_ptr_->timestamp_ns.store(timestamp);
            data_ptr_->sequence_number.fetch_add(1);
            data_ptr_->data_valid.store(true);
            
            sem_post(read_semaphore_);
            sem_post(write_semaphore_);
        }
    }
    
    // Consumer methods (Dyson Sphere components)
    bool read_data_begin(std::chrono::milliseconds timeout = std::chrono::milliseconds(100)) {
        if (mode_ != AccessMode::CONSUMER) {
            return false;
        }
        
        struct timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        ts.tv_nsec += timeout.count() * 1000000; // Convert ms to ns
        if (ts.tv_nsec >= 1000000000) {
            ts.tv_sec += ts.tv_nsec / 1000000000;
            ts.tv_nsec %= 1000000000;
        }
        
        return sem_timedwait(read_semaphore_, &ts) == 0;
    }
    
    void read_data_end() {
        if (mode_ == AccessMode::CONSUMER) {
            sem_post(read_semaphore_);
        }
    }
    
    // Data access methods
    SolarDataSharedMemory* get_data() {
        return data_ptr_;
    }
    
    const SolarDataSharedMemory* get_data() const {
        return data_ptr_;
    }
    
    // Utility methods
    bool is_data_valid() const {
        if (!data_ptr_) return false;
        
        bool valid = data_ptr_->data_valid.load();
        
        // Check if data is stale
        auto now = std::chrono::high_resolution_clock::now();
        auto current_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
            now.time_since_epoch()).count();
        auto data_time = data_ptr_->timestamp_ns.load();
        
        if (current_time - data_time > SharedMemoryConstants::STALE_DATA_THRESHOLD_NS) {
            valid = false;
        }
        
        return valid;
    }
    
    uint32_t get_sequence_number() const {
        return data_ptr_ ? data_ptr_->sequence_number.load() : 0;
    }
    
    double get_update_frequency() const {
        return data_ptr_ ? data_ptr_->update_frequency_hz.load() : 0.0;
    }
    
    uint64_t get_timestamp() const {
        return data_ptr_ ? data_ptr_->timestamp_ns.load() : 0;
    }
    
    // Static cleanup method for emergency situations
    static void cleanup_shared_resources() {
        shm_unlink(SharedMemoryConstants::SOLAR_DATA_SHM_NAME);
        sem_unlink(WRITE_SEM_NAME);
        sem_unlink(READ_SEM_NAME);
    }
    
    // Debug and monitoring methods
    void print_status() const {
        if (!data_ptr_) {
            std::cout << "SharedMemoryManager: No data pointer" << std::endl;
            return;
        }
        
        std::cout << "=== Shared Memory Status ===" << std::endl;
        std::cout << "Mode: " << (mode_ == AccessMode::PRODUCER ? "Producer" : "Consumer") << std::endl;
        std::cout << "Data Valid: " << (is_data_valid() ? "Yes" : "No") << std::endl;
        std::cout << "Sequence Number: " << get_sequence_number() << std::endl;
        std::cout << "Update Frequency: " << get_update_frequency() << " Hz" << std::endl;
        std::cout << "Timestamp: " << get_timestamp() << " ns" << std::endl;
        std::cout << "Solar Flux: " << data_ptr_->solar_flux_Wm2.load() << " W/m²" << std::endl;
        std::cout << "Distance from Sun: " << data_ptr_->distance_from_sun_m.load() / 1.496e11 << " AU" << std::endl;
        std::cout << "Solar Activity: " << data_ptr_->solar_activity_index.load() << std::endl;
        std::cout << "=============================" << std::endl;
    }
};

#endif // SHARED_MEMORY_MANAGER_H