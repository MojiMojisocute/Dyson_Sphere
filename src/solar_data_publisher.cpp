#include "solar_calculations.h"
#include "shared_memory_manager.h"
#include <thread>
#include <iostream>
#include <iomanip>
#include <signal.h>

class SolarDataPublisher {
private:
    SolarCalculations solar_calc_;
    SharedMemoryManager shm_manager_;
    std::atomic<bool> running_;
    std::thread publisher_thread_;
    
    // Current simulation parameters
    double current_distance_;        // Current distance from sun (m)
    double current_time_years_;      // Current simulation time (years)
    double sphere_mass_;             // Dyson sphere mass (kg)
    double sphere_radius_;           // Dyson sphere radius (m)
    double material_type_;           // Material type index
    std::string current_material_;   // Current material name

public:
    SolarDataPublisher() 
        : shm_manager_(SharedMemoryManager::AccessMode::PRODUCER),
          running_(false),
          current_distance_(1.495978707e11),  // 1 AU
          current_time_years_(0.0),
          sphere_mass_(1e15),                 // 1 petagram
          sphere_radius_(1e6),                // 1000 km
          material_type_(0),
          current_material_("aluminum") {
    }
    
    ~SolarDataPublisher() {
        stop();
    }
    
    void start(double update_frequency = 10.0) {
        if (running_.load()) {
            std::cout << "Publisher already running!" << std::endl;
            return;
        }
        
        running_.store(true);
        
        publisher_thread_ = std::thread([this, update_frequency]() {
            publish_loop(update_frequency);
        });
        
        std::cout << "Solar data publisher started at " << update_frequency << " Hz" << std::endl;
    }
    
    void stop() {
        if (!running_.load()) {
            return;
        }
        
        running_.store(false);
        
        if (publisher_thread_.joinable()) {
            publisher_thread_.join();
        }
        
        std::cout << "Solar data publisher stopped" << std::endl;
    }
    
    // Parameter setters for external control
    void set_distance(double distance_au) {
        current_distance_ = distance_au * PhysicsConstants::AU;
    }
    
    void set_simulation_time(double time_years) {
        current_time_years_ = time_years;
    }
    
    void set_sphere_properties(double mass_kg, double radius_m) {
        sphere_mass_ = mass_kg;
        sphere_radius_ = radius_m;
    }
    
    void set_material(const std::string& material) {
        current_material_ = material;
        // Map material to index for faster lookup
        if (material == "aluminum") material_type_ = 0;
        else if (material == "silicon") material_type_ = 1;
        else if (material == "iron") material_type_ = 2;
        else if (material == "carbon") material_type_ = 3;
        else if (material == "tungsten") material_type_ = 4;
        else material_type_ = 0; // Default to aluminum
    }

private:
    void publish_loop(double frequency) {
        const auto sleep_duration = std::chrono::microseconds(
            static_cast<int64_t>(1e6 / frequency));
        
        uint32_t iteration = 0;
        
        while (running_.load()) {
            auto start_time = std::chrono::high_resolution_clock::now();
            
            // Acquire write lock
            if (shm_manager_.write_data_begin()) {
                try {
                    update_solar_data();
                    shm_manager_.write_data_end();
                } catch (const std::exception& e) {
                    std::cerr << "Error updating solar data: " << e.what() << std::endl;
                    shm_manager_.write_data_end(); // Always release lock
                }
            } else {
                std::cerr << "Failed to acquire write lock" << std::endl;
            }
            
            // Print status periodically
            if (++iteration % static_cast<uint32_t>(frequency) == 0) {
                print_status();
            }
            
            // Sleep to maintain frequency
            auto end_time = std::chrono::high_resolution_clock::now();
            auto elapsed = end_time - start_time;
            
            if (elapsed < sleep_duration) {
                std::this_thread::sleep_for(sleep_duration - elapsed);
            }
            
            // Advance simulation time (assuming real-time scaling)
            current_time_years_ += 1.0 / (frequency * 365.25 * 24.0 * 3600.0);
        }
    }
    
    void update_solar_data() {
        auto* data = shm_manager_.get_data();
        if (!data) return;
        
        // Update frequency for monitoring
        data->update_frequency_hz.store(10.0);
        
        // Current simulation parameters
        data->distance_from_sun_m.store(current_distance_);
        data->solar_cycle_phase.store(std::fmod(current_time_years_ / 11.0, 1.0));
        
        // 1. Solar Radiation Data
        const double solar_flux = solar_calc_.radiation_flux_at_distance(current_distance_);
        data->solar_flux_Wm2.store(solar_flux);
        
        // Spectral band calculations
        auto spectral_bands = solar_calc_.detailed_spectral_bands();
        data->spectral_flux_UV.store(solar_flux * (spectral_bands[0] + spectral_bands[1] + spectral_bands[2]));
        data->spectral_flux_VIS.store(solar_flux * spectral_bands[3]);
        data->spectral_flux_IR.store(solar_flux * spectral_bands[4]);
        
        // Photon flux calculation
        const double avg_photon_energy = PhysicsConstants::h * PhysicsConstants::c / 550e-9; // Green light
        data->photon_flux.store(solar_flux / avg_photon_energy);
        
        // Radiation pressure
        data->radiation_pressure_Pa.store(solar_calc_.radiation_pressure(current_distance_));
        
        // 2. Solar Wind & Plasma Data
        const auto wind_data = solar_calc_.calculate_solar_wind(current_distance_);
        data->solar_wind_density.store(wind_data.density);
        data->solar_wind_velocity.store(wind_data.velocity);
        data->solar_wind_temperature.store(wind_data.temperature);
        data->solar_wind_pressure_Pa.store(wind_data.dynamic_pressure);
        
        // Particle fluxes (simplified distribution)
        data->proton_flux.store(wind_data.density * wind_data.velocity * 0.90); // 90% protons
        data->electron_flux.store(wind_data.density * wind_data.velocity * 0.98); // Charge balance
        data->alpha_particle_flux.store(wind_data.density * wind_data.velocity * 0.08); // 8% alphas
        
        // 3. Magnetic Field Data
        const double mag_field_strength = solar_calc_.magnetic_field_strength(current_distance_, current_time_years_);
        data->magnetic_field_strength_T.store(mag_field_strength);
        
        auto mag_field_vector = solar_calc_.magnetic_field_vector(current_distance_, 0.0, 0.0, current_time_years_);
        data->magnetic_field_radial.store(mag_field_vector[0]);
        data->magnetic_field_azimuthal.store(mag_field_vector[1]);
        data->magnetic_field_polar.store(mag_field_vector[2]);
        
        // Magnetic field fluctuations (simplified model)
        const double solar_activity = solar_calc_.interpolate_solar_activity(current_time_years_);
        data->magnetic_fluctuation.store(5.0 + 15.0 * solar_activity); // 5-20% fluctuation
        
        // 4. Solar Activity Events
        data->solar_activity_index.store(solar_activity);
        
        // Solar flare modeling
        const double flare_rate = solar_calc_.flare_occurrence_rate(solar_activity);
        const bool flare_active = (std::rand() / double(RAND_MAX)) < (flare_rate / 3600.0); // per second
        data->solar_flare_active.store(flare_active ? 1.0 : 0.0);
        
        if (flare_active) {
            auto flare_data = solar_calc_.generate_solar_flare();
            data->flare_energy_J.store(flare_data.energy);
            data->flare_duration_s.store(flare_data.duration);
        } else {
            data->flare_energy_J.store(0.0);
            data->flare_duration_s.store(0.0);
        }
        
        // CME modeling
        const bool cme_active = solar_calc_.cme_occurrence_check(1.0/10.0, solar_activity); // Check per 0.1 hour
        data->cme_active.store(cme_active ? 1.0 : 0.0);
        
        if (cme_active) {
            auto cme_data = solar_calc_.calculate_cme_properties(current_distance_);
            data->cme_velocity.store(cme_data.velocity);
            data->cme_mass_kg.store(cme_data.mass);
            data->cme_magnetic_field_T.store(cme_data.magnetic_field);
        } else {
            data->cme_velocity.store(0.0);
            data->cme_mass_kg.store(0.0);
            data->cme_magnetic_field_T.store(0.0);
        }
        
        // 5. Thermal & Environmental Data
        const double equilibrium_temp = solar_calc_.equilibrium_temperature(current_distance_);
        data->equilibrium_temperature_K.store(equilibrium_temp);
        
        // Maximum operating temperature (material dependent)
        double max_temp = 1000.0; // Default
        if (current_material_ == "aluminum") max_temp = 900.0;
        else if (current_material_ == "silicon") max_temp = 1200.0;
        else if (current_material_ == "iron") max_temp = 1400.0;
        else if (current_material_ == "carbon") max_temp = 3000.0;
        else if (current_material_ == "tungsten") max_temp = 3500.0;
        data->max_operating_temperature_K.store(max_temp);
        
        // Thermal stress
        const double temp_change = equilibrium_temp - 300.0; // Assume 300K base temperature
        const double thermal_stress = solar_calc_.thermal_expansion_stress(temp_change, current_material_);
        data->thermal_stress_MPa.store(thermal_stress / 1e6); // Convert Pa to MPa
        
        // Environmental hazards
        const double micromet_flux = solar_calc_.micrometeorite_impact_rate(current_distance_, 1.0);
        data->micrometeorite_flux.store(micromet_flux);
        
        // Simplified debris density model
        const double debris_base = 1e-12 * std::pow(PhysicsConstants::AU / current_distance_, 2);
        data->debris_density.store(debris_base);
        
        // 6. Radiation Damage Data
        const double radiation_dose = solar_calc_.total_radiation_dose_rate(current_distance_);
        data->radiation_dose_rate_Sv.store(radiation_dose);
        
        // Material-specific sputtering and erosion
        if (!wind_data.particles.empty()) {
            const double sputtering = solar_calc_.sputtering_yield(wind_data.particles[0], current_material_);
            data->sputtering_yield.store(sputtering);
            
            const double erosion_rate = solar_calc_.plasma_erosion_rate(current_distance_, current_material_, 1.0);
            data->erosion_rate_m_s.store(erosion_rate * 1e-10); // Atoms/s to m/s (rough conversion)
        } else {
            data->sputtering_yield.store(0.0);
            data->erosion_rate_m_s.store(0.0);
        }
        
        const double displacement = solar_calc_.atomic_displacement_rate(current_distance_, current_material_);
        data->displacement_rate.store(displacement);
        
        // Overall degradation rate (combination of all effects)
        const double base_degradation = 0.1; // 0.1% per year baseline
        const double radiation_factor = std::min(radiation_dose / 1e-6, 10.0); // Scale factor
        const double activity_factor = 1.0 + solar_activity;
        data->degradation_rate_per_year.store(base_degradation * radiation_factor * activity_factor);
        
        // 7. Orbital Dynamics Data
        const double orbital_perturb = solar_calc_.orbital_perturbation_acceleration(current_distance_, sphere_mass_);
        data->orbital_perturbation_accel.store(orbital_perturb);
        
        // Poynting-Robertson drag
        const double particle_density = 2700.0; // kg/m³, typical for aluminum
        auto pr_drag = solar_calc_.poynting_robertson_drag(current_distance_, sphere_radius_, particle_density);
        data->poynting_robertson_drag.store(std::sqrt(pr_drag[0]*pr_drag[0] + pr_drag[1]*pr_drag[1]));
        
        // Orbital stability (simplified model)
        const double stability = std::max(0.0, 1.0 - orbital_perturb / 1e-6);
        data->orbital_stability_index.store(stability);
        
        // Resonance risk (distance-dependent)
        const double resonance_distances[] = {0.5, 1.0, 1.5, 2.0, 3.0}; // AU
        double min_resonance_dist = 1000.0;
        for (double res_dist : resonance_distances) {
            const double dist_au = current_distance_ / PhysicsConstants::AU;
            const double diff = std::abs(dist_au - res_dist);
            min_resonance_dist = std::min(min_resonance_dist, diff);
        }
        data->resonance_risk.store(std::max(0.0, 1.0 - min_resonance_dist / 0.1));
        
        // 8. Control Parameters (can be modified by external systems)
        // These are typically set by the fuzzy controller or other control systems
        // For now, we'll compute some basic adaptive parameters
        
        // Energy boost factor (higher when solar flux is low)
        const double nominal_flux = 1361.0; // W/m² at 1 AU
        const double energy_boost = std::min(2.0, std::max(0.5, nominal_flux / solar_flux));
        data->energy_boost_factor.store(energy_boost);
        
        // Cooling power factor (higher when temperature is high)
        const double cooling_boost = std::min(2.0, std::max(0.5, equilibrium_temp / 400.0));
        data->cooling_power_factor.store(cooling_boost);
        
        // Orbit adjustment (positive moves outward, negative inward)
        double orbit_adjustment = 0.0;
        if (equilibrium_temp > max_temp * 0.9) orbit_adjustment = 0.05; // Move outward if too hot
        if (solar_flux < nominal_flux * 0.5) orbit_adjustment = -0.05; // Move inward if too little flux
        data->orbit_adjustment_delta.store(orbit_adjustment);
        
        // Attitude control (point toward sun for maximum collection)
        const double attitude = 0.0; // Simplified: always sun-pointing
        data->attitude_control.store(attitude);
        
        // Degradation mitigation (increase with activity and damage)
        const double mitigation = std::min(1.0, 0.3 + 0.4 * solar_activity + 0.3 * radiation_factor / 10.0);
        data->degradation_mitigation.store(mitigation);
    }
    
    void print_status() {
        std::cout << "\n=== Solar Data Publisher Status ===" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "Distance: " << current_distance_ / PhysicsConstants::AU << " AU" << std::endl;
        std::cout << "Simulation Time: " << current_time_years_ << " years" << std::endl;
        std::cout << "Material: " << current_material_ << std::endl;
        
        auto* data = shm_manager_.get_data();
        if (data) {
            std::cout << "Solar Flux: " << data->solar_flux_Wm2.load() << " W/m²" << std::endl;
            std::cout << "Temperature: " << data->equilibrium_temperature_K.load() << " K" << std::endl;
            std::cout << "Solar Activity: " << data->solar_activity_index.load() << std::endl;
            std::cout << "Sequence: " << data->sequence_number.load() << std::endl;
        }
        std::cout << "===================================" << std::endl;
    }
};

// Signal handler for graceful shutdown
std::unique_ptr<SolarDataPublisher> g_publisher;

void signal_handler(int signal) {
    std::cout << "\nReceived signal " << signal << ", shutting down..." << std::endl;
    if (g_publisher) {
        g_publisher->stop();
    }
    exit(0);
}

// Main function for testing
int main(int argc, char* argv[]) {
    // Setup signal handling
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);
    
    // Parse command line arguments
    double update_frequency = 10.0;
    double distance_au = 1.0;
    std::string material = "aluminum";
    
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-f" && i + 1 < argc) {
            update_frequency = std::stod(argv[++i]);
        } else if (std::string(argv[i]) == "-d" && i + 1 < argc) {
            distance_au = std::stod(argv[++i]);
        } else if (std::string(argv[i]) == "-m" && i + 1 < argc) {
            material = argv[++i];
        } else if (std::string(argv[i]) == "-h") {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "  -f frequency  Update frequency in Hz (default: 10)" << std::endl;
            std::cout << "  -d distance   Distance from sun in AU (default: 1.0)" << std::endl;
            std::cout << "  -m material   Material type (aluminum|silicon|iron|carbon|tungsten)" << std::endl;
            std::cout << "  -h           Show this help" << std::endl;
            return 0;
        }
    }
    
    try {
        // Create and configure publisher
        g_publisher = std::make_unique<SolarDataPublisher>();
        g_publisher->set_distance(distance_au);
        g_publisher->set_material(material);
        
        std::cout << "Starting Solar Data Publisher..." << std::endl;
        std::cout << "Frequency: " << update_frequency << " Hz" << std::endl;
        std::cout << "Distance: " << distance_au << " AU" << std::endl;
        std::cout << "Material: " << material << std::endl;
        std::cout << "Press Ctrl+C to stop" << std::endl;
        
        // Start publishing
        g_publisher->start(update_frequency);
        
        // Keep main thread alive
        while (true) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}