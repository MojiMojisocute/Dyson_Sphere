#ifndef SHARED_MEMORY_DATA_H
#define SHARED_MEMORY_DATA_H

#include <cstdint>
#include <atomic>
#include <chrono>

// Shared memory structure for solar data
struct SolarDataSharedMemory {
    // Header information
    std::atomic<uint64_t> timestamp_ns;        // Nanosecond timestamp
    std::atomic<uint32_t> sequence_number;     // Data sequence counter
    std::atomic<bool> data_valid;              // Data validity flag
    std::atomic<double> update_frequency_hz;   // Update rate
    
    // 1. Solar Radiation Data
    std::atomic<double> solar_flux_Wm2;        // Total solar flux (W/m²)
    std::atomic<double> spectral_flux_UV;      // UV flux (W/m²)
    std::atomic<double> spectral_flux_VIS;     // Visible light flux (W/m²)
    std::atomic<double> spectral_flux_IR;      // Infrared flux (W/m²)
    std::atomic<double> photon_flux;           // Photon flux (photons/m²/s)
    std::atomic<double> radiation_pressure_Pa; // Radiation pressure (Pa)
    
    // 2. Solar Wind & Plasma Data
    std::atomic<double> solar_wind_density;    // Particle density (#/m³)
    std::atomic<double> solar_wind_velocity;   // Solar wind velocity (m/s)
    std::atomic<double> solar_wind_temperature;// Solar wind temperature (K)
    std::atomic<double> solar_wind_pressure_Pa;// Solar wind pressure (Pa)
    std::atomic<double> proton_flux;           // Proton flux (#/m²/s)
    std::atomic<double> electron_flux;         // Electron flux (#/m²/s)
    std::atomic<double> alpha_particle_flux;   // Alpha particle flux (#/m²/s)
    
    // 3. Magnetic Field Data
    std::atomic<double> magnetic_field_strength_T;  // Magnetic field strength (T)
    std::atomic<double> magnetic_field_radial;      // Radial component (T)
    std::atomic<double> magnetic_field_azimuthal;   // Azimuthal component (T)
    std::atomic<double> magnetic_field_polar;       // Polar component (T)
    std::atomic<double> magnetic_fluctuation;       // Field fluctuation (%)
    
    // 4. Solar Activity Events
    std::atomic<double> solar_flare_active;    // Solar flare activity (0-1)
    std::atomic<double> flare_energy_J;        // Flare energy (J)
    std::atomic<double> flare_duration_s;      // Flare duration (s)
    std::atomic<double> cme_active;            // CME activity (0-1)
    std::atomic<double> cme_velocity;          // CME velocity (m/s)
    std::atomic<double> cme_mass_kg;           // CME mass (kg)
    std::atomic<double> cme_magnetic_field_T;  // CME magnetic field (T)
    
    // 5. Thermal & Environmental Data
    std::atomic<double> equilibrium_temperature_K;     // Equilibrium temperature (K)
    std::atomic<double> max_operating_temperature_K;   // Max operating temperature (K)
    std::atomic<double> thermal_stress_MPa;            // Thermal stress (MPa)
    std::atomic<double> micrometeorite_flux;           // Micrometeorite flux (kg/m²/s)
    std::atomic<double> debris_density;                // Debris density (#/m³)
    
    // 6. Radiation Damage Data
    std::atomic<double> radiation_dose_rate_Sv;    // Radiation dose rate (Sv/s)
    std::atomic<double> sputtering_yield;          // Sputtering yield (atoms/ion)
    std::atomic<double> erosion_rate_m_s;          // Erosion rate (m/s)
    std::atomic<double> displacement_rate;         // Displacement damage rate (#/m²/s)
    std::atomic<double> degradation_rate_per_year; // Degradation rate (%/year)
    
    // 7. Orbital Dynamics Data
    std::atomic<double> orbital_perturbation_accel;   // Orbital perturbation (m/s²)
    std::atomic<double> poynting_robertson_drag;      // P-R drag (m/s²)
    std::atomic<double> orbital_stability_index;     // Stability index (0-1)
    std::atomic<double> resonance_risk;              // Resonance risk (0-1)
    
    // 8. Control Parameters
    std::atomic<double> energy_boost_factor;      // Energy boost factor (0.5-2.0)
    std::atomic<double> cooling_power_factor;     // Cooling power factor (0.5-2.0)
    std::atomic<double> orbit_adjustment_delta;   // Orbit adjustment (-0.1 to +0.1)
    std::atomic<double> attitude_control;         // Attitude control (radians)
    std::atomic<double> degradation_mitigation;   // Degradation mitigation (0-1)
    
    // Additional derived parameters
    std::atomic<double> solar_activity_index;     // Overall solar activity (0-2)
    std::atomic<double> distance_from_sun_m;      // Current distance from sun (m)
    std::atomic<double> solar_cycle_phase;        // Solar cycle phase (0-1)
    
    // Constructor
    SolarDataSharedMemory() {
        // Initialize all atomic values to zero/false
        timestamp_ns.store(0);
        sequence_number.store(0);
        data_valid.store(false);
        update_frequency_hz.store(1.0);
        
        // Initialize all data fields to zero
        solar_flux_Wm2.store(0.0);
        spectral_flux_UV.store(0.0);
        spectral_flux_VIS.store(0.0);
        spectral_flux_IR.store(0.0);
        photon_flux.store(0.0);
        radiation_pressure_Pa.store(0.0);
        
        solar_wind_density.store(0.0);
        solar_wind_velocity.store(0.0);
        solar_wind_temperature.store(0.0);
        solar_wind_pressure_Pa.store(0.0);
        proton_flux.store(0.0);
        electron_flux.store(0.0);
        alpha_particle_flux.store(0.0);
        
        magnetic_field_strength_T.store(0.0);
        magnetic_field_radial.store(0.0);
        magnetic_field_azimuthal.store(0.0);
        magnetic_field_polar.store(0.0);
        magnetic_fluctuation.store(0.0);
        
        solar_flare_active.store(0.0);
        flare_energy_J.store(0.0);
        flare_duration_s.store(0.0);
        cme_active.store(0.0);
        cme_velocity.store(0.0);
        cme_mass_kg.store(0.0);
        cme_magnetic_field_T.store(0.0);
        
        equilibrium_temperature_K.store(0.0);
        max_operating_temperature_K.store(0.0);
        thermal_stress_MPa.store(0.0);
        micrometeorite_flux.store(0.0);
        debris_density.store(0.0);
        
        radiation_dose_rate_Sv.store(0.0);
        sputtering_yield.store(0.0);
        erosion_rate_m_s.store(0.0);
        displacement_rate.store(0.0);
        degradation_rate_per_year.store(0.0);
        
        orbital_perturbation_accel.store(0.0);
        poynting_robertson_drag.store(0.0);
        orbital_stability_index.store(1.0);
        resonance_risk.store(0.0);
        
        energy_boost_factor.store(1.0);
        cooling_power_factor.store(1.0);
        orbit_adjustment_delta.store(0.0);
        attitude_control.store(0.0);
        degradation_mitigation.store(0.5);
        
        solar_activity_index.store(1.0);
        distance_from_sun_m.store(1.495978707e11); // 1 AU
        solar_cycle_phase.store(0.0);
    }
};

// Constants for shared memory management
namespace SharedMemoryConstants {
    constexpr const char* SOLAR_DATA_SHM_NAME = "/solar_calculations_data";
    constexpr size_t SOLAR_DATA_SHM_SIZE = sizeof(SolarDataSharedMemory);
    constexpr int MAX_READERS = 10;                    // Maximum number of reader processes
    constexpr double DEFAULT_UPDATE_RATE = 10.0;      // Hz
    constexpr uint64_t STALE_DATA_THRESHOLD_NS = 1000000000; // 1 second in nanoseconds
}

#endif // SHARED_MEMORY_DATA_H