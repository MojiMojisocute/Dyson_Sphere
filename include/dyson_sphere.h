#ifndef DYSON_SPHERE_H
#define DYSON_SPHERE_H

#include "physicsConstants.h"
#include "shared_memory_manager.h"
#include "solar_calculations.h"
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <atomic>
#include <mutex>
#include <chrono>

// Dyson Sphere component types
enum class ComponentType {
    ENERGY_COLLECTOR,
    THERMAL_RADIATOR, 
    ATTITUDE_THRUSTER,
    COMMUNICATION_ARRAY,
    STRUCTURAL_SUPPORT,
    MAINTENANCE_SYSTEM
};

// Material properties structure
struct MaterialProperties {
    std::string name;
    double density;              // kg/m³
    double youngs_modulus;       // Pa
    double tensile_strength;     // Pa
    double melting_point;        // K
    double thermal_conductivity; // W/m·K
    double specific_heat;        // J/kg·K
    double thermal_expansion;    // 1/K
    double emissivity;          // dimensionless
    double absorptivity;        // dimensionless
    double radiation_resistance; // Sv⁻¹
};

// Structural element for finite element analysis
struct StructuralElement {
    uint32_t id;
    std::array<double, 3> position;    // m
    std::array<double, 3> velocity;    // m/s
    std::array<double, 3> force;       // N
    double mass;                       // kg
    double temperature;                // K
    double stress_von_mises;           // Pa
    double strain_energy;              // J
    bool is_failed;                    // failure flag
};

// Energy system component
struct EnergySystem {
    double energy_stored;              // J
    double energy_capacity_max;        // J
    double power_generation_rate;      // W
    double power_consumption_rate;     // W
    double conversion_efficiency;      // dimensionless
    double charge_rate_max;           // W
    double discharge_rate_max;        // W
    double self_discharge_rate;       // 1/s
    double temperature;               // K
};

// Thermal management system
struct ThermalSystem {
    double heat_capacity_total;       // J/K
    double cooling_power_max;         // W
    double temperature_current;       // K
    double temperature_optimal;       // K
    double temperature_critical;      // K
    double thermal_conductance;       // W/K
    double radiator_area;            // m²
    double radiator_emissivity;      // dimensionless
    std::vector<double> temperature_distribution; // K per element
};

// Orbital dynamics system
struct OrbitalSystem {
    double semi_major_axis;           // m
    double eccentricity;              // dimensionless
    double inclination;               // rad
    double longitude_ascending_node;  // rad
    double argument_periapsis;        // rad
    double true_anomaly;             // rad
    double orbital_period;           // s
    double orbital_velocity;         // m/s
    std::array<double, 3> position_heliocentric; // m
    std::array<double, 3> velocity_heliocentric; // m/s
    std::array<double, 3> acceleration_total;    // m/s²
};

// Degradation tracking system
struct DegradationSystem {
    double sputtering_damage;         // atoms removed per m²
    double radiation_damage;          // displacement per atom (dpa)
    double thermal_fatigue_cycles;    // number of cycles
    double micrometeorite_impacts;    // impacts per m²
    double structural_integrity;      // 0-1 scale
    double performance_degradation;   // 0-1 scale (1 = no degradation)
    double estimated_lifetime;        // s
    std::vector<double> local_damage; // damage per structural element
};

// Fuzzy control output structure
struct FuzzyControlOutput {
    double energy_boost_factor;       // 0.5-2.0
    double cooling_power_factor;      // 0.5-2.0
    double orbit_adjustment_delta;    // -0.1 to +0.1
    double attitude_control;          // radians
    double degradation_mitigation;    // 0-1
    double emergency_shutdown;        // 0-1 (boolean-like)
    double maintenance_priority;      // 0-1
};

class DysonSphere {
private:
    // Core parameters
    double sphere_radius_;            // m
    double sphere_mass_;              // kg
    double shell_thickness_;          // m
    double surface_area_;             // m²
    uint32_t num_elements_;           // number of structural elements
    
    // Material and construction
    MaterialProperties material_;
    std::vector<StructuralElement> structure_;
    
    // Physical systems
    EnergySystem energy_system_;
    ThermalSystem thermal_system_;
    OrbitalSystem orbital_system_;
    DegradationSystem degradation_system_;
    
    // External data interface
    std::unique_ptr<SharedMemoryManager> shm_manager_;
    SolarCalculations solar_calc_;
    
    // Time and simulation
    std::atomic<double> simulation_time_;     // s
    std::atomic<bool> simulation_running_;
    std::mutex state_mutex_;
    
    // Physics calculation methods
    void calculate_structural_forces();
    void calculate_thermal_distribution();
    void calculate_orbital_dynamics(double dt);
    void calculate_energy_balance(double dt);
    void calculate_degradation_effects(double dt);
    void update_material_properties();
    
    // Internal physics helpers
    double calculate_gravitational_force(const std::array<double, 3>& pos) const;
    double calculate_radiation_pressure_force(double distance) const;
    double calculate_solar_wind_pressure_force(double distance) const;
    std::array<double, 3> calculate_poynting_robertson_drag() const;
    double calculate_tidal_acceleration(const std::array<double, 3>& pos) const;
    
    // Thermal physics
    double calculate_equilibrium_temperature(double solar_flux) const;
    double calculate_radiative_cooling_power(double temperature, double area) const;
    double calculate_thermal_stress(double temp_gradient, double expansion_coeff) const;
    
    // Material science
    double calculate_von_mises_stress(const std::array<double, 6>& stress_tensor) const;
    double calculate_fatigue_damage(double stress_amplitude, int cycles) const;
    bool check_material_failure(double stress, double temperature) const;
    
    // Damage models
    double calculate_sputtering_rate(const SolarDataSharedMemory* solar_data) const;
    double calculate_radiation_damage_rate(const SolarDataSharedMemory* solar_data) const;
    double calculate_thermal_cycling_damage(double temp_range) const;
    
public:
    // Constructor
    DysonSphere(double radius, double shell_thickness, const std::string& material_name,
                uint32_t num_structural_elements = 10000);
    
    // Destructor
    ~DysonSphere();
    
    // Initialization methods
    bool initialize_sphere();
    bool initialize_shared_memory();
    void set_initial_orbit(double distance_from_sun, double eccentricity = 0.0);
    void set_material_properties(const std::string& material_name);
    
    // Main simulation methods
    void start_simulation(double time_step = 0.1); // s
    void stop_simulation();
    void update_physics(double dt);
    bool is_simulation_running() const;
    
    // Data access methods
    const EnergySystem& get_energy_system() const { return energy_system_; }
    const ThermalSystem& get_thermal_system() const { return thermal_system_; }
    const OrbitalSystem& get_orbital_system() const { return orbital_system_; }
    const DegradationSystem& get_degradation_system() const { return degradation_system_; }
    const MaterialProperties& get_material() const { return material_; }
    
    // Physical property getters
    double get_radius() const { return sphere_radius_; }
    double get_mass() const { return sphere_mass_; }
    double get_surface_area() const { return surface_area_; }
    double get_shell_thickness() const { return shell_thickness_; }
    double get_current_distance_from_sun() const;
    double get_current_solar_flux() const;
    
    // Performance metrics
    double get_energy_efficiency() const;           // W_out / W_in
    double get_mass_efficiency() const;             // W / kg
    double get_structural_integrity() const;        // 0-1 scale
    double get_thermal_stability() const;           // 0-1 scale
    double get_orbital_stability() const;           // 0-1 scale
    double get_overall_health() const;              // 0-1 scale
    
    // Critical status checks
    bool is_structurally_stable() const;
    bool is_thermally_stable() const;
    bool is_orbitally_stable() const;
    bool requires_emergency_action() const;
    double get_estimated_lifetime() const;          // years
    
    // Control interface methods (for Fuzzy Controller)
    void apply_control_inputs(const FuzzyControlOutput& control_output);
    void set_energy_boost_factor(double factor);   // 0.5 - 2.0
    void set_cooling_power_factor(double factor);  // 0.5 - 2.0
    void adjust_orbital_position(double delta);    // -0.1 to +0.1
    void set_attitude_control(double angle_rad);   // radians
    void activate_degradation_mitigation(double level); // 0-1
    
    // Fuzzy system interface - ข้อมูลที่จะส่งไปให้ Fuzzy Controller
    struct FuzzyInputs {
        // Energy parameters
        double energy_level_normalized;           // 0-1 (current/max)
        double power_generation_efficiency;      // 0-1
        double energy_balance_rate;              // W (positive = charging)
        
        // Thermal parameters  
        double temperature_normalized;           // 0-1 (current/critical)
        double thermal_gradient_max;            // K/m
        double cooling_system_efficiency;       // 0-1
        
        // Structural parameters
        double max_stress_normalized;           // 0-1 (current/failure)
        double structural_integrity;            // 0-1
        double fatigue_damage_level;           // 0-1
        
        // Environmental parameters
        double solar_flux_normalized;          // 0-1 (current/nominal)
        double solar_activity_level;          // 0-1
        double radiation_dose_rate;           // Sv/s
        
        // Orbital parameters
        double orbital_stability_index;        // 0-1
        double distance_from_optimal;         // AU (deviation)
        double perturbation_magnitude;        // m/s²
        
        // Degradation parameters
        double sputtering_rate_normalized;    // 0-1
        double overall_degradation_rate;      // %/year
        double estimated_lifetime_years;      // years
        
        // System status
        double emergency_level;               // 0-1 (0=normal, 1=critical)
        double maintenance_urgency;           // 0-1
        double performance_degradation;       // 0-1 (0=no degradation)
    };
    
    FuzzyInputs get_fuzzy_inputs();
    
    // Monitoring and diagnostics
    void print_status_report();
    void export_telemetry_data(const std::string& filename) const;
    std::vector<double> get_stress_distribution() const;
    std::vector<double> get_temperature_distribution() const;
    
    // Static utility methods
    static MaterialProperties get_material_properties(const std::string& material_name);
    static double calculate_sphere_mass(double radius, double thickness, double density);
    static double calculate_sphere_surface_area(double radius);
    static std::vector<std::string> get_available_materials();
    
    // Constants for physics calculations
    static constexpr double MIN_ORBITAL_DISTANCE = 0.1 * PhysicsConstants::AU;  // 0.1 AU
    static constexpr double MAX_ORBITAL_DISTANCE = 10.0 * PhysicsConstants::AU; // 10 AU
    static constexpr double SAFETY_FACTOR_STRUCTURAL = 2.0;  // 2x safety margin
    static constexpr double SAFETY_FACTOR_THERMAL = 1.5;     // 1.5x safety margin
    static constexpr double CRITICAL_STRESS_RATIO = 0.8;     // 80% of ultimate strength
    static constexpr double CRITICAL_TEMPERATURE_RATIO = 0.9; // 90% of melting point
};

// Helper function declarations
namespace DysonSphereUtils {
    double interpolate_material_property(double temperature, const std::vector<double>& temps,
                                       const std::vector<double>& values);
    std::array<double, 3> spherical_to_cartesian(double r, double theta, double phi);
    std::array<double, 3> cartesian_to_spherical(const std::array<double, 3>& pos);
    double calculate_orbital_period(double semi_major_axis, double central_mass);
    double calculate_escape_velocity(double mass, double radius);
}

#endif // DYSON_SPHERE_H