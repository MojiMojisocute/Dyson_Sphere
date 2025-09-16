#include "dyson_sphere.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <thread>
#include <numeric>
#include <mutex>

using namespace PhysicsConstants;

// Constructor
DysonSphere::DysonSphere(double radius, double shell_thickness, const std::string& material_name,
                        uint32_t num_structural_elements)
    : sphere_radius_(radius)
    , shell_thickness_(shell_thickness)
    , num_elements_(num_structural_elements)
    , simulation_time_(0.0)
    , simulation_running_(false) {
    
    // Set material properties
    set_material_properties(material_name);
    
    // Calculate basic sphere properties
    sphere_mass_ = calculate_sphere_mass(radius, shell_thickness, material_.density);
    surface_area_ = calculate_sphere_surface_area(radius);
    
    // Initialize systems with default values
    energy_system_ = {};
    thermal_system_ = {};
    orbital_system_ = {};
    degradation_system_ = {};
    
    // Initialize shared memory manager
    try {
        shm_manager_ = std::make_unique<SharedMemoryManager>(SharedMemoryManager::AccessMode::CONSUMER);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to initialize shared memory: " << e.what() << std::endl;
        shm_manager_ = nullptr;
    }
}

// Destructor
DysonSphere::~DysonSphere() {
    stop_simulation();
}

// Initialize sphere structure
bool DysonSphere::initialize_sphere() {
    std::lock_guard<std::mutex> lock(state_mutex_);
    
    // Initialize structural elements in spherical grid
    structure_.clear();
    structure_.reserve(num_elements_);
    
    // Create approximately uniform distribution on sphere surface
    const double golden_angle = M_PI * (3.0 - std::sqrt(5.0)); // Golden angle in radians
    
    for (uint32_t i = 0; i < num_elements_; ++i) {
        StructuralElement element;
        element.id = i;
        
        // Golden spiral distribution for uniform coverage
        double y = 1.0 - 2.0 * double(i) / double(num_elements_ - 1); // y from 1 to -1
        double radius_at_y = std::sqrt(1.0 - y * y);
        double theta = golden_angle * i;
        
        // Position on sphere surface (outer surface)
        element.position[0] = radius_at_y * std::cos(theta) * sphere_radius_;
        element.position[1] = radius_at_y * std::sin(theta) * sphere_radius_;
        element.position[2] = y * sphere_radius_;
        
        // Initialize other properties
        element.velocity = {0.0, 0.0, 0.0};
        element.force = {0.0, 0.0, 0.0};
        element.mass = sphere_mass_ / num_elements_; // Uniform mass distribution
        element.temperature = 300.0; // Room temperature initial
        element.stress_von_mises = 0.0;
        element.strain_energy = 0.0;
        element.is_failed = false;
        
        structure_.push_back(element);
    }
    
    // Initialize energy system
    energy_system_.energy_capacity_max = 1e24; // 1 YJ (from spec)
    energy_system_.energy_stored = energy_system_.energy_capacity_max * 0.5; // Start half charged
    energy_system_.conversion_efficiency = 0.3; // 30% solar panel efficiency
    energy_system_.charge_rate_max = 0.1 * SolarParams::luminosity;
    energy_system_.discharge_rate_max = 0.2 * SolarParams::luminosity;
    energy_system_.self_discharge_rate = 1e-6; // 0.0001% per second
    energy_system_.temperature = 300.0;
    
    // Initialize thermal system
    thermal_system_.heat_capacity_total = sphere_mass_ * material_.specific_heat;
    thermal_system_.temperature_current = 300.0;
    thermal_system_.temperature_optimal = 400.0;
    thermal_system_.temperature_critical = material_.melting_point * CRITICAL_TEMPERATURE_RATIO;
    thermal_system_.radiator_area = surface_area_;
    thermal_system_.radiator_emissivity = material_.emissivity;
    thermal_system_.cooling_power_max = sigma * surface_area_ * 
                                       std::pow(thermal_system_.temperature_critical, 4);
    thermal_system_.temperature_distribution.resize(num_elements_, 300.0);
    
    // Initialize degradation system
    degradation_system_.structural_integrity = 1.0;
    degradation_system_.performance_degradation = 1.0;
    degradation_system_.estimated_lifetime = 50.0 * year; // 50 years initial estimate
    degradation_system_.local_damage.resize(num_elements_, 0.0);
    
    return true;
}

// Set initial orbital parameters
void DysonSphere::set_initial_orbit(double distance_from_sun, double eccentricity) {
    std::lock_guard<std::mutex> lock(state_mutex_);
    
    orbital_system_.semi_major_axis = distance_from_sun;
    orbital_system_.eccentricity = eccentricity;
    orbital_system_.inclination = 0.0; // Ecliptic plane
    orbital_system_.longitude_ascending_node = 0.0;
    orbital_system_.argument_periapsis = 0.0;
    orbital_system_.true_anomaly = 0.0;
    
    // Calculate derived orbital parameters
    orbital_system_.orbital_period = DysonSphereUtils::calculate_orbital_period(
        distance_from_sun, SolarParams::mass);
    orbital_system_.orbital_velocity = 2.0 * M_PI * distance_from_sun / orbital_system_.orbital_period;
    
    // Set initial position (at aphelion for simplicity)
    orbital_system_.position_heliocentric[0] = distance_from_sun;
    orbital_system_.position_heliocentric[1] = 0.0;
    orbital_system_.position_heliocentric[2] = 0.0;
    
    // Set initial velocity (circular orbit approximation)
    orbital_system_.velocity_heliocentric[0] = 0.0;
    orbital_system_.velocity_heliocentric[1] = orbital_system_.orbital_velocity;
    orbital_system_.velocity_heliocentric[2] = 0.0;
    
    orbital_system_.acceleration_total = {0.0, 0.0, 0.0};
}

// Set material properties
void DysonSphere::set_material_properties(const std::string& material_name) {
    material_ = get_material_properties(material_name);
}

// Main physics update
void DysonSphere::update_physics(double dt) {
    if (!simulation_running_) return;
    
    std::lock_guard<std::mutex> lock(state_mutex_);
    
    // Get latest solar data
    const SolarDataSharedMemory* solar_data = nullptr;
    if (shm_manager_ && shm_manager_->is_data_valid()) {
        solar_data = shm_manager_->get_data();
    }
    
    // Update all physics subsystems
    calculate_orbital_dynamics(dt);
    calculate_energy_balance(dt);
    calculate_thermal_distribution();
    calculate_structural_forces();
    calculate_degradation_effects(dt);
    update_material_properties();
    
    double old_time = simulation_time_.load();
    simulation_time_.store(old_time + dt);
}

// Calculate orbital dynamics
void DysonSphere::calculate_orbital_dynamics(double dt) {
    // Get current position
    auto& pos = orbital_system_.position_heliocentric;
    auto& vel = orbital_system_.velocity_heliocentric;
    auto& acc = orbital_system_.acceleration_total;
    
    double distance = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    
    // Gravitational acceleration (primary)
    double grav_acc_magnitude = G * SolarParams::mass / (distance * distance);
    acc[0] = -grav_acc_magnitude * pos[0] / distance;
    acc[1] = -grav_acc_magnitude * pos[1] / distance;
    acc[2] = -grav_acc_magnitude * pos[2] / distance;
    
    // Add radiation pressure acceleration
    double rad_press_acc = calculate_radiation_pressure_force(distance) / sphere_mass_;
    double unit_vec[3] = {pos[0]/distance, pos[1]/distance, pos[2]/distance};
    acc[0] += rad_press_acc * unit_vec[0];
    acc[1] += rad_press_acc * unit_vec[1];
    acc[2] += rad_press_acc * unit_vec[2];
    
    // Add solar wind pressure
    double wind_press_acc = calculate_solar_wind_pressure_force(distance) / sphere_mass_;
    acc[0] += wind_press_acc * unit_vec[0];
    acc[1] += wind_press_acc * unit_vec[1];
    acc[2] += wind_press_acc * unit_vec[2];
    
    // Add Poynting-Robertson drag
    auto pr_drag = calculate_poynting_robertson_drag();
    acc[0] += pr_drag[0];
    acc[1] += pr_drag[1];
    acc[2] += pr_drag[2];
    
    // Integrate using Verlet integration for better stability
    pos[0] += vel[0] * dt + 0.5 * acc[0] * dt * dt;
    pos[1] += vel[1] * dt + 0.5 * acc[1] * dt * dt;
    pos[2] += vel[2] * dt + 0.5 * acc[2] * dt * dt;
    
    vel[0] += acc[0] * dt;
    vel[1] += acc[1] * dt;
    vel[2] += acc[2] * dt;
    
    // Update derived parameters
    distance = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    orbital_system_.orbital_velocity = std::sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
}

// Calculate energy balance
void DysonSphere::calculate_energy_balance(double dt) {
    double distance = get_current_distance_from_sun();
    double solar_flux = solar_calc_.radiation_flux_at_distance(distance);
    
    // Energy generation from solar collection
    double power_input = solar_flux * surface_area_ * energy_system_.conversion_efficiency;
    
    // Energy consumption (basic systems)
    double power_consumption = 1e12; // 1 TW baseline consumption for station-keeping, etc.
    
    // Thermal losses
    double thermal_power_loss = calculate_radiative_cooling_power(
        thermal_system_.temperature_current, surface_area_);
    
    // Net power balance
    double net_power = power_input - power_consumption - thermal_power_loss;
    
    // Update energy storage
    if (net_power > 0) {
        double charge_power = std::min(net_power, energy_system_.charge_rate_max);
        energy_system_.energy_stored += charge_power * dt;
    } else {
        double discharge_power = std::min(-net_power, energy_system_.discharge_rate_max);
        energy_system_.energy_stored -= discharge_power * dt;
    }
    
    // Apply self-discharge
    energy_system_.energy_stored *= (1.0 - energy_system_.self_discharge_rate * dt);
    
    // Clamp to capacity limits
    energy_system_.energy_stored = std::max(0.0, 
        std::min(energy_system_.energy_stored, energy_system_.energy_capacity_max));
    
    // Update power rates
    energy_system_.power_generation_rate = power_input;
    energy_system_.power_consumption_rate = power_consumption + thermal_power_loss;
}

// Calculate thermal distribution
void DysonSphere::calculate_thermal_distribution() {
    double distance = get_current_distance_from_sun();
    double solar_flux = solar_calc_.radiation_flux_at_distance(distance);
    
    // Calculate equilibrium temperature
    double T_equilibrium = calculate_equilibrium_temperature(solar_flux);
    
    // Simple thermal model - exponential approach to equilibrium
    double thermal_time_constant = thermal_system_.heat_capacity_total / 
                                  (sigma * surface_area_ * material_.emissivity * 4.0 * 
                                   std::pow(thermal_system_.temperature_current, 3));
    
    double dt = 0.1; // Assume 0.1s time step for thermal calculation
    double alpha = dt / thermal_time_constant;
    
    thermal_system_.temperature_current = thermal_system_.temperature_current * (1.0 - alpha) +
                                         T_equilibrium * alpha;
    
    // Update temperature distribution (simplified: uniform with small variations)
    for (size_t i = 0; i < thermal_system_.temperature_distribution.size(); ++i) {
        double local_variation = 0.05 * std::sin(2.0 * M_PI * i / thermal_system_.temperature_distribution.size());
        thermal_system_.temperature_distribution[i] = thermal_system_.temperature_current * (1.0 + local_variation);
        
        // Update structural element temperatures
        if (i < structure_.size()) {
            structure_[i].temperature = thermal_system_.temperature_distribution[i];
        }
    }
}

// Calculate structural forces and stresses
void DysonSphere::calculate_structural_forces() {
    double distance = get_current_distance_from_sun();
    
    // Calculate forces on each structural element
    for (auto& element : structure_) {
        // Reset forces
        element.force = {0.0, 0.0, 0.0};
        
        // Gravitational force from sun
        double grav_force = calculate_gravitational_force(element.position);
        double r = std::sqrt(element.position[0]*element.position[0] + 
                           element.position[1]*element.position[1] + 
                           element.position[2]*element.position[2]);
        
        element.force[0] -= grav_force * element.position[0] / r;
        element.force[1] -= grav_force * element.position[1] / r;
        element.force[2] -= grav_force * element.position[2] / r;
        
        // Radiation pressure force
        double rad_force = calculate_radiation_pressure_force(distance) / num_elements_;
        element.force[0] += rad_force * element.position[0] / r;
        element.force[1] += rad_force * element.position[1] / r;
        element.force[2] += rad_force * element.position[2] / r;
        
        // Thermal stress
        double temp_gradient = std::abs(element.temperature - thermal_system_.temperature_optimal);
        double thermal_stress = calculate_thermal_stress(temp_gradient, material_.thermal_expansion);
        
        // Convert thermal stress to equivalent force (simplified)
        double thermal_force = thermal_stress * shell_thickness_ * sphere_radius_ / num_elements_;
        
        // Calculate von Mises stress (simplified - hoop stress dominant)
        double hoop_stress = (rad_force + thermal_force) / (shell_thickness_ * 2.0 * M_PI * sphere_radius_ / num_elements_);
        element.stress_von_mises = std::abs(hoop_stress);
        
        // Check for failure
        element.is_failed = check_material_failure(element.stress_von_mises, element.temperature);
        
        // Calculate strain energy
        element.strain_energy = 0.5 * element.stress_von_mises * element.stress_von_mises / material_.youngs_modulus;
    }
}

// Calculate degradation effects
void DysonSphere::calculate_degradation_effects(double dt) {
    if (!shm_manager_ || !shm_manager_->is_data_valid()) return;
    
    const auto* solar_data = shm_manager_->get_data();
    
    // Sputtering damage
    double sputtering_rate = calculate_sputtering_rate(solar_data);
    degradation_system_.sputtering_damage += sputtering_rate * dt;
    
    // Radiation damage
    double radiation_damage_rate = calculate_radiation_damage_rate(solar_data);
    degradation_system_.radiation_damage += radiation_damage_rate * dt;
    
    // Thermal cycling damage
    double temp_range = thermal_system_.temperature_current - thermal_system_.temperature_optimal;
    double thermal_cycling_damage = calculate_thermal_cycling_damage(std::abs(temp_range));
    degradation_system_.thermal_fatigue_cycles += thermal_cycling_damage * dt;
    
    // Micrometeorite damage
    double micromet_rate = solar_data->micrometeorite_flux.load() * surface_area_;
    degradation_system_.micrometeorite_impacts += micromet_rate * dt;
    
    // Update overall structural integrity
    double total_damage = degradation_system_.sputtering_damage / 1e15 +  // Normalize
                         degradation_system_.radiation_damage / 1e10 +
                         degradation_system_.thermal_fatigue_cycles / 1e6 +
                         degradation_system_.micrometeorite_impacts / 1e12;
    
    degradation_system_.structural_integrity = std::max(0.0, 1.0 - total_damage);
    degradation_system_.performance_degradation = std::max(0.1, degradation_system_.structural_integrity);
    
    // Update estimated lifetime
    if (total_damage > 0) {
        double damage_rate_per_year = total_damage / (simulation_time_ / year);
        degradation_system_.estimated_lifetime = (1.0 - total_damage) / damage_rate_per_year * year;
    }
    
    // Update local damage for each element
    for (size_t i = 0; i < degradation_system_.local_damage.size() && i < structure_.size(); ++i) {
        double local_factor = 1.0 + 0.2 * std::sin(2.0 * M_PI * i / degradation_system_.local_damage.size());
        degradation_system_.local_damage[i] = total_damage * local_factor;
        
        // Mark severely damaged elements as failed
        if (degradation_system_.local_damage[i] > 0.8) {
            structure_[i].is_failed = true;
        }
    }
}

// Physics helper methods
double DysonSphere::calculate_gravitational_force(const std::array<double, 3>& pos) const {
    double r_squared = pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2];
    return G * SolarParams::mass * (sphere_mass_ / num_elements_) / r_squared;
}

double DysonSphere::calculate_radiation_pressure_force(double distance) const {
    double solar_flux = solar_calc_.radiation_flux_at_distance(distance);
    double radiation_pressure = solar_flux / c; // Assuming perfect absorption
    return radiation_pressure * surface_area_;
}

double DysonSphere::calculate_solar_wind_pressure_force(double distance) const {
    auto wind_data = solar_calc_.calculate_solar_wind(distance);
    return wind_data.dynamic_pressure * surface_area_;
}

std::array<double, 3> DysonSphere::calculate_poynting_robertson_drag() const {
    double distance = get_current_distance_from_sun();
    auto pr_drag = solar_calc_.poynting_robertson_drag(distance, sphere_radius_, material_.density);
    return pr_drag;
}

double DysonSphere::calculate_equilibrium_temperature(double solar_flux) const {
    double absorbed_power = solar_flux * surface_area_ * material_.absorptivity;
    double T_eq_4 = absorbed_power / (material_.emissivity * sigma * surface_area_);
    return std::pow(T_eq_4, 0.25);
}

double DysonSphere::calculate_radiative_cooling_power(double temperature, double area) const {
    return material_.emissivity * sigma * area * std::pow(temperature, 4);
}

double DysonSphere::calculate_thermal_stress(double temp_gradient, double expansion_coeff) const {
    return material_.youngs_modulus * expansion_coeff * temp_gradient;
}

double DysonSphere::calculate_von_mises_stress(const std::array<double, 6>& stress_tensor) const {
    double sx = stress_tensor[0], sy = stress_tensor[1], sz = stress_tensor[2];
    double txy = stress_tensor[3], txz = stress_tensor[4], tyz = stress_tensor[5];
    
    return std::sqrt(0.5 * ((sx-sy)*(sx-sy) + (sy-sz)*(sy-sz) + (sz-sx)*(sz-sx) + 
                           6.0 * (txy*txy + txz*txz + tyz*tyz)));
}

bool DysonSphere::check_material_failure(double stress, double temperature) const {
    // Temperature-dependent strength reduction
    double temp_factor = 1.0;
    if (temperature > 0.7 * material_.melting_point) {
        temp_factor = 1.0 - (temperature - 0.7 * material_.melting_point) / (0.3 * material_.melting_point);
    }
    
    double effective_strength = material_.tensile_strength * temp_factor / SAFETY_FACTOR_STRUCTURAL;
    return stress > effective_strength;
}

double DysonSphere::calculate_sputtering_rate(const SolarDataSharedMemory* solar_data) const {
    if (!solar_data) return 0.0;
    
    double sputtering_yield = solar_data->sputtering_yield.load();
    double proton_flux = solar_data->proton_flux.load();
    double alpha_flux = solar_data->alpha_particle_flux.load();
    
    return (proton_flux + alpha_flux * 4.0) * sputtering_yield; // Alpha particles ~4x more effective
}

double DysonSphere::calculate_radiation_damage_rate(const SolarDataSharedMemory* solar_data) const {
    if (!solar_data) return 0.0;
    
    double dose_rate = solar_data->radiation_dose_rate_Sv.load();
    double displacement_rate = solar_data->displacement_rate.load();
    
    // Convert to displacements per atom (dpa)
    double atoms_per_m3 = material_.density / (material_.density > 5000 ? 60.0 * m_u : 27.0 * m_u); // Rough estimate
    return displacement_rate / atoms_per_m3;
}

double DysonSphere::calculate_thermal_cycling_damage(double temp_range) const {
    if (temp_range < 50.0) return 0.0; // No significant damage below 50K cycles
    
    // Coffin-Manson relationship for thermal fatigue
    double cycles_to_failure = 1e6 * std::pow(temp_range / 100.0, -2.5);
    return 1.0 / cycles_to_failure; // Damage per cycle
}

// Control interface methods
void DysonSphere::apply_control_inputs(const FuzzyControlOutput& control_output) {
    std::lock_guard<std::mutex> lock(state_mutex_);
    
    set_energy_boost_factor(control_output.energy_boost_factor);
    set_cooling_power_factor(control_output.cooling_power_factor);
    adjust_orbital_position(control_output.orbit_adjustment_delta);
    set_attitude_control(control_output.attitude_control);
    activate_degradation_mitigation(control_output.degradation_mitigation);
    
    // Handle emergency shutdown
    if (control_output.emergency_shutdown > 0.5) {
        // Implement emergency procedures
        energy_system_.power_generation_rate *= 0.1; // Reduce power generation
        thermal_system_.cooling_power_max *= 2.0;    // Emergency cooling
    }
}

void DysonSphere::set_energy_boost_factor(double factor) {
    factor = std::max(0.5, std::min(2.0, factor));
    energy_system_.conversion_efficiency *= factor;
    energy_system_.charge_rate_max *= factor;
}

void DysonSphere::set_cooling_power_factor(double factor) {
    factor = std::max(0.5, std::min(2.0, factor));
    thermal_system_.cooling_power_max *= factor;
}

void DysonSphere::adjust_orbital_position(double delta) {
    delta = std::max(-0.1, std::min(0.1, delta));
    double current_distance = get_current_distance_from_sun();
    double new_distance = current_distance * (1.0 + delta);
    
    // Clamp to safe orbital range
    new_distance = std::max(MIN_ORBITAL_DISTANCE, std::min(MAX_ORBITAL_DISTANCE, new_distance));
    
    // Apply orbital adjustment (simplified - instantaneous)
    double scale_factor = new_distance / current_distance;
    orbital_system_.position_heliocentric[0] *= scale_factor;
    orbital_system_.position_heliocentric[1] *= scale_factor;
    orbital_system_.position_heliocentric[2] *= scale_factor;
    
    // Adjust velocity to maintain orbit
    double new_orbital_velocity = std::sqrt(G * SolarParams::mass / new_distance);
    double velocity_scale = new_orbital_velocity / orbital_system_.orbital_velocity;
    orbital_system_.velocity_heliocentric[0] *= velocity_scale;
    orbital_system_.velocity_heliocentric[1] *= velocity_scale;
    orbital_system_.velocity_heliocentric[2] *= velocity_scale;
}

void DysonSphere::set_attitude_control(double angle_rad) {
    // Simplified attitude control - just store the value
    // In real implementation, this would adjust solar panel orientations, etc.
    (void)angle_rad; // Suppress unused parameter warning
}

void DysonSphere::activate_degradation_mitigation(double level) {
    level = std::max(0.0, std::min(1.0, level));
    
    // Reduce degradation rates based on mitigation level
    degradation_system_.structural_integrity += level * 0.001; // Small recovery
    degradation_system_.structural_integrity = std::min(1.0, degradation_system_.structural_integrity);
    
    // Mitigation costs energy
    double mitigation_power = level * 1e11; // 100 GW at full mitigation
    energy_system_.power_consumption_rate += mitigation_power;
}

// Fuzzy system interface
DysonSphere::FuzzyInputs DysonSphere::get_fuzzy_inputs() {
    std::lock_guard<std::mutex> lock(state_mutex_);
    
    FuzzyInputs inputs;
    
    // Energy parameters
    inputs.energy_level_normalized = energy_system_.energy_stored / energy_system_.energy_capacity_max;
    inputs.power_generation_efficiency = energy_system_.power_generation_rate / 
                                        (get_current_solar_flux() * surface_area_);
    inputs.energy_balance_rate = energy_system_.power_generation_rate - energy_system_.power_consumption_rate;
    
    // Thermal parameters
    inputs.temperature_normalized = thermal_system_.temperature_current / thermal_system_.temperature_critical;
    
    // Calculate maximum temperature gradient
    if (!thermal_system_.temperature_distribution.empty()) {
        auto minmax = std::minmax_element(thermal_system_.temperature_distribution.begin(),
                                        thermal_system_.temperature_distribution.end());
        inputs.thermal_gradient_max = (*minmax.second - *minmax.first) / sphere_radius_;
    } else {
        inputs.thermal_gradient_max = 0.0;
    }
    
    inputs.cooling_system_efficiency = std::min(1.0, thermal_system_.cooling_power_max / 
                                      (sigma * surface_area_ * std::pow(thermal_system_.temperature_critical, 4)));
    
    // Structural parameters
    double max_stress = 0.0;
    for (const auto& element : structure_) {
        max_stress = std::max(max_stress, element.stress_von_mises);
    }
    inputs.max_stress_normalized = max_stress / (material_.tensile_strength / SAFETY_FACTOR_STRUCTURAL);
    inputs.structural_integrity = degradation_system_.structural_integrity;
    inputs.fatigue_damage_level = std::min(1.0, degradation_system_.thermal_fatigue_cycles / 1e6);
    
    // Environmental parameters
    double current_solar_flux = get_current_solar_flux();
    inputs.solar_flux_normalized = current_solar_flux / 1361.0; // Solar constant at 1 AU
    
    // Get solar activity if available
    if (shm_manager_ && shm_manager_->is_data_valid()) {
        const auto* solar_data = shm_manager_->get_data();
        inputs.solar_activity_level = solar_data->solar_activity_index.load() / 2.0; // Normalize 0-2 to 0-1
        inputs.radiation_dose_rate = solar_data->radiation_dose_rate_Sv.load();
    } else {
        inputs.solar_activity_level = 1.0;
        inputs.radiation_dose_rate = 1e-8; // Default background
    }
    
    // Orbital parameters
    inputs.orbital_stability_index = get_orbital_stability();
    inputs.distance_from_optimal = (get_current_distance_from_sun() - AU) / AU; // Deviation from 1 AU
    inputs.perturbation_magnitude = std::sqrt(orbital_system_.acceleration_total[0] * orbital_system_.acceleration_total[0] +
                                            orbital_system_.acceleration_total[1] * orbital_system_.acceleration_total[1] +
                                            orbital_system_.acceleration_total[2] * orbital_system_.acceleration_total[2]);
    
    // Degradation parameters
    inputs.sputtering_rate_normalized = std::min(1.0, degradation_system_.sputtering_damage / 1e15);
    inputs.overall_degradation_rate = (1.0 - degradation_system_.performance_degradation) * 100.0; // %/year
    inputs.estimated_lifetime_years = degradation_system_.estimated_lifetime / year;
    
    // System status
    inputs.emergency_level = 0.0;
    if (!is_structurally_stable()) inputs.emergency_level = std::max(inputs.emergency_level, 0.7);
    if (!is_thermally_stable()) inputs.emergency_level = std::max(inputs.emergency_level, 0.8);
    if (!is_orbitally_stable()) inputs.emergency_level = std::max(inputs.emergency_level, 0.6);
    
    inputs.maintenance_urgency = 1.0 - degradation_system_.structural_integrity;
    inputs.performance_degradation = 1.0 - degradation_system_.performance_degradation;
    
    return inputs;
}

// Simulation control
void DysonSphere::start_simulation(double time_step) {
    if (simulation_running_) return;
    
    simulation_running_ = true;
    
    // Start simulation thread
    std::thread([this, time_step]() {
        while (simulation_running_) {
            auto start_time = std::chrono::high_resolution_clock::now();
            
            update_physics(time_step);
            
            // Maintain real-time or faster execution
            auto end_time = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            auto target_duration = std::chrono::milliseconds(static_cast<int>(time_step * 1000));
            
            if (elapsed < target_duration) {
                std::this_thread::sleep_for(target_duration - elapsed);
            }
        }
    }).detach();
}

void DysonSphere::stop_simulation() {
    simulation_running_ = false;
}

// Status and diagnostic methods
double DysonSphere::get_current_distance_from_sun() const {
    const auto& pos = orbital_system_.position_heliocentric;
    return std::sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
}

double DysonSphere::get_current_solar_flux() const {
    return solar_calc_.radiation_flux_at_distance(get_current_distance_from_sun());
}

double DysonSphere::get_energy_efficiency() const {
    if (energy_system_.power_consumption_rate <= 0) return 0.0;
    return energy_system_.power_generation_rate / energy_system_.power_consumption_rate;
}

double DysonSphere::get_mass_efficiency() const {
    return energy_system_.power_generation_rate / sphere_mass_;
}

double DysonSphere::get_structural_integrity() const {
    return degradation_system_.structural_integrity;
}

double DysonSphere::get_thermal_stability() const {
    if (thermal_system_.temperature_current > thermal_system_.temperature_critical) return 0.0;
    return 1.0 - (thermal_system_.temperature_current / thermal_system_.temperature_critical);
}

double DysonSphere::get_orbital_stability() const {
    double distance = get_current_distance_from_sun();
    double perturbation_ratio = std::sqrt(orbital_system_.acceleration_total[0] * orbital_system_.acceleration_total[0] +
                                        orbital_system_.acceleration_total[1] * orbital_system_.acceleration_total[1] +
                                        orbital_system_.acceleration_total[2] * orbital_system_.acceleration_total[2]) /
                               (G * SolarParams::mass / (distance * distance));
    
    return std::max(0.0, 1.0 - perturbation_ratio * 1000.0); // Scale factor for sensitivity
}

double DysonSphere::get_overall_health() const {
    double weights[4] = {0.3, 0.25, 0.25, 0.2}; // Energy, thermal, structural, orbital
    return weights[0] * std::min(1.0, get_energy_efficiency()) +
           weights[1] * get_thermal_stability() +
           weights[2] * get_structural_integrity() +
           weights[3] * get_orbital_stability();
}

bool DysonSphere::is_structurally_stable() const {
    return degradation_system_.structural_integrity > 0.7 && 
           !std::any_of(structure_.begin(), structure_.end(), 
                       [](const StructuralElement& e) { return e.is_failed; });
}

bool DysonSphere::is_thermally_stable() const {
    return thermal_system_.temperature_current < thermal_system_.temperature_critical &&
           thermal_system_.temperature_current > 200.0; // Above minimum operating temp
}

bool DysonSphere::is_orbitally_stable() const {
    double distance = get_current_distance_from_sun();
    return distance > MIN_ORBITAL_DISTANCE && distance < MAX_ORBITAL_DISTANCE &&
           get_orbital_stability() > 0.8;
}

bool DysonSphere::requires_emergency_action() const {
    return !is_structurally_stable() || !is_thermally_stable() || !is_orbitally_stable() ||
           energy_system_.energy_stored < 0.1 * energy_system_.energy_capacity_max;
}

double DysonSphere::get_estimated_lifetime() const {
    return degradation_system_.estimated_lifetime / year;
}

void DysonSphere::print_status_report() {
    std::lock_guard<std::mutex> lock(state_mutex_);
    
    std::cout << "\n=== DYSON SPHERE STATUS REPORT ===" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    
    // Basic properties
    std::cout << "Radius: " << sphere_radius_ / 1e6 << " km" << std::endl;
    std::cout << "Mass: " << sphere_mass_ / 1e15 << " Pg" << std::endl;
    std::cout << "Material: " << material_.name << std::endl;
    std::cout << "Simulation Time: " << simulation_time_ / year << " years" << std::endl;
    
    // Orbital status
    std::cout << "\n--- ORBITAL STATUS ---" << std::endl;
    std::cout << "Distance from Sun: " << get_current_distance_from_sun() / AU << " AU" << std::endl;
    std::cout << "Orbital Velocity: " << orbital_system_.orbital_velocity / 1000.0 << " km/s" << std::endl;
    std::cout << "Orbital Stability: " << get_orbital_stability() * 100.0 << "%" << std::endl;
    
    // Energy status
    std::cout << "\n--- ENERGY STATUS ---" << std::endl;
    std::cout << "Energy Stored: " << energy_system_.energy_stored / 1e24 << " YJ" << std::endl;
    std::cout << "Energy Level: " << (energy_system_.energy_stored / energy_system_.energy_capacity_max) * 100.0 << "%" << std::endl;
    std::cout << "Power Generation: " << energy_system_.power_generation_rate / 1e12 << " TW" << std::endl;
    std::cout << "Power Consumption: " << energy_system_.power_consumption_rate / 1e12 << " TW" << std::endl;
    std::cout << "Energy Efficiency: " << get_energy_efficiency() * 100.0 << "%" << std::endl;
    
    // Thermal status
    std::cout << "\n--- THERMAL STATUS ---" << std::endl;
    std::cout << "Current Temperature: " << thermal_system_.temperature_current << " K" << std::endl;
    std::cout << "Critical Temperature: " << thermal_system_.temperature_critical << " K" << std::endl;
    std::cout << "Thermal Stability: " << get_thermal_stability() * 100.0 << "%" << std::endl;
    
    // Structural status
    std::cout << "\n--- STRUCTURAL STATUS ---" << std::endl;
    std::cout << "Structural Integrity: " << degradation_system_.structural_integrity * 100.0 << "%" << std::endl;
    std::cout << "Failed Elements: " << std::count_if(structure_.begin(), structure_.end(),
                                                     [](const StructuralElement& e) { return e.is_failed; }) << std::endl;
    std::cout << "Performance Degradation: " << (1.0 - degradation_system_.performance_degradation) * 100.0 << "%" << std::endl;
    std::cout << "Estimated Lifetime: " << get_estimated_lifetime() << " years" << std::endl;
    
    // Overall health
    std::cout << "\n--- OVERALL HEALTH ---" << std::endl;
    std::cout << "Overall Health Score: " << get_overall_health() * 100.0 << "%" << std::endl;
    std::cout << "Structurally Stable: " << (is_structurally_stable() ? "YES" : "NO") << std::endl;
    std::cout << "Thermally Stable: " << (is_thermally_stable() ? "YES" : "NO") << std::endl;
    std::cout << "Orbitally Stable: " << (is_orbitally_stable() ? "YES" : "NO") << std::endl;
    std::cout << "Emergency Action Required: " << (requires_emergency_action() ? "YES" : "NO") << std::endl;
    
    std::cout << "==================================\n" << std::endl;
}

// Static utility methods
MaterialProperties DysonSphere::get_material_properties(const std::string& material_name) {
    MaterialProperties props;
    props.name = material_name;
    
    if (material_name == "aluminum") {
        props.density = MaterialConstants::rho_aluminum;
        props.youngs_modulus = MaterialConstants::E_aluminum;
        props.tensile_strength = 200e6;
        props.melting_point = MaterialConstants::T_melt_aluminum;
        props.thermal_conductivity = 237.0;
        props.specific_heat = 897.0;
        props.thermal_expansion = MaterialConstants::alpha_aluminum;
        props.emissivity = 0.05;
        props.absorptivity = 0.15;
        props.radiation_resistance = 1e5;
    } else if (material_name == "silicon") {
        props.density = MaterialConstants::rho_silicon;
        props.youngs_modulus = MaterialConstants::E_silicon;
        props.tensile_strength = 150e6;
        props.melting_point = MaterialConstants::T_melt_silicon;
        props.thermal_conductivity = 148.0;
        props.specific_heat = 700.0;
        props.thermal_expansion = MaterialConstants::alpha_silicon;
        props.emissivity = 0.65;
        props.absorptivity = 0.7;
        props.radiation_resistance = 1e3;
    } else if (material_name == "iron") {
        props.density = MaterialConstants::rho_iron;
        props.youngs_modulus = MaterialConstants::E_iron;
        props.tensile_strength = 400e6;
        props.melting_point = MaterialConstants::T_melt_iron;
        props.thermal_conductivity = 80.0;
        props.specific_heat = 449.0;
        props.thermal_expansion = MaterialConstants::alpha_iron;
        props.emissivity = 0.28;
        props.absorptivity = 0.4;
        props.radiation_resistance = 5e4;
    } else if (material_name == "carbon") {
        props.density = MaterialConstants::rho_carbon;
        props.youngs_modulus = MaterialConstants::E_carbon;
        props.tensile_strength = 2e9;
        props.melting_point = MaterialConstants::T_melt_carbon;
        props.thermal_conductivity = 1000.0;
        props.specific_heat = 709.0;
        props.thermal_expansion = MaterialConstants::alpha_carbon;
        props.emissivity = 0.81;
        props.absorptivity = 0.85;
        props.radiation_resistance = 1e6;
    } else if (material_name == "tungsten") {
        props.density = MaterialConstants::rho_tungsten;
        props.youngs_modulus = MaterialConstants::E_tungsten;
        props.tensile_strength = 500e6;
        props.melting_point = MaterialConstants::T_melt_tungsten;
        props.thermal_conductivity = 173.0;
        props.specific_heat = 132.0;
        props.thermal_expansion = MaterialConstants::alpha_tungsten;
        props.emissivity = 0.35;
        props.absorptivity = 0.4;
        props.radiation_resistance = 1e7;
    } else {
        // Default to aluminum
        return get_material_properties("aluminum");
    }
    
    return props;
}

double DysonSphere::calculate_sphere_mass(double radius, double thickness, double density) {
    double outer_volume = (4.0 / 3.0) * M_PI * std::pow(radius, 3);
    double inner_volume = (4.0 / 3.0) * M_PI * std::pow(radius - thickness, 3);
    return (outer_volume - inner_volume) * density;
}

double DysonSphere::calculate_sphere_surface_area(double radius) {
    return 4.0 * M_PI * radius * radius;
}

std::vector<std::string> DysonSphere::get_available_materials() {
    return {"aluminum", "silicon", "iron", "carbon", "tungsten"};
}

// Utility namespace implementation
namespace DysonSphereUtils {
    double calculate_orbital_period(double semi_major_axis, double central_mass) {
        return 2.0 * M_PI * std::sqrt(std::pow(semi_major_axis, 3) / (G * central_mass));
    }
    
    double calculate_escape_velocity(double mass, double radius) {
        return std::sqrt(2.0 * G * mass / radius);
    }
    
    std::array<double, 3> spherical_to_cartesian(double r, double theta, double phi) {
        return {r * std::sin(theta) * std::cos(phi),
                r * std::sin(theta) * std::sin(phi),
                r * std::cos(theta)};
    }
    
    std::array<double, 3> cartesian_to_spherical(const std::array<double, 3>& pos) {
        double r = std::sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        double theta = std::acos(pos[2] / r);
        double phi = std::atan2(pos[1], pos[0]);
        return {r, theta, phi};
    }
}