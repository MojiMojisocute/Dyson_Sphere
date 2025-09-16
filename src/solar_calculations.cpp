#include "solar_calculations.h"

using namespace PhysicsConstants;
using namespace SolarParams;

SolarCalculations::SolarCalculations() : rng_(std::random_device{}()) {}

// Basic solar properties
double SolarCalculations::surface_gravity() const {
    return G * mass / (radius * radius);
}

double SolarCalculations::escape_velocity() const {
    return std::sqrt(2.0 * G * mass / radius);
}

double SolarCalculations::schwarzschild_radius() const {
    return 2.0 * G * mass / (c * c);
}

double SolarCalculations::mass_loss_rate_fusion() const {
    return luminosity / (c * c);
}

// Internal spectral calculation helpers
double SolarCalculations::planck_function(double wavelength, double temperature) const {
    const double hc_lambda = h * c / wavelength;
    const double exp_term = std::exp(hc_lambda / (k_B * temperature)) - 1.0;
    return (2.0 * h * c * c) / (std::pow(wavelength, 5) * exp_term);
}

double SolarCalculations::wien_displacement_law(double temperature) const {
    return 2.8977729e-3 / temperature; // Wien's displacement constant
}

double SolarCalculations::rayleigh_jeans_approx(double wavelength, double temperature) const {
    return (2.0 * c * k_B * temperature) / std::pow(wavelength, 4);
}

// Improved spectral calculations
SpectralData SolarCalculations::calculate_planckian_spectrum(int num_bins, double lambda_min, double lambda_max) const {
    SpectralData spectrum;
    spectrum.wavelengths.resize(num_bins);
    spectrum.intensities.resize(num_bins);
    spectrum.photon_flux.resize(num_bins);
    
    const double log_lambda_min = std::log(lambda_min);
    const double log_lambda_max = std::log(lambda_max);
    const double d_log_lambda = (log_lambda_max - log_lambda_min) / (num_bins - 1);
    
    for (int i = 0; i < num_bins; ++i) {
        const double log_lambda = log_lambda_min + i * d_log_lambda;
        const double wavelength = std::exp(log_lambda);
        
        spectrum.wavelengths[i] = wavelength;
        spectrum.intensities[i] = planck_function(wavelength, temperature);
        
        // Calculate photon flux: N_photons = Intensity * wavelength / (h * c)
        spectrum.photon_flux[i] = spectrum.intensities[i] * wavelength / (h * c);
    }
    
    return spectrum;
}

std::array<double, 5> SolarCalculations::detailed_spectral_bands() const {
    // More accurate spectral distribution based on solar spectrum observations
    // UV-C (100-280nm), UV-B (280-315nm), UV-A (315-400nm), Visible (400-700nm), IR (700nm+)
    std::array<double, 5> bands;
    
    // Calculate integrated intensities for each band
    const auto spectrum = calculate_planckian_spectrum(10000, 1e-7, 1e-4);
    
    std::array<std::pair<double, double>, 5> band_limits = {{
        {100e-9, 280e-9},   // UV-C
        {280e-9, 315e-9},   // UV-B
        {315e-9, 400e-9},   // UV-A
        {400e-9, 700e-9},   // Visible
        {700e-9, 1e-4}      // IR
    }};
    
    double total_intensity = 0.0;
    for (const auto& intensity : spectrum.intensities) {
        total_intensity += intensity;
    }
    
    for (size_t band = 0; band < 5; ++band) {
        double band_intensity = 0.0;
        for (size_t i = 0; i < spectrum.wavelengths.size(); ++i) {
            if (spectrum.wavelengths[i] >= band_limits[band].first && 
                spectrum.wavelengths[i] <= band_limits[band].second) {
                band_intensity += spectrum.intensities[i];
            }
        }
        bands[band] = band_intensity / total_intensity;
    }
    
    return bands;
}

double SolarCalculations::spectral_irradiance_at_distance(double distance, double wavelength) const {
    const double intensity_at_sun = planck_function(wavelength, temperature);
    const double solid_angle = M_PI * radius * radius / (distance * distance);
    return intensity_at_sun * solid_angle;
}

// Enhanced radiation calculations
double SolarCalculations::radiation_flux_at_distance(double distance) const {
    return luminosity / (4.0 * M_PI * distance * distance);
}

double SolarCalculations::radiation_pressure(double distance, double absorption_coeff, double reflection_coeff) const {
    const double flux = radiation_flux_at_distance(distance);
    return flux * (absorption_coeff + reflection_coeff) / c;
}

double SolarCalculations::photon_momentum_transfer(double distance, double cross_section) const {
    const double flux = radiation_flux_at_distance(distance);
    return flux * cross_section / c;
}

// Plasma physics helpers
double SolarCalculations::debye_length(double density, double temperature) const {
    return std::sqrt(epsilon_0 * k_B * temperature / (density * e * e));
}

double SolarCalculations::plasma_frequency(double density) const {
    return std::sqrt(density * e * e / (epsilon_0 * m_e));
}

double SolarCalculations::cyclotron_frequency(double magnetic_field, double mass, double charge) const {
    return charge * magnetic_field / mass;
}

// MHD calculations
double SolarCalculations::alfven_velocity(double magnetic_field, double density) const {
    const double mass_density = density * m_p; // Assuming proton plasma
    return magnetic_field / std::sqrt(mu_0 * mass_density);
}

double SolarCalculations::magnetic_pressure(double magnetic_field) const {
    return magnetic_field * magnetic_field / (2.0 * mu_0);
}

double SolarCalculations::magnetic_reynolds_number(double velocity, double length_scale) const {
    const double magnetic_diffusivity = 1.0 / (mu_0 * 1e-4); // Typical solar plasma conductivity
    return velocity * length_scale / magnetic_diffusivity;
}

// Advanced solar wind modeling
SolarWindData SolarCalculations::calculate_solar_wind(double distance, double solar_activity) const {
    SolarWindData wind_data;
    
    // Parker wind model with solar activity modulation
    const double base_density = 5e6; // particles/m³ at 1 AU
    const double base_velocity = 400e3; // m/s
    const double base_temperature = 1e5; // K
    
    // Distance scaling
    wind_data.density = base_density * solar_activity * AU * AU / (distance * distance);
    wind_data.velocity = base_velocity * (1.0 + 0.5 * solar_activity);
    wind_data.temperature = base_temperature * std::pow(distance / AU, -0.7);
    
    // Magnetic field (Parker spiral)
    const double omega_sun = 2.0 * M_PI / rotation_period;
    const double spiral_angle = std::atan(omega_sun * distance / wind_data.velocity);
    wind_data.magnetic_field = 5e-9 * (AU / distance) * (AU / distance) * solar_activity;
    
    // Dynamic pressure
    wind_data.dynamic_pressure = wind_data.density * m_p * wind_data.velocity * wind_data.velocity;
    
    // Particle distributions
    PlasmaParticle proton, alpha, electron;
    
    // Protons (90% of particles)
    proton.mass = m_p;
    proton.charge = e;
    proton.energy = 1.5 * k_B * wind_data.temperature;
    proton.velocity = std::sqrt(2.0 * proton.energy / proton.mass);
    proton.cross_section = 1e-19; // Typical proton cross-section
    
    // Alpha particles (8% of particles)
    alpha.mass = 4.0 * m_p;
    alpha.charge = 2.0 * e;
    alpha.energy = 1.5 * k_B * wind_data.temperature;
    alpha.velocity = std::sqrt(2.0 * alpha.energy / alpha.mass);
    alpha.cross_section = 4e-19;
    
    // Electrons (balance charge)
    electron.mass = m_e;
    electron.charge = -e;
    electron.energy = 1.5 * k_B * wind_data.temperature;
    electron.velocity = std::sqrt(2.0 * electron.energy / electron.mass);
    electron.cross_section = 1e-24;
    
    wind_data.particles = {proton, alpha, electron};
    
    return wind_data;
}

double SolarCalculations::parker_wind_velocity(double distance) const {
    const double critical_radius = G * mass / (2.0 * k_B * corona_temp / m_p);
    if (distance < critical_radius) {
        return std::sqrt(k_B * corona_temp / m_p);
    } else {
        return std::sqrt(2.0 * k_B * corona_temp / m_p) * std::sqrt(1.0 - critical_radius / distance);
    }
}

double SolarCalculations::solar_wind_dynamic_pressure(double distance) const {
    const auto wind_data = calculate_solar_wind(distance);
    return wind_data.dynamic_pressure;
}

// Improved sputtering calculations
double SolarCalculations::sputtering_yield(const PlasmaParticle& particle, const std::string& material) const {
    // Realistic sputtering yield based on particle energy and material properties
    std::map<std::string, std::array<double, 3>> material_props = {
        {"aluminum", {25.0, 2.7, 0.15}},  // threshold_eV, surface_binding_eV, yield_factor
        {"silicon", {30.0, 4.7, 0.12}},
        {"iron", {35.0, 4.3, 0.10}},
        {"carbon", {40.0, 7.4, 0.08}},
        {"tungsten", {45.0, 8.9, 0.05}}
    };
    
    auto it = material_props.find(material);
    if (it == material_props.end()) {
        return 0.01; // Default low yield
    }
    
    const auto& props = it->second;
    const double threshold_energy = props[0] * e; // Convert eV to J
    const double binding_energy = props[1] * e;
    const double yield_factor = props[2];
    
    if (particle.energy < threshold_energy) {
        return 0.0;
    }
    
    // Sigmund theory for sputtering yield
    const double reduced_energy = particle.energy / binding_energy;
    const double yield = yield_factor * reduced_energy * std::exp(-std::sqrt(reduced_energy));
    
    return std::min(yield, 10.0); // Cap at 10 atoms/ion
}

double SolarCalculations::plasma_erosion_rate(double distance, const std::string& material, double surface_area) const {
    const auto wind_data = calculate_solar_wind(distance);
    double total_erosion = 0.0;
    
    for (const auto& particle : wind_data.particles) {
        const double flux = wind_data.density * particle.velocity / 3.0; // Assume equal distribution
        const double yield = sputtering_yield(particle, material);
        total_erosion += flux * yield * surface_area;
    }
    
    return total_erosion; // atoms/s
}

std::vector<double> SolarCalculations::energy_dependent_cross_sections(double energy, const std::string& interaction) const {
    std::vector<double> cross_sections;
    
    if (interaction == "photoionization") {
        // Photoionization cross-sections for different photon energies
        const double threshold = 13.6 * e; // Hydrogen ionization
        if (energy < threshold) {
            cross_sections.push_back(0.0);
        } else {
            const double sigma_0 = 6.3e-22; // m²
            cross_sections.push_back(sigma_0 * std::pow(threshold / energy, 3.5));
        }
    } else if (interaction == "compton") {
        // Compton scattering cross-section (Klein-Nishina)
        const double alpha = energy / (m_e * c * c);
        const double sigma_thomson = 6.652e-29; // m²
        const double term1 = 1.0 + alpha;
        const double term2 = (2.0 * alpha * (1.0 + alpha)) / (1.0 + 2.0 * alpha);
        const double term3 = std::log(1.0 + 2.0 * alpha) / alpha;
        const double term4 = (1.0 + 3.0 * alpha) / ((1.0 + 2.0 * alpha) * (1.0 + 2.0 * alpha));
        
        cross_sections.push_back(sigma_thomson * (term2 - term3 + term4) / term1);
    }
    
    return cross_sections;
}

// Magnetic field modeling with solar cycle
double SolarCalculations::magnetic_field_strength(double distance, double time_years) const {
    const double cycle_modulation = solar_cycle_modulation(time_years);
    const double base_field = 5e-9 * cycle_modulation; // Tesla at 1 AU
    return base_field * AU * AU / (distance * distance);
}

double SolarCalculations::solar_cycle_modulation(double time_years) const {
    const double cycle_period = 11.0; // years
    const double phase = 2.0 * M_PI * time_years / cycle_period;
    return 0.7 + 0.3 * std::sin(phase); // Varies between 0.4 and 1.0
}

std::array<double, 3> SolarCalculations::magnetic_field_vector(double distance, double latitude, double longitude, double time_years) const {
    const double field_strength = magnetic_field_strength(distance, time_years);
    const double omega_sun = 2.0 * M_PI / rotation_period;
    const double wind_velocity = parker_wind_velocity(distance);
    
    // Parker spiral magnetic field
    const double spiral_angle = std::atan(omega_sun * distance * std::cos(latitude) / wind_velocity);
    
    std::array<double, 3> B_field;
    B_field[0] = field_strength * std::cos(spiral_angle) * std::cos(latitude); // Radial
    B_field[1] = -field_strength * std::sin(spiral_angle) * std::cos(latitude); // Azimuthal
    B_field[2] = field_strength * std::sin(latitude); // Polar
    
    return B_field;
}

// Enhanced flare modeling
SolarFlareData SolarCalculations::generate_solar_flare() const {
    SolarFlareData flare;
    
    // Power law distribution with exponential cutoff for flare energy
    std::exponential_distribution<double> exp_dist(1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    const double power_index = -1.8; // Observed power law index
    const double energy_cutoff = 1e26; // J, maximum observed flare energy
    const double energy_min = 1e20; // J, minimum detectable flare
    
    double energy;
    do {
        const double u = uniform(rng_);
        energy = energy_min * std::pow(u, 1.0 / (power_index + 1.0));
    } while (energy > energy_cutoff);
    
    flare.energy = energy;
    
    // Empirical relationships for flare properties
    flare.duration = 100.0 * std::pow(energy / 1e22, 0.3); // seconds
    flare.peak_flux = energy / (4.0 * M_PI * AU * AU * flare.duration); // W/m²
    flare.proton_flux = std::pow(energy / 1e20, 1.5) * 1e4; // particles/cm²/s
    
    return flare;
}

double SolarCalculations::flare_occurrence_rate(double solar_activity) const {
    // More realistic flare occurrence rates based on solar observations
    const double base_rate = 0.01; // flares per hour during solar minimum
    return base_rate * solar_activity * solar_activity; // Quadratic dependence on activity
}

double SolarCalculations::flare_energy_distribution() const {
    // Generate flare energy following observed power law distribution
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const double u = uniform(rng_);
    const double power_index = -1.8;
    const double energy_min = 1e20; // J
    const double energy_max = 1e26; // J
    
    return energy_min * std::pow(energy_max / energy_min, u);
}

// Improved CME modeling
CMEData SolarCalculations::calculate_cme_properties(double distance) const {
    CMEData cme;
    
    // CME occurrence probability (more realistic)
    cme.probability = 0.5 / (24.0 * 7.0); // ~0.5 CMEs per week
    
    // CME velocity distribution (typical range 200-2000 km/s)
    std::lognormal_distribution<double> velocity_dist(6.2, 0.5);
    cme.velocity = std::min(velocity_dist(rng_), 3000e3); // Cap at 3000 km/s
    
    // CME mass (typical range 1e12 - 1e16 kg)
    std::lognormal_distribution<double> mass_dist(32.0, 2.0);
    cme.mass = mass_dist(rng_);
    
    // Embedded magnetic field
    cme.magnetic_field = 1e-8 * std::pow(AU / distance, 2.0); // Tesla
    
    // Event duration at given distance
    const double cme_width = 0.2 * AU; // Typical CME width
    cme.duration = cme_width / cme.velocity;
    
    return cme;
}

bool SolarCalculations::cme_occurrence_check(double time_hours, double solar_activity) const {
    const double base_rate = 0.5 / (24.0 * 7.0); // CMEs per hour
    const double rate = base_rate * solar_activity;
    
    std::poisson_distribution<int> poisson(rate * time_hours);
    return poisson(rng_) > 0;
}

double SolarCalculations::cme_travel_time(double distance) const {
    const double typical_velocity = 500e3; // m/s, typical CME velocity
    return distance / typical_velocity;
}

// Thermal calculations
double SolarCalculations::equilibrium_temperature(double distance, double absorptivity, double emissivity) const {
    const double flux = radiation_flux_at_distance(distance);
    const double absorbed_power = absorptivity * flux;
    const double T_eq_4 = absorbed_power / (emissivity * sigma);
    return std::pow(T_eq_4, 0.25);
}

double SolarCalculations::greenhouse_effect_temperature(double distance, double atmosphere_opacity) const {
    const double T_eq = equilibrium_temperature(distance);
    const double greenhouse_factor = 1.0 + 0.75 * atmosphere_opacity;
    return T_eq * std::pow(greenhouse_factor, 0.25);
}

// Material properties and interactions
double SolarCalculations::thermal_expansion_stress(double temperature_change, const std::string& material) const {
    // Thermal expansion coefficients (1/K) and Young's moduli (Pa)
    std::map<std::string, std::pair<double, double>> material_props = {
        {"aluminum", {23.1e-6, 70e9}},
        {"silicon", {2.6e-6, 130e9}},
        {"iron", {11.8e-6, 200e9}},
        {"carbon", {7.1e-6, 500e9}},
        {"tungsten", {4.5e-6, 411e9}}
    };
    
    auto it = material_props.find(material);
    if (it == material_props.end()) {
        return 0.0; // Unknown material
    }
    
    const double alpha = it->second.first;   // Thermal expansion coefficient
    const double E = it->second.second;      // Young's modulus
    
    return alpha * E * temperature_change;   // Thermal stress (Pa)
}

double SolarCalculations::radiation_damage_rate(double distance, const std::string& material) const {
    const double flux = radiation_flux_at_distance(distance);
    const double photon_energy = h * c / (500e-9); // Assuming visible light
    const double photon_flux = flux / photon_energy;
    
    // Material-dependent damage cross-sections (m²)
    std::map<std::string, double> damage_cross_sections = {
        {"silicon", 1e-21},    // Semiconductor materials more sensitive
        {"aluminum", 1e-23},
        {"iron", 5e-24},
        {"carbon", 1e-24},     // Carbon materials very resistant
        {"tungsten", 1e-25}
    };
    
    auto it = damage_cross_sections.find(material);
    const double cross_section = (it != damage_cross_sections.end()) ? it->second : 1e-23;
    
    return photon_flux * cross_section; // Damage events per m² per second
}

double SolarCalculations::atomic_displacement_rate(double distance, const std::string& material) const {
    const auto wind_data = calculate_solar_wind(distance);
    double displacement_rate = 0.0;
    
    for (const auto& particle : wind_data.particles) {
        if (particle.charge > 0) { // Only charged particles cause displacements
            const double flux = wind_data.density * particle.velocity / 3.0;
            const double threshold_energy = 25.0 * e; // Typical displacement threshold
            
            if (particle.energy > threshold_energy) {
                const double displacement_cross_section = 1e-24; // m², typical value
                displacement_rate += flux * displacement_cross_section;
            }
        }
    }
    
    return displacement_rate; // Displacements per m² per second
}

// Advanced dynamics
double SolarCalculations::orbital_perturbation_acceleration(double distance, double mass_object) const {
    // Solar radiation pressure perturbation
    const double radiation_pressure = this->radiation_pressure(distance);
    const double radiation_acceleration = radiation_pressure / (mass_object / (4.0 * M_PI));
    
    // Solar wind perturbation
    const double wind_pressure = solar_wind_dynamic_pressure(distance);
    const double wind_acceleration = wind_pressure / (mass_object / (4.0 * M_PI));
    
    return radiation_acceleration + wind_acceleration;
}

std::array<double, 3> SolarCalculations::poynting_robertson_drag(double distance, double particle_radius, 
                                                               double particle_density) const {
    const double particle_mass = (4.0/3.0) * M_PI * std::pow(particle_radius, 3) * particle_density;
    const double cross_section = M_PI * particle_radius * particle_radius;
    const double flux = radiation_flux_at_distance(distance);
    
    // Poynting-Robertson drag coefficient
    const double beta = cross_section * flux / (4.0 * M_PI * c * particle_mass);
    const double orbital_velocity = std::sqrt(G * mass / distance);
    
    std::array<double, 3> drag_acceleration;
    drag_acceleration[0] = -beta * orbital_velocity;  // Radial (inward)
    drag_acceleration[1] = -beta * orbital_velocity;  // Tangential (retrograde)
    drag_acceleration[2] = 0.0;                       // Normal (negligible)
    
    return drag_acceleration;
}

// Environmental hazard assessment
double SolarCalculations::total_radiation_dose_rate(double distance, double shielding_thickness) const {
    const double flux = radiation_flux_at_distance(distance);
    const auto wind_data = calculate_solar_wind(distance);
    
    // Solar electromagnetic radiation dose
    const double em_dose_rate = flux * 1e-6; // Convert to Gy/s (simplified)
    
    // Solar energetic particle dose
    double particle_dose_rate = 0.0;
    for (const auto& particle : wind_data.particles) {
        if (particle.charge != 0) {
            const double flux_particles = wind_data.density * particle.velocity / 3.0;
            const double energy_deposit = particle.energy * 1e-15; // J to MeV conversion factor
            particle_dose_rate += flux_particles * energy_deposit;
        }
    }
    
    // Shielding attenuation (exponential approximation)
    const double attenuation_length = 10.0; // g/cm² typical for aluminum
    const double shielding_factor = std::exp(-shielding_thickness / attenuation_length);
    
    return (em_dose_rate + particle_dose_rate) * shielding_factor;
}

double SolarCalculations::micrometeorite_impact_rate(double distance, double surface_area) const {
    // Zodiacal dust distribution model
    const double base_flux = 1e-12; // kg/m²/s at 1 AU
    const double distance_scaling = std::pow(AU / distance, 1.3); // Empirical scaling
    const double flux_density = base_flux * distance_scaling;
    
    return flux_density * surface_area; // kg/s total impact mass
}

// Utility functions
void SolarCalculations::set_random_seed(unsigned int seed) {
    rng_.seed(seed);
}

double SolarCalculations::interpolate_solar_activity(double time_years) const {
    // Solar cycle model with irregular variations
    const double cycle_period = 11.0;
    const double base_cycle = solar_cycle_modulation(time_years);
    
    // Add some randomness for realistic variations
    std::normal_distribution<double> noise(0.0, 0.1);
    const double variation = noise(rng_);
    
    return std::max(0.1, std::min(2.0, base_cycle + variation));
}

std::vector<double> SolarCalculations::time_series_prediction(double start_time, double end_time, 
                                                            double dt, const std::string& parameter) const {
    std::vector<double> time_series;
    const int num_points = static_cast<int>((end_time - start_time) / dt) + 1;
    time_series.reserve(num_points);
    
    for (int i = 0; i < num_points; ++i) {
        const double time = start_time + i * dt;
        double value = 0.0;
        
        if (parameter == "solar_activity") {
            value = interpolate_solar_activity(time);
        } else if (parameter == "magnetic_field") {
            value = magnetic_field_strength(AU, time);
        } else if (parameter == "solar_wind_velocity") {
            const auto wind_data = calculate_solar_wind(AU, interpolate_solar_activity(time));
            value = wind_data.velocity;
        } else if (parameter == "flare_rate") {
            value = flare_occurrence_rate(interpolate_solar_activity(time));
        } else if (parameter == "cme_rate") {
            const double activity = interpolate_solar_activity(time);
            value = 0.5 / (24.0 * 7.0) * activity; // CMEs per hour
        }
        
        time_series.push_back(value);
    }
    
    return time_series;
}