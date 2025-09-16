#ifndef SOLAR_CALCULATIONS_H
#define SOLAR_CALCULATIONS_H

#include "physicsConstants.h"

// Solar Parameters
namespace SolarParams {
    constexpr double mass = 1.98847e30;      // Solar mass (kg)
    constexpr double radius = 6.957e8;       // Solar radius (m)
    constexpr double temperature = 5778.0;   // Solar effective temperature (K)
    constexpr double luminosity = 3.828e26;  // Solar luminosity (W)
    constexpr double rotation_period = 2.14e6; // Solar rotation period (s) ~25 days
    constexpr double magnetic_dipole = 8.0e22;  // Solar magnetic dipole moment (A·m²)
    constexpr double corona_temp = 2e6;      // Coronal temperature (K)
}

// Spectral data structure
struct SpectralData {
    std::vector<double> wavelengths;  // Wavelength bins (m)
    std::vector<double> intensities;  // Spectral intensity (W/m²/m)
    std::vector<double> photon_flux;  // Photon flux (photons/m²/s/m)
};

// Plasma particle data
struct PlasmaParticle {
    double energy;      // Kinetic energy (J)
    double mass;        // Particle mass (kg)
    double charge;      // Particle charge (C)
    double velocity;    // Particle velocity (m/s)
    double cross_section; // Interaction cross-section (m²)
};

// Solar wind data structure
struct SolarWindData {
    double density;     // Number density (particles/m³)
    double velocity;    // Bulk velocity (m/s)
    double temperature; // Temperature (K)
    double magnetic_field; // Magnetic field strength (T)
    double dynamic_pressure; // Dynamic pressure (Pa)
    std::vector<PlasmaParticle> particles; // Particle distributions
};

// CME data structure
struct CMEData {
    double probability;  // Occurrence probability
    double velocity;     // CME velocity (m/s)
    double mass;         // CME mass (kg)
    double magnetic_field; // Embedded magnetic field (T)
    double duration;     // Event duration (s)
};

// Solar flare data structure
struct SolarFlareData {
    double energy;       // Flare energy (J)
    double duration;     // Flare duration (s)
    double peak_flux;    // Peak X-ray flux (W/m²)
    double proton_flux;  // Energetic proton flux (particles/cm²/s)
};

class SolarCalculations {
private:
    // Random number generator
    mutable std::mt19937 rng_;
    
    // Internal calculation helpers
    double planck_function(double wavelength, double temperature) const;
    double wien_displacement_law(double temperature) const;
    double rayleigh_jeans_approx(double wavelength, double temperature) const;
    
    // Plasma physics helpers
    double debye_length(double density, double temperature) const;
    double plasma_frequency(double density) const;
    double cyclotron_frequency(double magnetic_field, double mass, double charge) const;
    
    // MHD calculations
    double alfven_velocity(double magnetic_field, double density) const;
    double magnetic_pressure(double magnetic_field) const;
    double magnetic_reynolds_number(double velocity, double length_scale) const;

public:
    SolarCalculations();
    
    // Basic solar properties
    double surface_gravity() const;
    double escape_velocity() const;
    double schwarzschild_radius() const;
    double mass_loss_rate_fusion() const;
    
    // Improved spectral calculations
    SpectralData calculate_planckian_spectrum(int num_bins = 1000, 
                                            double lambda_min = 1e-8, 
                                            double lambda_max = 1e-4) const;
    
    std::array<double, 5> detailed_spectral_bands() const; // UV-C, UV-B, UV-A, Visible, IR
    double spectral_irradiance_at_distance(double distance, double wavelength) const;
    
    // Enhanced radiation calculations
    double radiation_flux_at_distance(double distance) const;
    double radiation_pressure(double distance, double absorption_coeff = 1.0, 
                            double reflection_coeff = 0.0) const;
    double photon_momentum_transfer(double distance, double cross_section) const;
    
    // Advanced solar wind modeling
    SolarWindData calculate_solar_wind(double distance, double solar_activity = 1.0) const;
    double parker_wind_velocity(double distance) const;
    double solar_wind_dynamic_pressure(double distance) const;
    
    // Plasma interaction improvements
    double sputtering_yield(const PlasmaParticle& particle, const std::string& material) const;
    double plasma_erosion_rate(double distance, const std::string& material, 
                              double surface_area) const;
    std::vector<double> energy_dependent_cross_sections(double energy, 
                                                       const std::string& interaction) const;
    
    // Magnetic field modeling with solar cycle
    double magnetic_field_strength(double distance, double time_years = 0.0) const;
    double solar_cycle_modulation(double time_years) const;
    std::array<double, 3> magnetic_field_vector(double distance, double latitude, 
                                               double longitude, double time_years = 0.0) const;
    
    // Enhanced flare modeling
    SolarFlareData generate_solar_flare() const;
    double flare_occurrence_rate(double solar_activity = 1.0) const;
    double flare_energy_distribution() const; // Using power law + exponential cutoff
    
    // Improved CME modeling
    CMEData calculate_cme_properties(double distance) const;
    bool cme_occurrence_check(double time_hours, double solar_activity = 1.0) const;
    double cme_travel_time(double distance) const;
    
    // Thermal equilibrium calculations
    double equilibrium_temperature(double distance, double absorptivity = 0.7, 
                                 double emissivity = 0.9) const;
    double greenhouse_effect_temperature(double distance, double atmosphere_opacity = 0.0) const;
    
    // Material interaction models
    double thermal_expansion_stress(double temperature_change, const std::string& material) const;
    double radiation_damage_rate(double distance, const std::string& material) const;
    double atomic_displacement_rate(double distance, const std::string& material) const;
    
    // Advanced dynamics
    double orbital_perturbation_acceleration(double distance, double mass_object) const;
    std::array<double, 3> poynting_robertson_drag(double distance, double particle_radius, 
                                                 double particle_density) const;
    
    // Environmental hazard assessment
    double total_radiation_dose_rate(double distance, double shielding_thickness = 0.0) const;
    double micrometeorite_impact_rate(double distance, double surface_area) const;
    
    // Utility functions
    void set_random_seed(unsigned int seed);
    double interpolate_solar_activity(double time_years) const;
    std::vector<double> time_series_prediction(double start_time, double end_time, 
                                             double dt, const std::string& parameter) const;
};

#endif // SOLAR_CALCULATIONS_H