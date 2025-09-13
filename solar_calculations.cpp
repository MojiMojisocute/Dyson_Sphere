#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <thread>
#include <signal.h>

namespace Constants {
    constexpr double G = 6.67430e-11;           // m³/kg⋅s²
    constexpr double c = 2.99792458e8;          // m/s
    constexpr double h = 6.62607015e-34;       // J⋅s
    constexpr double k_B = 1.380649e-23;       // J/K
    constexpr double sigma = 5.670374419e-8;   // W/m²⋅K⁴ Stefan-Boltzmann
    constexpr double m_p = 1.67262192369e-27;  // kg proton mass
    constexpr double m_e = 9.1093837015e-31;   // kg electron mass
    constexpr double e = 1.602176634e-19;      // C elementary charge
    constexpr double AU = 1.495978707e11;      // m
}

namespace Solar {
    constexpr double mass = 1.989e30;          // kg
    constexpr double radius = 6.96e8;          // m
    constexpr double T_surface = 5778.0;       // K
    constexpr double T_core = 1.57e7;          // K
    constexpr double luminosity = 3.828e26;    // W
    constexpr double wind_speed = 4e5;         // m/s
    constexpr double wind_density = 5e6;       // particles/m³ at 1 AU
    constexpr double cycle_period = 11.0;      // years
    constexpr double surface_B = 1e-4;         // T (Tesla)
}

struct SpectralBand {
    double wavelength_min, wavelength_max;     // m
    double flux_fraction;                      // fraction of total flux
    double opacity;                            // absorption coefficient
};

struct RadiationPressure {
    double total;                              // Pa
    double photon_momentum;                    // kg⋅m/s per m²⋅s
    double absorption_factor;                  // 0-1
    double reflection_factor;                  // 0-1
};

struct ThermalProperties {
    double heat_capacity;                      // J/kg⋅K
    double thermal_conductivity;               // W/m⋅K
    double emissivity;                         // 0-1
    double max_operating_temp;                 // K
};

struct PlasmaInteraction {
    double charge_deposition_rate;             // C/m²⋅s
    double sputtering_yield;                   // atoms/ion
    double erosion_rate;                       // m/s
    double degradation_factor;                 // efficiency loss per year
};

volatile sig_atomic_t keep_running = 1;

void signal_handler(int signal) {
    keep_running = 0;
    std::cout << "\nShutdown requested. Stopping gracefully...\n";
}

class SolarCalculator {
private:
    std::mt19937 rng;
    std::vector<SpectralBand> spectral_bands;
    double simulation_time;                   
    double timestep;                          
    
    double flare_probability_per_hour;
    double cme_probability_per_day;
    
public:
    SolarCalculator() : rng(std::chrono::steady_clock::now().time_since_epoch().count()),
                       simulation_time(0.0), timestep(60.0), 
                       flare_probability_per_hour(0.1),
                       cme_probability_per_day(0.05) {
        initialize_spectral_bands();
    }
    
    void initialize_spectral_bands() {
        spectral_bands = {
            {1e-8, 4e-7, 0.07, 0.8},    // UV: 10nm-400nm, 7% flux
            {4e-7, 7e-7, 0.43, 0.1},    // Visible: 400-700nm, 43% flux
            {7e-7, 1e-4, 0.50, 0.05}    // IR: 700nm-100μm, 50% flux
        };
    }

    double calculate_surface_gravity() const {
        return Constants::G * Solar::mass / (Solar::radius * Solar::radius);
    }
    
    double calculate_escape_velocity() const {
        return sqrt(2.0 * Constants::G * Solar::mass / Solar::radius);
    }
    
    double calculate_schwarzschild_radius() const {
        return 2.0 * Constants::G * Solar::mass / (Constants::c * Constants::c);
    }

    double calculate_fusion_mass_loss_rate() const {
        // From E=mc², luminosity = mass_loss * c²
        return Solar::luminosity / (Constants::c * Constants::c); // kg/s
    }
    
    double calculate_solar_wind_mass_loss() const {
        double surface_area = 4.0 * M_PI * Solar::radius * Solar::radius;
        return Solar::wind_density * Constants::m_p * Solar::wind_speed * surface_area; // kg/s
    }
    
    // Radiation calculations
    double calculate_flux_at_distance(double distance) const {
        return Solar::luminosity / (4.0 * M_PI * distance * distance); // W/m²
    }
    
    std::vector<double> calculate_spectral_flux(double distance) const {
        std::vector<double> fluxes;
        double total_flux = calculate_flux_at_distance(distance);
        
        for (const auto& band : spectral_bands) {
            fluxes.push_back(total_flux * band.flux_fraction);
        }
        return fluxes;
    }

    RadiationPressure calculate_radiation_pressure(double distance, 
                                                  double absorption = 0.9,
                                                  double reflection = 0.1) const {
        RadiationPressure rp;
        double flux = calculate_flux_at_distance(distance);
        
        rp.absorption_factor = absorption;
        rp.reflection_factor = reflection;
        rp.photon_momentum = flux / Constants::c;
        rp.total = (absorption + reflection) * flux / Constants::c; // Pa
        
        return rp;
    }

    double calculate_magnetic_field_strength(double solar_cycle_phase) const {
        return Solar::surface_B * (1.0 + 0.5 * sin(2.0 * M_PI * solar_cycle_phase));
    }

    bool check_solar_flare(double dt_hours) const {
        double probability = 1.0 - exp(-flare_probability_per_hour * dt_hours);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(const_cast<std::mt19937&>(rng)) < probability;
    }
    
    bool check_cme_event(double dt_days) const {
        double probability = 1.0 - exp(-cme_probability_per_day * dt_days);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(const_cast<std::mt19937&>(rng)) < probability;
    }
    
    double calculate_flare_energy_enhancement() const {
        std::lognormal_distribution<double> dist(0.0, 1.0);
        return 1.0 + dist(const_cast<std::mt19937&>(rng)); // 1x to ~10x enhancement
    }
    
    PlasmaInteraction calculate_plasma_interaction(double distance) const {
        PlasmaInteraction pi;
        
        double wind_flux = Solar::wind_density * Solar::wind_speed / (distance/Solar::radius);
        
        pi.charge_deposition_rate = wind_flux * Constants::e * 0.1; // C/m²⋅s
        pi.sputtering_yield = std::min(0.1, wind_flux * 1e-15); // atoms/ion
        pi.erosion_rate = pi.sputtering_yield * wind_flux * 1e-20; // m/s
        pi.degradation_factor = pi.erosion_rate * 365.25 * 24 * 3600 / 1e-6; // per year
        
        return pi;
    }

    ThermalProperties calculate_thermal_properties(double distance, 
                                                  double absorption_efficiency = 0.9) const {
        ThermalProperties tp;
        double flux = calculate_flux_at_distance(distance);
        
        tp.heat_capacity = 900.0; // J/kg⋅K (typical for solar panel materials)
        tp.thermal_conductivity = 150.0; // W/m⋅K
        tp.emissivity = 0.85;
        
        double absorbed_power = flux * absorption_efficiency;
        double equilibrium_temp = pow(absorbed_power / (Constants::sigma * tp.emissivity), 0.25);
        
        tp.max_operating_temp = equilibrium_temp * 1.2; // 20% safety margin
        
        return tp;
    }

    std::map<std::string, double> calculate_comprehensive_state(double distance = Constants::AU,
                                                               double current_time = 0.0) const {
        std::map<std::string, double> state;

        state["mass_kg"] = Solar::mass;
        state["radius_m"] = Solar::radius;
        state["surface_temperature_K"] = Solar::T_surface;
        state["core_temperature_K"] = Solar::T_core;
        state["luminosity_W"] = Solar::luminosity;

        state["surface_gravity_ms2"] = calculate_surface_gravity();
        state["escape_velocity_ms"] = calculate_escape_velocity();
        state["schwarzschild_radius_m"] = calculate_schwarzschild_radius();

        state["fusion_mass_loss_rate_kgs"] = calculate_fusion_mass_loss_rate();
        state["wind_mass_loss_rate_kgs"] = calculate_solar_wind_mass_loss();
        state["total_mass_loss_rate_kgs"] = state["fusion_mass_loss_rate_kgs"] + 
                                           state["wind_mass_loss_rate_kgs"];

        state["distance_m"] = distance;
        state["flux_Wm2"] = calculate_flux_at_distance(distance);

        double cycle_phase = fmod(current_time / (365.25 * 24 * 3600), Solar::cycle_period) / Solar::cycle_period;
        state["solar_cycle_phase"] = cycle_phase;
        state["magnetic_field_T"] = calculate_magnetic_field_strength(cycle_phase);

        auto rp = calculate_radiation_pressure(distance);
        state["radiation_pressure_Pa"] = rp.total;
        state["photon_momentum_kgms_per_m2s"] = rp.photon_momentum;

        auto pi = calculate_plasma_interaction(distance);
        state["charge_deposition_Cm2s"] = pi.charge_deposition_rate;
        state["sputtering_yield"] = pi.sputtering_yield;
        state["erosion_rate_ms"] = pi.erosion_rate;
        state["degradation_factor_per_year"] = pi.degradation_factor;

        auto tp = calculate_thermal_properties(distance);
        state["equilibrium_temperature_K"] = tp.max_operating_temp / 1.2;
        state["max_operating_temp_K"] = tp.max_operating_temp;

        state["wind_speed_ms"] = Solar::wind_speed;
        state["wind_density_particles_m3"] = Solar::wind_density;
        state["wind_pressure_Pa"] = Solar::wind_density * Constants::m_p * 
                                   Solar::wind_speed * Solar::wind_speed;
        
        return state;
    }

    void export_realtime_data(const std::map<std::string, double>& state) const {
        std::ofstream json_file("solar_realtime.json");
        json_file << std::setprecision(10) << std::scientific;
        json_file << "{\n";
        
        bool first = true;
        for (const auto& [key, value] : state) {
            if (!first) json_file << ",\n";
            json_file << "  \"" << key << "\": " << value;
            first = false;
        }
        json_file << "\n}\n";
        json_file.close();

        std::ofstream compact_file("solar_stream.txt");
        compact_file << std::fixed << std::setprecision(6)
                     << state.at("timestamp") << " "
                     << state.at("flux_Wm2") << " "
                     << state.at("radiation_pressure_Pa") << " "
                     << state.at("degradation_factor_per_year") << " "
                     << state.at("equilibrium_temperature_K") << " "
                     << state.at("charge_deposition_Cm2s") << " "
                     << state.at("erosion_rate_ms") << " "
                     << state.at("solar_flare_active") << " "
                     << state.at("cme_active") << " "
                     << state.at("flux_multiplier") << "\n";
        compact_file.close();
    }
};

int main(int argc, char* argv[]) {
    signal(SIGINT, signal_handler);
    
    SolarCalculator calc;
    
    double distance = Constants::AU; // Default 1 AU
    double update_interval = 1.0;    // Default 1 second
    
    // Parse command line arguments
    if (argc > 1) distance = std::stod(argv[1]) * Constants::AU;
    if (argc > 2) update_interval = std::stod(argv[2]);
    
    std::cout << "SOLAR PHYSICS CALCULATOR - STREAMING MODE\n";
    std::cout << "==========================================\n";
    std::cout << "Distance: " << distance/Constants::AU << " AU\n";
    std::cout << "Update interval: " << update_interval << " seconds\n";
    std::cout << "Output files: solar_realtime.json, solar_stream.txt\n";
    std::cout << "Press Ctrl+C to stop...\n\n";
    
    int iteration = 0;
    auto start_program = std::chrono::high_resolution_clock::now();
    
    while (keep_running) {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        double current_time = iteration * update_interval;
        auto state = calc.calculate_comprehensive_state(distance, current_time);
        
        bool flare = calc.check_solar_flare(update_interval / 3600.0);
        bool cme = calc.check_cme_event(update_interval / 86400.0);
        double flux_multiplier = 1.0;
        
        if (flare) {
            flux_multiplier = calc.calculate_flare_energy_enhancement();
            state["flux_Wm2"] *= flux_multiplier;
            state["radiation_pressure_Pa"] *= flux_multiplier;
        }

        state["timestamp"] = current_time;
        state["iteration"] = iteration;
        state["flux_multiplier"] = flux_multiplier;
        state["solar_flare_active"] = flare ? 1.0 : 0.0;
        state["cme_active"] = cme ? 1.0 : 0.0;

        calc.export_realtime_data(state);

        std::cout << "T+" << std::setw(8) << std::fixed << std::setprecision(1) << current_time 
                  << "s | Flux: " << std::scientific << std::setprecision(3) << state["flux_Wm2"] 
                  << " W/m² | Pressure: " << state["radiation_pressure_Pa"] << " Pa";
        
        if (flare) std::cout << " | FLARE x" << std::fixed << std::setprecision(1) << flux_multiplier;
        if (cme) std::cout << " | CME";
        
        std::cout << "                    \r" << std::flush;
        
        iteration++;

        auto end_time = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        auto sleep_time = std::chrono::milliseconds(static_cast<int>(update_interval * 1000)) - elapsed;
        
        if (sleep_time.count() > 0) {
            std::this_thread::sleep_for(sleep_time);
        }
    }
    
    auto end_program = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::seconds>(end_program - start_program);
    
    std::cout << "\n\nFinal Statistics:\n";
    std::cout << "- Total iterations: " << iteration << "\n";
    std::cout << "- Runtime: " << total_time.count() << " seconds\n";
    std::cout << "- Average FPS: " << iteration / double(total_time.count()) << "\n";
    std::cout << "- Data files updated " << iteration << " times\n";
    std::cout << "Ready for fuzzy_controller.py and dyson_model.py integration!\n";
    
    return 0;
}