#include "fuzzy_controller.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <cassert>

using namespace PhysicsConstants;

// English language fuzzy set names
namespace FuzzyUtils::EnglishSets {
    const std::vector<std::string> ENERGY_LEVELS = {"Very Low", "Low", "Medium", "High", "Very High"};
    const std::vector<std::string> TEMPERATURE_LEVELS = {"Cold", "Comfortable", "Hot", "Critical"};
    const std::vector<std::string> STRUCTURAL_LEVELS = {"Severely Damaged", "Degraded", "Normal", "Strong"};
    const std::vector<std::string> ORBITAL_LEVELS = {"Unstable", "Risky", "Stable", "Excellent"};
    const std::vector<std::string> SOLAR_LEVELS = {"Calm", "Moderate", "Turbulent", "Severe"};
    
    const std::vector<std::string> ENERGY_BOOST = {"Reduced Output", "Normal", "Slight Increase", "Significant Increase", "Maximum Power"};
    const std::vector<std::string> COOLING_POWER = {"Off", "Low", "Medium", "High", "Maximum"};
    const std::vector<std::string> ORBIT_ADJUST = {"Move Far Out", "Move Out", "Hold Position", "Move In", "Move Much Closer"};
    const std::vector<std::string> MITIGATION_LEVELS = {"None", "Minor", "Moderate", "Intense", "Maximum"};
}

// FuzzySet implementation
double FuzzySet::membership(double x) const {
    switch (type) {
        case TRIANGULAR: {
            if (parameters.size() != 3) return 0.0;
            double a = parameters[0], b = parameters[1], c = parameters[2];
            
            if (x <= a || x >= c) return 0.0;
            if (x == b) return 1.0;
            if (x < b) return (x - a) / (b - a);
            else return (c - x) / (c - b);
        }
        
        case TRAPEZOIDAL: {
            if (parameters.size() != 4) return 0.0;
            double a = parameters[0], b = parameters[1], c = parameters[2], d = parameters[3];
            
            if (x <= a || x >= d) return 0.0;
            if (x >= b && x <= c) return 1.0;
            if (x < b) return (x - a) / (b - a);
            else return (d - x) / (d - c);
        }
        
        case GAUSSIAN: {
            if (parameters.size() != 2) return 0.0;
            double center = parameters[0], sigma = parameters[1];
            double diff = x - center;
            return std::exp(-(diff * diff) / (2.0 * sigma * sigma));
        }
    }
    return 0.0;
}

// LinguisticVariable implementation
std::vector<double> LinguisticVariable::get_memberships(double value) const {
    std::vector<double> memberships;
    memberships.reserve(sets.size());
    
    for (const auto& set : sets) {
        memberships.push_back(set.membership(value));
    }
    
    return memberships;
}

// FuzzyRule implementation
double FuzzyRule::calculate_activation(const std::map<std::string, double>& inputs,
                                     const std::map<std::string, LinguisticVariable>& variables) const {
    double activation = 1.0;
    
    for (const auto& antecedent : antecedents) {
        const std::string& var_name = antecedent.first;
        const std::string& set_name = antecedent.second;
        
        auto var_it = variables.find(var_name);
        if (var_it == variables.end()) return 0.0;
        
        auto input_it = inputs.find(var_name);
        if (input_it == inputs.end()) return 0.0;
        
        // Find the fuzzy set
        const auto& sets = var_it->second.sets;
        auto set_it = std::find_if(sets.begin(), sets.end(),
                                  [&set_name](const FuzzySet& s) { return s.name == set_name; });
        
        if (set_it == sets.end()) return 0.0;
        
        double membership = set_it->membership(input_it->second);
        activation = std::min(activation, membership); // AND operation (minimum)
    }
    
    return activation * weight;
}

// FuzzyController implementation
FuzzyController::FuzzyController() 
    : initialized_(false), rule_counter_(0) {
    try {
        shm_manager_ = std::make_unique<SharedMemoryManager>(SharedMemoryManager::AccessMode::CONSUMER);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Could not initialize shared memory: " << e.what() << std::endl;
        shm_manager_ = nullptr;
    }
}

FuzzyController::~FuzzyController() = default;

bool FuzzyController::initialize() {
    if (initialized_) return true;
    
    try {
        initialize_input_variables();
        initialize_output_variables();
        initialize_rule_base();
        
        initialized_ = true;
        
        std::cout << "Fuzzy Controller initialized successfully!" << std::endl;
        std::cout << "Input variables: " << input_variables_.size() << std::endl;
        std::cout << "Output variables: " << output_variables_.size() << std::endl;
        std::cout << "Rules: " << rules_.size() << std::endl;
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Failed to initialize Fuzzy Controller: " << e.what() << std::endl;
        return false;
    }
}

void FuzzyController::initialize_input_variables() {
    // 1. Energy Level
    input_variables_["energy"] = FuzzyUtils::create_five_level_variable(
        "energy", 0.0, 1.0, FuzzyUtils::EnglishSets::ENERGY_LEVELS);
    
    // 2. Temperature
    LinguisticVariable temp_var;
    temp_var.name = "temperature";
    temp_var.min_range = 0.0;
    temp_var.max_range = 1.0;
    temp_var.sets = {
        FuzzyUtils::create_trapezoidal_set("Cold", 0.0, 0.0, 0.3, 0.5),
        FuzzyUtils::create_triangular_set("Comfortable", 0.3, 0.5, 0.7),
        FuzzyUtils::create_triangular_set("Hot", 0.5, 0.7, 0.9),
        FuzzyUtils::create_trapezoidal_set("Critical", 0.7, 0.9, 1.0, 1.0)
    };
    input_variables_["temperature"] = temp_var;
    
    // 3. Structural Integrity
    LinguisticVariable struct_var;
    struct_var.name = "structural";
    struct_var.min_range = 0.0;
    struct_var.max_range = 1.0;
    struct_var.sets = {
        FuzzyUtils::create_trapezoidal_set("Severely Damaged", 0.0, 0.0, 0.2, 0.4),
        FuzzyUtils::create_triangular_set("Degraded", 0.2, 0.4, 0.6),
        FuzzyUtils::create_triangular_set("Normal", 0.4, 0.6, 0.8),
        FuzzyUtils::create_trapezoidal_set("Strong", 0.6, 0.8, 1.0, 1.0)
    };
    input_variables_["structural"] = struct_var;
    
    // 4. Orbital Stability
    LinguisticVariable orbital_var;
    orbital_var.name = "orbital";
    orbital_var.min_range = 0.0;
    orbital_var.max_range = 1.0;
    orbital_var.sets = {
        FuzzyUtils::create_trapezoidal_set("Unstable", 0.0, 0.0, 0.25, 0.45),
        FuzzyUtils::create_triangular_set("Risky", 0.25, 0.45, 0.65),
        FuzzyUtils::create_triangular_set("Stable", 0.45, 0.65, 0.85),
        FuzzyUtils::create_trapezoidal_set("Excellent", 0.65, 0.85, 1.0, 1.0)
    };
    input_variables_["orbital"] = orbital_var;
    
    // 5. Solar Activity
    LinguisticVariable solar_var;
    solar_var.name = "solar_activity";
    solar_var.min_range = 0.0;
    solar_var.max_range = 1.0;
    solar_var.sets = {
        FuzzyUtils::create_trapezoidal_set("Calm", 0.0, 0.0, 0.2, 0.4),
        FuzzyUtils::create_triangular_set("Moderate", 0.2, 0.4, 0.6),
        FuzzyUtils::create_triangular_set("Turbulent", 0.4, 0.6, 0.8),
        FuzzyUtils::create_trapezoidal_set("Severe", 0.6, 0.8, 1.0, 1.0)
    };
    input_variables_["solar_activity"] = solar_var;
    
    // Additional inputs
    input_variables_["thermal_gradient"] = FuzzyUtils::create_low_medium_high_variable(
        "thermal_gradient", 0.0, 1.0);
    input_variables_["radiation"] = FuzzyUtils::create_low_medium_high_variable(
        "radiation", 0.0, 1.0);
    input_variables_["degradation"] = FuzzyUtils::create_low_medium_high_variable(
        "degradation", 0.0, 1.0);
}

void FuzzyController::initialize_output_variables() {
    // 1. Energy Boost (mapped to 0.5-2.0 range)
    output_variables_["energy_boost"] = FuzzyUtils::create_five_level_variable(
        "energy_boost", 0.0, 1.0, FuzzyUtils::EnglishSets::ENERGY_BOOST);
    
    // 2. Cooling Power (mapped to 0.5-2.0 range)
    output_variables_["cooling_power"] = FuzzyUtils::create_five_level_variable(
        "cooling_power", 0.0, 1.0, FuzzyUtils::EnglishSets::COOLING_POWER);
    
    // 3. Orbit Adjustment (mapped to -0.1 to +0.1 range)
    output_variables_["orbit_adjustment"] = FuzzyUtils::create_five_level_variable(
        "orbit_adjustment", 0.0, 1.0, FuzzyUtils::EnglishSets::ORBIT_ADJUST);
    
    // 4. Degradation Mitigation
    output_variables_["mitigation"] = FuzzyUtils::create_five_level_variable(
        "mitigation", 0.0, 1.0, FuzzyUtils::EnglishSets::MITIGATION_LEVELS);
    
    // 5. Emergency Shutdown
    LinguisticVariable emergency_var;
    emergency_var.name = "emergency";
    emergency_var.min_range = 0.0;
    emergency_var.max_range = 1.0;
    emergency_var.sets = {
        FuzzyUtils::create_trapezoidal_set("Normal", 0.0, 0.0, 0.3, 0.5),
        FuzzyUtils::create_triangular_set("Alert", 0.3, 0.5, 0.7),
        FuzzyUtils::create_trapezoidal_set("Emergency", 0.5, 0.7, 1.0, 1.0)
    };
    output_variables_["emergency"] = emergency_var;
    
    // 6. Maintenance Priority
    output_variables_["maintenance"] = FuzzyUtils::create_low_medium_high_variable(
        "maintenance", 0.0, 1.0);
}

void FuzzyController::initialize_rule_base() {
    rule_counter_ = 0;
    rules_.clear();
    
    add_energy_rules();
    add_thermal_rules();
    add_structural_rules();
    add_orbital_rules();
    add_emergency_rules();
    add_maintenance_rules();
}

void FuzzyController::add_energy_rules() {
    // Energy management rules
    
    // If energy is low and temperature is comfortable -> increase energy slightly
    add_rule({{"energy", "Low"}, {"temperature", "Comfortable"}},
             {{"energy_boost", "Slight Increase"}}, 0.8);
    
    // If energy is very low and temperature is cold -> increase energy significantly
    add_rule({{"energy", "Very Low"}, {"temperature", "Cold"}},
             {{"energy_boost", "Maximum Power"}}, 0.9);
    
    // If energy is high and temperature is hot -> reduce energy output
    add_rule({{"energy", "High"}, {"temperature", "Hot"}},
             {{"energy_boost", "Reduced Output"}}, 0.7);
    
    // If energy is very high and temperature is critical -> reduce energy output immediately
    add_rule({{"energy", "Very High"}, {"temperature", "Critical"}},
             {{"energy_boost", "Reduced Output"}}, 1.0);
    
    // If solar activity is severe and energy is medium -> increase energy boost
    add_rule({{"solar_activity", "Severe"}, {"energy", "Medium"}},
             {{"energy_boost", "Significant Increase"}}, 0.6);
}

void FuzzyController::add_thermal_rules() {
    // Thermal management rules
    
    // If temperature is hot -> activate high cooling
    add_rule({{"temperature", "Hot"}},
             {{"cooling_power", "High"}}, 0.8);
    
    // If temperature is critical -> activate maximum cooling
    add_rule({{"temperature", "Critical"}},
             {{"cooling_power", "Maximum"}}, 1.0);
    
    // If temperature is cold and energy is low -> turn off cooling
    add_rule({{"temperature", "Cold"}, {"energy", "Low"}},
             {{"cooling_power", "Off"}}, 0.7);
    
    // If thermal gradient is high and temperature is comfortable -> activate medium cooling
    add_rule({{"thermal_gradient", "High"}, {"temperature", "Comfortable"}},
             {{"cooling_power", "Medium"}}, 0.6);
    
    // If solar activity is severe -> prepare high cooling
    add_rule({{"solar_activity", "Severe"}},
             {{"cooling_power", "High"}}, 0.7);
}

void FuzzyController::add_structural_rules() {
    // Structural integrity rules
    
    // If structural integrity is severely damaged -> move far out from sun
    add_rule({{"structural", "Severely Damaged"}},
             {{"orbit_adjustment", "Move Far Out"}, {"mitigation", "Maximum"}}, 0.9);
    
    // If structural integrity is degraded and solar activity is calm -> perform maintenance
    add_rule({{"structural", "Degraded"}, {"solar_activity", "Calm"}},
             {{"mitigation", "Intense"}, {"maintenance", "High"}}, 0.8);
    
    // If structural integrity is strong and energy is low -> move closer to sun
    add_rule({{"structural", "Strong"}, {"energy", "Low"}},
             {{"orbit_adjustment", "Move In"}}, 0.6);
    
    // If radiation is high -> activate mitigation systems
    add_rule({{"radiation", "High"}},
             {{"mitigation", "Intense"}}, 0.7);
}

void FuzzyController::add_orbital_rules() {
    // Orbital dynamics rules
    
    // If orbit is unstable and structural integrity is strong -> adjust orbit outward
    add_rule({{"orbital", "Unstable"}, {"structural", "Strong"}},
             {{"orbit_adjustment", "Move Out"}}, 0.8);
    
    // If orbit is risky and temperature is hot -> move away from sun
    add_rule({{"orbital", "Risky"}, {"temperature", "Hot"}},
             {{"orbit_adjustment", "Move Out"}}, 0.7);
    
    // If orbit is stable and energy is low -> move closer to sun
    add_rule({{"orbital", "Stable"}, {"energy", "Low"}},
             {{"orbit_adjustment", "Move In"}}, 0.6);
    
    // If orbit is excellent and everything is normal -> hold position
    add_rule({{"orbital", "Excellent"}, {"energy", "Medium"}},
             {{"orbit_adjustment", "Hold Position"}}, 0.5);
}

void FuzzyController::add_emergency_rules() {
    // Emergency situation rules - highest priority
    
    // If temperature is critical -> emergency shutdown
    add_rule({{"temperature", "Critical"}},
             {{"emergency", "Emergency"}, {"cooling_power", "Maximum"}, {"energy_boost", "Reduced Output"}}, 1.0);
    
    // If structural integrity is severely damaged -> emergency mode
    add_rule({{"structural", "Severely Damaged"}},
             {{"emergency", "Emergency"}, {"mitigation", "Maximum"}}, 1.0);
    
    // If orbit is unstable and structural integrity is degraded -> emergency mode
    add_rule({{"orbital", "Unstable"}, {"structural", "Degraded"}},
             {{"emergency", "Emergency"}, {"orbit_adjustment", "Move Far Out"}}, 0.9);
    
    // If solar activity is severe and structural integrity is degraded -> alert mode
    add_rule({{"solar_activity", "Severe"}, {"structural", "Degraded"}},
             {{"emergency", "Alert"}, {"mitigation", "Intense"}}, 0.8);
    
    // If radiation is high and energy is low -> alert mode
    add_rule({{"radiation", "High"}, {"energy", "Low"}},
             {{"emergency", "Alert"}}, 0.7);
}

void FuzzyController::add_maintenance_rules() {
    // Maintenance priority rules
    
    // If structural integrity is degraded -> high maintenance priority
    add_rule({{"structural", "Degraded"}},
             {{"maintenance", "High"}}, 0.8);
    
    // If solar activity is calm and degradation is high -> schedule maintenance
    add_rule({{"solar_activity", "Calm"}, {"degradation", "High"}},
             {{"maintenance", "High"}}, 0.7);
    
    // If everything is normal -> low maintenance priority
    add_rule({{"structural", "Normal"}, {"energy", "Medium"}, {"temperature", "Comfortable"}},
             {{"maintenance", "Low"}}, 0.5);
    
    // If multiple degradation factors present -> increase maintenance priority
    add_rule({{"degradation", "High"}, {"radiation", "High"}},
             {{"maintenance", "High"}}, 0.8);
}

FuzzyControlOutput FuzzyController::compute_control_output(const FuzzyControlInput& input) {
    if (!initialized_) {
        throw std::runtime_error("Fuzzy Controller not initialized");
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Step 1: Fuzzification
    std::map<std::string, double> input_map = {
        {"energy", input.energy_level_normalized},
        {"temperature", input.temperature_normalized},
        {"structural", input.structural_integrity},
        {"orbital", input.orbital_stability_index},
        {"solar_activity", input.solar_activity_level},
        {"thermal_gradient", normalize_input(input.thermal_gradient_max, 0.0, 100.0)},
        {"radiation", normalize_input(input.radiation_dose_rate, 1e-10, 1e-3)},
        {"degradation", normalize_input(input.degradation_rate, 0.0, 10.0)}
    };
    
    auto fuzzy_inputs = fuzzification(input);
    
    // Step 2: Inference
    auto fuzzy_outputs = inference(fuzzy_inputs);
    
    // Step 3: Defuzzification
    auto output = defuzzification(fuzzy_outputs);
    
    // Record performance metrics
    auto end_time = std::chrono::high_resolution_clock::now();
    last_metrics_.computation_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    return output;
}

std::map<std::string, std::vector<double>> FuzzyController::fuzzification(const FuzzyControlInput& input) {
    std::map<std::string, std::vector<double>> fuzzy_inputs;
    
    std::map<std::string, double> input_values = {
        {"energy", input.energy_level_normalized},
        {"temperature", input.temperature_normalized},
        {"structural", input.structural_integrity},
        {"orbital", input.orbital_stability_index},
        {"solar_activity", input.solar_activity_level},
        {"thermal_gradient", normalize_input(input.thermal_gradient_max, 0.0, 100.0)},
        {"radiation", normalize_input(input.radiation_dose_rate, 1e-10, 1e-3)},
        {"degradation", normalize_input(input.degradation_rate, 0.0, 10.0)}
    };
    
    for (const auto& [var_name, var] : input_variables_) {
        if (input_values.find(var_name) != input_values.end()) {
            fuzzy_inputs[var_name] = var.get_memberships(input_values[var_name]);
        }
    }
    
    return fuzzy_inputs;
}

std::map<std::string, std::vector<double>> FuzzyController::inference(
    const std::map<std::string, std::vector<double>>& fuzzy_inputs) {
    
    std::map<std::string, std::vector<double>> fuzzy_outputs;
    
    // Initialize output membership values
    for (const auto& [var_name, var] : output_variables_) {
        fuzzy_outputs[var_name] = std::vector<double>(var.sets.size(), 0.0);
    }
    
    last_metrics_.rules_fired = 0;
    last_metrics_.max_rule_activation = 0.0;
    last_metrics_.min_rule_activation = 1.0;
    
    // Convert fuzzy_inputs to simple map for rule evaluation
    std::map<std::string, double> crisp_inputs;
    for (const auto& [var_name, memberships] : fuzzy_inputs) {
        // Use center of mass for crisp value
        double sum = 0.0, weighted_sum = 0.0;
        for (size_t i = 0; i < memberships.size(); ++i) {
            weighted_sum += memberships[i] * (double(i) / (memberships.size() - 1));
            sum += memberships[i];
        }
        crisp_inputs[var_name] = (sum > 0) ? weighted_sum / sum : 0.5;
    }
    
    // Process each rule
    for (const auto& rule : rules_) {
        double activation = rule.calculate_activation(crisp_inputs, input_variables_);
        
        if (activation > 0.01) {  // Only process significantly activated rules
            last_metrics_.rules_fired++;
            last_metrics_.max_rule_activation = std::max(last_metrics_.max_rule_activation, activation);
            last_metrics_.min_rule_activation = std::min(last_metrics_.min_rule_activation, activation);
            
            // Apply rule to outputs (Mamdani inference)
            for (const auto& consequent : rule.consequents) {
                const std::string& output_var_name = consequent.first;
                const std::string& output_set_name = consequent.second;
                
                auto output_var_it = output_variables_.find(output_var_name);
                if (output_var_it != output_variables_.end()) {
                    const auto& sets = output_var_it->second.sets;
                    
                    for (size_t i = 0; i < sets.size(); ++i) {
                        if (sets[i].name == output_set_name) {
                            fuzzy_outputs[output_var_name][i] = std::max(
                                fuzzy_outputs[output_var_name][i], activation);
                            break;
                        }
                    }
                }
            }
        }
    }
    
    return fuzzy_outputs;
}

FuzzyControlOutput FuzzyController::defuzzification(
    const std::map<std::string, std::vector<double>>& fuzzy_outputs) {
    
    FuzzyControlOutput output;
    
    // Defuzzify each output variable using centroid method
    double energy_boost = 0.0;
    if (fuzzy_outputs.find("energy_boost") != fuzzy_outputs.end()) {
        energy_boost = centroid_defuzzify(fuzzy_outputs.at("energy_boost"), 
                                        output_variables_.at("energy_boost"));
    }
    
    double cooling_power = 0.0;
    if (fuzzy_outputs.find("cooling_power") != fuzzy_outputs.end()) {
        cooling_power = centroid_defuzzify(fuzzy_outputs.at("cooling_power"), 
                                         output_variables_.at("cooling_power"));
    }
    
    double orbit_adjustment = 0.0;
    if (fuzzy_outputs.find("orbit_adjustment") != fuzzy_outputs.end()) {
        orbit_adjustment = centroid_defuzzify(fuzzy_outputs.at("orbit_adjustment"), 
                                            output_variables_.at("orbit_adjustment"));
    }
    
    double mitigation = 0.0;
    if (fuzzy_outputs.find("mitigation") != fuzzy_outputs.end()) {
        mitigation = centroid_defuzzify(fuzzy_outputs.at("mitigation"), 
                                      output_variables_.at("mitigation"));
    }
    
    double emergency = 0.0;
    if (fuzzy_outputs.find("emergency") != fuzzy_outputs.end()) {
        emergency = centroid_defuzzify(fuzzy_outputs.at("emergency"), 
                                     output_variables_.at("emergency"));
    }
    
    double maintenance = 0.0;
    if (fuzzy_outputs.find("maintenance") != fuzzy_outputs.end()) {
        maintenance = centroid_defuzzify(fuzzy_outputs.at("maintenance"), 
                                       output_variables_.at("maintenance"));
    }
    
    // Convert to output ranges
    output.energy_boost_factor = denormalize_output(energy_boost, 0.5, 2.0);
    output.cooling_power_factor = denormalize_output(cooling_power, 0.5, 2.0);
    output.orbit_adjustment_delta = denormalize_output(orbit_adjustment, -0.1, 0.1);
    output.degradation_mitigation = mitigation;
    output.emergency_shutdown = emergency;
    output.maintenance_priority = maintenance;
    
    // Calculate confidence levels
    last_metrics_.output_confidence["energy_boost"] = 
        *std::max_element(fuzzy_outputs.at("energy_boost").begin(), 
                         fuzzy_outputs.at("energy_boost").end());
    last_metrics_.output_confidence["cooling_power"] = 
        *std::max_element(fuzzy_outputs.at("cooling_power").begin(), 
                         fuzzy_outputs.at("cooling_power").end());
    
    return output;
}

double FuzzyController::centroid_defuzzify(const std::vector<double>& membership_values,
                                         const LinguisticVariable& variable) const {
    double numerator = 0.0;
    double denominator = 0.0;
    
    const int resolution = 1000;  // Resolution for numerical integration
    const double step = (variable.max_range - variable.min_range) / resolution;
    
    for (int i = 0; i < resolution; ++i) {
        double x = variable.min_range + i * step;
        double membership = 0.0;
        
        // Calculate aggregated membership at point x
        for (size_t j = 0; j < variable.sets.size(); ++j) {
            membership = std::max(membership, 
                                std::min(membership_values[j], variable.sets[j].membership(x)));
        }
        
        numerator += x * membership;
        denominator += membership;
    }
    
    return (denominator > 0) ? numerator / denominator : 0.5;
}

double FuzzyController::normalize_input(double value, double min_val, double max_val) const {
    if (max_val <= min_val) return 0.5;
    return std::max(0.0, std::min(1.0, (value - min_val) / (max_val - min_val)));
}

double FuzzyController::denormalize_output(double value, double min_val, double max_val) const {
    return min_val + value * (max_val - min_val);
}

void FuzzyController::add_rule(const std::vector<std::pair<std::string, std::string>>& antecedents,
                              const std::vector<std::pair<std::string, std::string>>& consequents,
                              double weight) {
    FuzzyRule rule;
    rule.id = rule_counter_++;
    rule.antecedents = antecedents;
    rule.consequents = consequents;
    rule.weight = weight;
    
    rules_.push_back(rule);
}

void FuzzyController::print_fuzzy_sets() const {
    std::cout << "\n=== INPUT VARIABLES ===" << std::endl;
    for (const auto& [name, var] : input_variables_) {
        std::cout << name << " (" << var.min_range << "-" << var.max_range << "):" << std::endl;
        for (const auto& set : var.sets) {
            std::cout << "  " << set.name << " ";
            for (double param : set.parameters) {
                std::cout << param << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    std::cout << "\n=== OUTPUT VARIABLES ===" << std::endl;
    for (const auto& [name, var] : output_variables_) {
        std::cout << name << " (" << var.min_range << "-" << var.max_range << "):" << std::endl;
        for (const auto& set : var.sets) {
            std::cout << "  " << set.name << " ";
            for (double param : set.parameters) {
                std::cout << param << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void FuzzyController::print_rules() const {
    std::cout << "\n=== FUZZY RULES (" << rules_.size() << " total) ===" << std::endl;
    
    for (const auto& rule : rules_) {
        std::cout << "Rule " << rule.id << " (weight=" << rule.weight << "): ";
        
        // Print antecedents
        std::cout << "IF ";
        for (size_t i = 0; i < rule.antecedents.size(); ++i) {
            if (i > 0) std::cout << " AND ";
            std::cout << rule.antecedents[i].first << " IS \"" << rule.antecedents[i].second << "\"";
        }
        
        // Print consequents
        std::cout << " THEN ";
        for (size_t i = 0; i < rule.consequents.size(); ++i) {
            if (i > 0) std::cout << " AND ";
            std::cout << rule.consequents[i].first << " IS \"" << rule.consequents[i].second << "\"";
        }
        
        std::cout << std::endl;
    }
}

void FuzzyController::print_inference_trace(const FuzzyControlInput& input) const {
    std::cout << "\n=== INFERENCE TRACE ===" << std::endl;
    
    // Print inputs
    std::cout << "INPUTS:" << std::endl;
    std::cout << "  Energy Level: " << input.energy_level_normalized << std::endl;
    std::cout << "  Temperature: " << input.temperature_normalized << std::endl;
    std::cout << "  Structural: " << input.structural_integrity << std::endl;
    std::cout << "  Orbital: " << input.orbital_stability_index << std::endl;
    std::cout << "  Solar Activity: " << input.solar_activity_level << std::endl;
    
    // Compute and print output
    auto output = const_cast<FuzzyController*>(this)->compute_control_output(input);
    
    std::cout << "\nOUTPUTS:" << std::endl;
    std::cout << "  Energy Boost Factor: " << output.energy_boost_factor << std::endl;
    std::cout << "  Cooling Power Factor: " << output.cooling_power_factor << std::endl;
    std::cout << "  Orbit Adjustment: " << output.orbit_adjustment_delta << std::endl;
    std::cout << "  Mitigation: " << output.degradation_mitigation << std::endl;
    std::cout << "  Emergency: " << output.emergency_shutdown << std::endl;
    std::cout << "  Maintenance: " << output.maintenance_priority << std::endl;
    
    std::cout << "\nPERFORMANCE:" << std::endl;
    std::cout << "  Rules Fired: " << last_metrics_.rules_fired << std::endl;
    std::cout << "  Computation Time: " << last_metrics_.computation_time_ms << " ms" << std::endl;
}

// Utility functions implementation
namespace FuzzyUtils {
    FuzzySet create_triangular_set(const std::string& name, double a, double b, double c) {
        FuzzySet set;
        set.name = name;
        set.type = FuzzySet::TRIANGULAR;
        set.parameters = {a, b, c};
        return set;
    }
    
    FuzzySet create_trapezoidal_set(const std::string& name, double a, double b, double c, double d) {
        FuzzySet set;
        set.name = name;
        set.type = FuzzySet::TRAPEZOIDAL;
        set.parameters = {a, b, c, d};
        return set;
    }
    
    FuzzySet create_gaussian_set(const std::string& name, double center, double sigma) {
        FuzzySet set;
        set.name = name;
        set.type = FuzzySet::GAUSSIAN;
        set.parameters = {center, sigma};
        return set;
    }
    
    LinguisticVariable create_low_medium_high_variable(const std::string& name, 
                                                      double min_range, double max_range) {
        LinguisticVariable var;
        var.name = name;
        var.min_range = min_range;
        var.max_range = max_range;
        
        var.sets = {
            create_trapezoidal_set("Low", 0.0, 0.0, 0.3, 0.5),
            create_triangular_set("Medium", 0.2, 0.5, 0.8),
            create_trapezoidal_set("High", 0.5, 0.7, 1.0, 1.0)
        };
        
        return var;
    }
    
    LinguisticVariable create_five_level_variable(const std::string& name,
                                                 double min_range, double max_range,
                                                 const std::vector<std::string>& level_names) {
        LinguisticVariable var;
        var.name = name;
        var.min_range = min_range;
        var.max_range = max_range;
        
        if (level_names.size() >= 5) {
            var.sets = {
                create_trapezoidal_set(level_names[0], 0.0, 0.0, 0.15, 0.3),   // Very Low/Reduced Output/Off/Move Far Out/None
                create_triangular_set(level_names[1], 0.15, 0.3, 0.45),        // Low/Normal/Low/Move Out/Minor
                create_triangular_set(level_names[2], 0.3, 0.5, 0.7),          // Medium/Slight Increase/Medium/Hold Position/Moderate
                create_triangular_set(level_names[3], 0.55, 0.7, 0.85),        // High/Significant Increase/High/Move In/Intense
                create_trapezoidal_set(level_names[4], 0.7, 0.85, 1.0, 1.0)    // Very High/Maximum Power/Maximum/Move Much Closer/Maximum
            };
        } else {
            // Fallback to 3 levels
            var.sets = {
                create_triangular_set(level_names.size() > 0 ? level_names[0] : "Low", 0.0, 0.0, 0.5),
                create_triangular_set(level_names.size() > 1 ? level_names[1] : "Medium", 0.2, 0.5, 0.8),
                create_triangular_set(level_names.size() > 2 ? level_names[2] : "High", 0.5, 1.0, 1.0)
            };
        }
        
        return var;
    }
}