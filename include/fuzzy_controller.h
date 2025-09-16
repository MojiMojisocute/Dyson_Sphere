#ifndef FUZZY_CONTROLLER_H
#define FUZZY_CONTROLLER_H

#include "shared_memory_manager.h"
#include "physicsConstants.h"
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <functional>

// Fuzzy set definition
struct FuzzySet {
    std::string name;
    std::vector<double> parameters;  // For triangular/trapezoidal: [a, b, c, d]
    enum Type { TRIANGULAR, TRAPEZOIDAL, GAUSSIAN } type;
    
    // Calculate membership value
    double membership(double x) const;
};

// Linguistic variable definition
struct LinguisticVariable {
    std::string name;
    double min_range;
    double max_range;
    std::vector<FuzzySet> sets;
    
    // Get membership values for all sets
    std::vector<double> get_memberships(double value) const;
};

// Fuzzy rule definition
struct FuzzyRule {
    uint32_t id;
    std::vector<std::pair<std::string, std::string>> antecedents; // variable_name, set_name
    std::vector<std::pair<std::string, std::string>> consequents; // variable_name, set_name
    double weight;
    
    // Calculate rule activation strength
    double calculate_activation(const std::map<std::string, double>& inputs,
                               const std::map<std::string, LinguisticVariable>& variables) const;
};

// Input/Output structures
struct FuzzyControlInput {
    double energy_level_normalized;      // 0-1
    double temperature_normalized;       // 0-1 
    double structural_integrity;         // 0-1
    double orbital_stability_index;      // 0-1
    double solar_activity_level;         // 0-1
    double thermal_gradient_max;         // 0-1 (normalized)
    double radiation_dose_rate;          // 0-1 (normalized)
    double degradation_rate;             // 0-1 (normalized)
};

struct FuzzyControlOutput {
    double energy_boost_factor;          // 0.5-2.0
    double cooling_power_factor;         // 0.5-2.0
    double orbit_adjustment_delta;       // -0.1 to +0.1
    double degradation_mitigation;       // 0-1
    double emergency_shutdown;           // 0-1
    double maintenance_priority;         // 0-1
};

class FuzzyController {
private:
    // Linguistic variables for inputs
    std::map<std::string, LinguisticVariable> input_variables_;
    std::map<std::string, LinguisticVariable> output_variables_;
    
    // Rule base
    std::vector<FuzzyRule> rules_;
    
    // Shared memory interface
    std::unique_ptr<SharedMemoryManager> shm_manager_;
    
    // Internal state
    bool initialized_;
    uint32_t rule_counter_;
    
    // Private methods
    void initialize_input_variables();
    void initialize_output_variables();
    void initialize_rule_base();
    
    // Helper methods for rule creation
    void add_energy_rules();
    void add_thermal_rules();
    void add_structural_rules();
    void add_orbital_rules();
    void add_emergency_rules();
    void add_maintenance_rules();
    
    // Inference engine methods
    std::map<std::string, std::vector<double>> fuzzification(const FuzzyControlInput& input);
    std::map<std::string, std::vector<double>> inference(
        const std::map<std::string, std::vector<double>>& fuzzy_inputs);
    FuzzyControlOutput defuzzification(
        const std::map<std::string, std::vector<double>>& fuzzy_outputs);
    
    // Utility methods
    double centroid_defuzzify(const std::vector<double>& membership_values,
                             const LinguisticVariable& variable) const;
    double normalize_input(double value, double min_val, double max_val) const;
    double denormalize_output(double value, double min_val, double max_val) const;
    
    // Rule addition helper
    void add_rule(const std::vector<std::pair<std::string, std::string>>& antecedents,
                  const std::vector<std::pair<std::string, std::string>>& consequents,
                  double weight = 1.0);

public:
    FuzzyController();
    ~FuzzyController();
    
    // Initialization
    bool initialize();
    bool is_initialized() const { return initialized_; }
    
    // Main control method
    FuzzyControlOutput compute_control_output(const FuzzyControlInput& input);
    
    // Rule management
    void add_custom_rule(const FuzzyRule& rule);
    void clear_rules();
    size_t get_rule_count() const { return rules_.size(); }
    
    // Variable access
    const LinguisticVariable* get_input_variable(const std::string& name) const;
    const LinguisticVariable* get_output_variable(const std::string& name) const;
    
    // Diagnostic methods
    void print_fuzzy_sets() const;
    void print_rules() const;
    void print_inference_trace(const FuzzyControlInput& input) const;
    
    // Tuning methods
    void set_rule_weight(uint32_t rule_id, double weight);
    void adjust_membership_function(const std::string& variable_name,
                                   const std::string& set_name,
                                   const std::vector<double>& new_parameters);
    
    // Performance monitoring
    struct PerformanceMetrics {
        double computation_time_ms;
        uint32_t rules_fired;
        double max_rule_activation;
        double min_rule_activation;
        std::map<std::string, double> output_confidence;
    };
    
    PerformanceMetrics get_performance_metrics() const { return last_metrics_; }
    
private:
    mutable PerformanceMetrics last_metrics_;
};

// Utility functions for membership function creation
namespace FuzzyUtils {
    FuzzySet create_triangular_set(const std::string& name, double a, double b, double c);
    FuzzySet create_trapezoidal_set(const std::string& name, double a, double b, double c, double d);
    FuzzySet create_gaussian_set(const std::string& name, double center, double sigma);
    
    // Common membership function shapes
    LinguisticVariable create_low_medium_high_variable(const std::string& name, 
                                                      double min_range, double max_range);
    LinguisticVariable create_five_level_variable(const std::string& name,
                                                 double min_range, double max_range,
                                                 const std::vector<std::string>& level_names);
    
    namespace ThaiSets {
        extern const std::vector<std::string> ENERGY_LEVELS;      // ["Very Low", "Low", "Medium", "High", "Very High"]
        extern const std::vector<std::string> TEMPERATURE_LEVELS; // ["Cold", "Comfortable", "Hot", "Critical"]
        extern const std::vector<std::string> STRUCTURAL_LEVELS;  // ["Severely Damaged", "Degraded", "Normal", "Strong"]
        extern const std::vector<std::string> ORBITAL_LEVELS;     // ["Unstable", "Risky", "Stable", "Excellent"]
        extern const std::vector<std::string> SOLAR_LEVELS;       // ["Calm", "Moderate", "Turbulent", "Severe"]

        extern const std::vector<std::string> ENERGY_BOOST;       // ["Reduced Output", "Normal", "Slight Increase", "Significant Increase", "Maximum Power"]
        extern const std::vector<std::string> COOLING_POWER;      // ["Off", "Low", "Medium", "High", "Maximum"]
        extern const std::vector<std::string> ORBIT_ADJUST;       // ["Move Far Out", "Move Out", "Hold Position", "Move In", "Move Much Closer"]
        extern const std::vector<std::string> MITIGATION_LEVELS;  // ["None", "Minor", "Moderate", "Intense", "Maximum"]
    }
}

#endif // FUZZY_CONTROLLER_H