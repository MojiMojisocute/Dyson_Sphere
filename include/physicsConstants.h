#ifndef PHYSICS_CONSTANTS_H
#define PHYSICS_CONSTANTS_H

// Standard library includes
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <map>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <string>

// Physical Constants (CODATA 2018 values)
namespace PhysicsConstants {
    // Fundamental constants
    constexpr double G = 6.67430e-11;        // Gravitational constant (m³/kg·s²)
    constexpr double c = 299792458.0;        // Speed of light in vacuum (m/s)
    constexpr double h = 6.62607015e-34;     // Planck constant (J·s)
    constexpr double h_bar = h / (2.0 * M_PI); // Reduced Planck constant (J·s)
    constexpr double k_B = 1.380649e-23;     // Boltzmann constant (J/K)
    constexpr double N_A = 6.02214076e23;    // Avogadro constant (1/mol)
    constexpr double R = 8.314462618;        // Gas constant (J/mol·K)
    
    // Electromagnetic constants
    constexpr double e = 1.602176634e-19;    // Elementary charge (C)
    constexpr double epsilon_0 = 8.8541878128e-12; // Vacuum permittivity (F/m)
    constexpr double mu_0 = 1.25663706212e-6;      // Vacuum permeability (H/m)
    constexpr double alpha = 7.2573532e-3;   // Fine structure constant (dimensionless)
    
    // Particle masses
    constexpr double m_e = 9.1093837015e-31;  // Electron mass (kg)
    constexpr double m_p = 1.672621898e-27;   // Proton mass (kg)
    constexpr double m_n = 1.674927498e-27;   // Neutron mass (kg)
    constexpr double m_u = 1.66053906660e-27; // Atomic mass unit (kg)
    constexpr double m_alpha = 6.6446573357e-27; // Alpha particle mass (kg)
    
    // Radiation constants
    constexpr double sigma = 5.670374419e-8;  // Stefan-Boltzmann constant (W/m²·K⁴)
    constexpr double a = 7.5657e-16;          // Radiation density constant (J/m³·K⁴)
    constexpr double Wien_b = 2.8977729e-3;   // Wien's displacement constant (m·K)
    
    // Astronomical constants
    constexpr double AU = 1.495978707e11;     // Astronomical Unit (m)
    constexpr double pc = 3.0857e16;          // Parsec (m)
    constexpr double ly = 9.4607e15;          // Light year (m)
    constexpr double year = 3.15576e7;        // Year in seconds (s)
    constexpr double day = 86400.0;           // Day in seconds (s)
    
    // Mathematical constants (for convenience)
    constexpr double pi = M_PI;
    constexpr double e_euler = M_E;
    constexpr double sqrt_2 = M_SQRT2;
    constexpr double sqrt_pi = 1.7724538509;
    
    // Unit conversions
    constexpr double eV_to_J = 1.602176634e-19;  // Electron volt to joule
    constexpr double amu_to_kg = 1.66053906660e-27; // Atomic mass unit to kg
    constexpr double barn = 1e-28;               // Barn (nuclear cross-section unit) in m²
    constexpr double gauss_to_tesla = 1e-4;      // Gauss to Tesla conversion
    constexpr double deg_to_rad = M_PI / 180.0;  // Degrees to radians
    
    // Derived constants for solar physics
    constexpr double classical_electron_radius = e * e / (4.0 * pi * epsilon_0 * m_e * c * c); // m
    constexpr double thomson_cross_section = (8.0 * pi / 3.0) * classical_electron_radius * classical_electron_radius; // m²
    constexpr double bohr_radius = 4.0 * pi * epsilon_0 * h_bar * h_bar / (m_e * e * e); // m
    constexpr double rydberg_energy = m_e * e * e * e * e / (8.0 * epsilon_0 * epsilon_0 * h * h); // J
}

// Material property constants
namespace MaterialConstants {
    // Thermal expansion coefficients (1/K)
    constexpr double alpha_aluminum = 23.1e-6;
    constexpr double alpha_silicon = 2.6e-6;
    constexpr double alpha_iron = 11.8e-6;
    constexpr double alpha_carbon = 7.1e-6;
    constexpr double alpha_tungsten = 4.5e-6;
    
    // Young's moduli (Pa)
    constexpr double E_aluminum = 70e9;
    constexpr double E_silicon = 130e9;
    constexpr double E_iron = 200e9;
    constexpr double E_carbon = 500e9;
    constexpr double E_tungsten = 411e9;
    
    // Densities (kg/m³)
    constexpr double rho_aluminum = 2700.0;
    constexpr double rho_silicon = 2330.0;
    constexpr double rho_iron = 7874.0;
    constexpr double rho_carbon = 2267.0;  // Graphite
    constexpr double rho_tungsten = 19250.0;
    
    // Melting points (K)
    constexpr double T_melt_aluminum = 933.47;
    constexpr double T_melt_silicon = 1687.0;
    constexpr double T_melt_iron = 1811.0;
    constexpr double T_melt_carbon = 4098.0;  // Graphite sublimation
    constexpr double T_melt_tungsten = 3695.0;
}

// Atomic data constants
namespace AtomicConstants {
    // Atomic numbers
    constexpr int Z_hydrogen = 1;
    constexpr int Z_helium = 2;
    constexpr int Z_carbon = 6;
    constexpr int Z_oxygen = 8;
    constexpr int Z_silicon = 14;
    constexpr int Z_iron = 26;
    constexpr int Z_tungsten = 74;
    
    // Atomic masses (u)
    constexpr double A_hydrogen = 1.008;
    constexpr double A_helium = 4.003;
    constexpr double A_carbon = 12.011;
    constexpr double A_oxygen = 15.999;
    constexpr double A_silicon = 28.085;
    constexpr double A_iron = 55.845;
    constexpr double A_tungsten = 183.84;
    
    // Ionization energies (eV)
    constexpr double I_hydrogen = 13.6;
    constexpr double I_helium = 24.6;
    constexpr double I_carbon = 11.3;
    constexpr double I_oxygen = 13.6;
    constexpr double I_silicon = 8.2;
    constexpr double I_iron = 7.9;
    constexpr double I_tungsten = 7.98;
}

#endif // PHYSICS_CONSTANTS_H