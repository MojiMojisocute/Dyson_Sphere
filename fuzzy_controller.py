#!/usr/bin/env python3
"""
Improved Fuzzy Logic Controller for Dyson Swarm Management
Fixes: Real control application, clearer visualization, energy accumulation
"""

import numpy as np
import json
import time
import subprocess
import signal
import sys
import os
import threading
from pathlib import Path
from typing import Dict, List, Tuple, Optional

try:
    import skfuzzy as fuzz
    from skfuzzy import control as ctrl
    import matplotlib.pyplot as plt
    FUZZY_AVAILABLE = True
    PLOT_AVAILABLE = True
except ImportError:
    print("Warning: skfuzzy/matplotlib not available. Using rule-based fallback.")
    FUZZY_AVAILABLE = False
    PLOT_AVAILABLE = False

class ImprovedDysonController:
    def __init__(self):
        self.running = True
        self.cpp_processes = []
        self.simulation_data = {}
        self.control_history = []
        
        self.total_energy_collected = 0.0  
        self.energy_collection_rate = 0.0  
        self.peak_power = 0.0
        self.avg_efficiency = 0.0
        self.start_time = time.time()

        self.display_metrics = {
            'power_output_MW': 0.0,
            'total_energy_GWh': 0.0,
            'panels_active': 0,
            'panels_total': 0,
            'avg_temperature_C': 0.0,
            'system_efficiency_percent': 0.0,
            'time_running_hours': 0.0
        }

        self.fuzzy_system = None
        self.fuzzy_simulation = None
        self.last_control_outputs = {
            'energy_boost': 1.0,
            'cooling_power': 1.0,
            'orbit_adjustment': 0.0
        }
        
        if FUZZY_AVAILABLE:
            self.setup_fuzzy_system()
        else:
            print("Using rule-based control fallback")
            
        signal.signal(signal.SIGINT, self.signal_handler)
    
    def signal_handler(self, signum, frame):
        print("\nShutting down Improved Dyson Controller...")
        self.save_final_report()
        self.running = False
        self.cleanup()
        sys.exit(0)
    
    def cleanup(self):
        for proc in self.cpp_processes:
            if proc.poll() is None:
                proc.terminate()
                proc.wait(timeout=5)
    
    def setup_fuzzy_system(self):
        try:
            print("Initializing improved fuzzy control system...")

            panel_efficiency = ctrl.Antecedent(np.arange(0, 101, 1), 'panel_efficiency')  # %
            temperature_C = ctrl.Antecedent(np.arange(-50, 501, 1), 'temperature_C')  # Celsius
            power_demand = ctrl.Antecedent(np.arange(0, 201, 1), 'power_demand')  # % of capacity
            solar_activity = ctrl.Antecedent(np.arange(0, 301, 1), 'solar_activity')  # % normal
            system_health = ctrl.Antecedent(np.arange(0, 101, 1), 'system_health')  # %

            energy_boost = ctrl.Consequent(np.arange(50, 151, 1), 'energy_boost')  # % normal rate
            cooling_power = ctrl.Consequent(np.arange(0, 201, 1), 'cooling_power')  # % max cooling
            orbit_adjustment = ctrl.Consequent(np.arange(-10, 11, 1), 'orbit_adjustment')  # % closer/farther

            panel_efficiency['poor'] = fuzz.trimf(panel_efficiency.universe, [0, 0, 30])
            panel_efficiency['fair'] = fuzz.trimf(panel_efficiency.universe, [20, 50, 80])
            panel_efficiency['excellent'] = fuzz.trimf(panel_efficiency.universe, [70, 100, 100])
            
            temperature_C['cold'] = fuzz.trimf(temperature_C.universe, [-50, -50, 50])
            temperature_C['optimal'] = fuzz.trimf(temperature_C.universe, [25, 75, 125])
            temperature_C['hot'] = fuzz.trimf(temperature_C.universe, [100, 200, 300])
            temperature_C['critical'] = fuzz.trimf(temperature_C.universe, [250, 500, 500])
            
            power_demand['low'] = fuzz.trimf(power_demand.universe, [0, 0, 50])
            power_demand['medium'] = fuzz.trimf(power_demand.universe, [30, 100, 170])
            power_demand['high'] = fuzz.trimf(power_demand.universe, [150, 200, 200])
            
            solar_activity['low'] = fuzz.trimf(solar_activity.universe, [0, 0, 80])
            solar_activity['normal'] = fuzz.trimf(solar_activity.universe, [60, 100, 140])
            solar_activity['high'] = fuzz.trimf(solar_activity.universe, [120, 200, 300])
            
            system_health['poor'] = fuzz.trimf(system_health.universe, [0, 0, 40])
            system_health['good'] = fuzz.trimf(system_health.universe, [30, 70, 90])
            system_health['excellent'] = fuzz.trimf(system_health.universe, [80, 100, 100])
            
            energy_boost['reduce'] = fuzz.trimf(energy_boost.universe, [50, 50, 80])
            energy_boost['normal'] = fuzz.trimf(energy_boost.universe, [70, 100, 130])
            energy_boost['increase'] = fuzz.trimf(energy_boost.universe, [120, 150, 150])
            
            cooling_power['minimal'] = fuzz.trimf(cooling_power.universe, [0, 0, 40])
            cooling_power['moderate'] = fuzz.trimf(cooling_power.universe, [30, 100, 170])
            cooling_power['maximum'] = fuzz.trimf(cooling_power.universe, [150, 200, 200])
            
            orbit_adjustment['closer'] = fuzz.trimf(orbit_adjustment.universe, [-10, -5, 0])
            orbit_adjustment['maintain'] = fuzz.trimf(orbit_adjustment.universe, [-2, 0, 2])
            orbit_adjustment['farther'] = fuzz.trimf(orbit_adjustment.universe, [0, 5, 10])
            
            rules = [
                ctrl.Rule(temperature_C['critical'], 
                         (cooling_power['maximum'], energy_boost['reduce'], orbit_adjustment['farther'])),
ล
                ctrl.Rule(temperature_C['hot'], 
                         (cooling_power['moderate'], energy_boost['reduce'])),
                
                ctrl.Rule(temperature_C['cold'], 
                         (cooling_power['minimal'], energy_boost['increase'])),

                ctrl.Rule(panel_efficiency['poor'] & power_demand['high'], 
                         (energy_boost['increase'], orbit_adjustment['closer'])),
                
                ctrl.Rule(panel_efficiency['excellent'] & solar_activity['high'], 
                         (energy_boost['increase'], cooling_power['moderate'])),

                ctrl.Rule(system_health['poor'], 
                         (energy_boost['reduce'], cooling_power['moderate'])),

                ctrl.Rule(temperature_C['optimal'] & panel_efficiency['fair'] & system_health['good'], 
                         (energy_boost['normal'], cooling_power['moderate'], orbit_adjustment['maintain'])),

                ctrl.Rule(power_demand['high'] & system_health['excellent'], 
                         (energy_boost['increase'], cooling_power['moderate'])),
                
                ctrl.Rule(solar_activity['high'] & system_health['good'], 
                         (energy_boost['increase'], cooling_power['maximum']))
            ]
            
            self.fuzzy_system = ctrl.ControlSystem(rules)
            self.fuzzy_simulation = ctrl.ControlSystemSimulation(self.fuzzy_system)
            print(f"Fuzzy system initialized with {len(rules)} rules")
            
        except Exception as e:
            print(f"Error setting up fuzzy system: {e}")
            self.fuzzy_system = None
            self.fuzzy_simulation = None
    
    def start_cpp_simulations(self):
        print("Starting C++ simulations...")

        if not os.path.exists('solar_calculations.cpp'):
            print("Warning: solar_calculations.cpp not found")
            return False
            
        if not os.path.exists('dyson_sphere.cpp'):
            print("Warning: dyson_sphere.cpp not found") 
            return False
        
        try:
            print("Compiling solar calculator...")
            result = subprocess.run(['g++', '-std=c++17', '-O2', '-o', 'solar_calculator', 
                                   'solar_calculations.cpp'], 
                                  capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Solar calculator compilation failed: {result.stderr}")
                return False
                
            print("Compiling Dyson simulation...")
            result = subprocess.run(['g++', '-std=c++17', '-O2', '-o', 'dyson_simulation',
                                   'dyson_sphere.cpp'], 
                                  capture_output=True, text=True)
            if result.returncode != 0:
                print(f"Dyson simulation compilation failed: {result.stderr}")
                return False
                
            print("Starting simulations...")

            solar_proc = subprocess.Popen(['./solar_calculator', '1.0', '1.0'],
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True)
            self.cpp_processes.append(solar_proc)
            
            dyson_proc = subprocess.Popen(['./dyson_simulation', '1.0', '0.01', '500'],
                                        stdout=subprocess.PIPE, 
                                        stderr=subprocess.PIPE,
                                        universal_newlines=True)
            self.cpp_processes.append(dyson_proc)
            
            time.sleep(2)
            return True
            
        except Exception as e:
            print(f"Failed to start simulations: {e}")
            return False
    
    def read_and_process_data(self) -> Dict:
        data = {}

        try:
            with open('solar_realtime.json', 'r') as f:
                data['solar'] = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            data['solar'] = {
                'flux_Wm2': 1367,
                'radiation_pressure_Pa': 4.56e-6,
                'timestamp': time.time(),
                'solar_flare_active': 0.0,
                'flux_multiplier': 1.0
            }
        
        try:
            with open('dyson_status.json', 'r') as f:
                data['dyson'] = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            data['dyson'] = {
                'total_power_absorbed_W': 1e12,  # 1 TW
                'thermodynamic_efficiency': 0.3,
                'avg_exterior_temp_K': 350,
                'operational_components': 500,
                'num_components': 1000,
                'simulation_time_years': 0.0,
                'avg_degradation': 0.05
            }
        
        current_power = data['dyson'].get('total_power_absorbed_W', 0)
        dt = 1.0 
        self.total_energy_collected += current_power * dt  # Joules
        self.energy_collection_rate = current_power
        self.peak_power = max(self.peak_power, current_power)

        self.display_metrics.update({
            'power_output_MW': current_power / 1e6,  # Convert to MW
            'total_energy_GWh': self.total_energy_collected / 3.6e12,  # Convert to GWh
            'panels_active': data['dyson'].get('operational_components', 0),
            'panels_total': data['dyson'].get('num_components', 1000),
            'avg_temperature_C': data['dyson'].get('avg_exterior_temp_K', 300) - 273.15,
            'system_efficiency_percent': data['dyson'].get('thermodynamic_efficiency', 0) * 100,
            'time_running_hours': (time.time() - self.start_time) / 3600
        })
        
        return data
    
    def compute_fuzzy_control(self, data: Dict) -> Dict:
        if not FUZZY_AVAILABLE or self.fuzzy_simulation is None:
            return self.compute_rule_based_control(data)
        
        try:
            dyson = data.get('dyson', {})
            solar = data.get('solar', {})
            
            panel_eff = (dyson.get('operational_components', 500) / 
                        max(1, dyson.get('num_components', 1000))) * 100
            
            temp_C = dyson.get('avg_exterior_temp_K', 350) - 273.15
            
            power_demand_pct = min(200, (self.energy_collection_rate / 1e12) * 100)  # % of 1TW
            
            solar_activity_pct = solar.get('flux_multiplier', 1.0) * 100
            
            system_health_pct = (1.0 - dyson.get('avg_degradation', 0.1)) * 100
            
            self.fuzzy_simulation.input['panel_efficiency'] = max(0, min(100, panel_eff))
            self.fuzzy_simulation.input['temperature_C'] = max(-50, min(500, temp_C))
            self.fuzzy_simulation.input['power_demand'] = max(0, min(200, power_demand_pct))
            self.fuzzy_simulation.input['solar_activity'] = max(0, min(300, solar_activity_pct))
            self.fuzzy_simulation.input['system_health'] = max(0, min(100, system_health_pct))
            
            self.fuzzy_simulation.compute()
            
            controls = {
                'energy_boost': self.fuzzy_simulation.output.get('energy_boost', 100) / 100.0,
                'cooling_power': self.fuzzy_simulation.output.get('cooling_power', 100) / 100.0,
                'orbit_adjustment': self.fuzzy_simulation.output.get('orbit_adjustment', 0) / 100.0
            }
            
            self.last_control_outputs = controls
            
            return controls
            
        except Exception as e:
            print(f"Fuzzy computation error: {e}")
            return self.compute_rule_based_control(data)
    
    def compute_rule_based_control(self, data: Dict) -> Dict:
        dyson = data.get('dyson', {})
        
        temp_K = dyson.get('avg_exterior_temp_K', 350)
        efficiency = dyson.get('thermodynamic_efficiency', 0.3)
        degradation = dyson.get('avg_degradation', 0.1)
        
        controls = {
            'energy_boost': 1.0,
            'cooling_power': 1.0,
            'orbit_adjustment': 0.0
        }

        if temp_K > 450:
            controls['cooling_power'] = 2.0
            controls['energy_boost'] = 0.7
            controls['orbit_adjustment'] = 0.05 

        elif efficiency > 0.4 and degradation < 0.2:
            controls['energy_boost'] = 1.3
            controls['cooling_power'] = 1.2
        
        elif degradation > 0.5:
            controls['energy_boost'] = 0.8
            controls['cooling_power'] = 1.5
        
        self.last_control_outputs = controls
        return controls
    
    def apply_control_to_simulation(self, controls: Dict):
        if not controls:
            return
        
        control_data = {
            'timestamp': time.time(),
            'energy_multiplier': controls.get('energy_boost', 1.0),
            'cooling_multiplier': controls.get('cooling_power', 1.0),
            'orbit_delta': controls.get('orbit_adjustment', 0.0),
            'apply_controls': True
        }

        try:
            with open('control_commands.json', 'w') as f:
                json.dump(control_data, f, indent=2)
                
            with open('control_stream.txt', 'w') as f:
                f.write(f"{control_data['energy_multiplier']:.6f} ")
                f.write(f"{control_data['cooling_multiplier']:.6f} ")
                f.write(f"{control_data['orbit_delta']:.6f}\n")
                
        except Exception as e:
            print(f"Warning: Could not write control commands: {e}")
    
    def display_status(self, data: Dict, controls: Dict):
        metrics = self.display_metrics
        
        status_line = (
            f"\r🔋 {metrics['power_output_MW']:.1f}MW | "
            f"⚡ {metrics['total_energy_GWh']:.2f}GWh | "
            f"🌡️  {metrics['avg_temperature_C']:.0f}°C | "
            f"📊 {metrics['panels_active']}/{metrics['panels_total']} panels | "
            f"⭐ {metrics['system_efficiency_percent']:.1f}% | "
            f"🕒 {metrics['time_running_hours']:.1f}h | "
            f"🎛️  E:{controls.get('energy_boost', 1.0):.2f} "
            f"C:{controls.get('cooling_power', 1.0):.2f} "
            f"O:{controls.get('orbit_adjustment', 0.0):+.3f}"
        )

        solar = data.get('solar', {})
        if solar.get('solar_flare_active', 0) > 0.5:
            status_line += " 🔆FLARE"
        if solar.get('flux_multiplier', 1.0) > 1.1:
            status_line += f" ⚡{solar.get('flux_multiplier', 1.0):.1f}x"
        
        status_line += "                    "  # Clear end of line
        
        print(status_line, end="", flush=True)
    
    def save_final_report(self):
        report = {
            'simulation_summary': {
                'total_runtime_hours': (time.time() - self.start_time) / 3600,
                'total_energy_collected_GWh': self.total_energy_collected / 3.6e12,
                'peak_power_MW': self.peak_power / 1e6,
                'average_power_MW': (self.total_energy_collected / (time.time() - self.start_time)) / 1e6,
                'control_iterations': len(self.control_history)
            },
            'final_metrics': self.display_metrics,
            'control_history': self.control_history[-100:] if self.control_history else []
        }
        
        try:
            with open('dyson_final_report.json', 'w') as f:
                json.dump(report, f, indent=2)
            print(f"\n📊 Final report saved: dyson_final_report.json")
            print(f"⚡ Total energy collected: {report['simulation_summary']['total_energy_collected_GWh']:.2f} GWh")
            print(f"🔋 Peak power: {report['simulation_summary']['peak_power_MW']:.1f} MW")
        except Exception as e:
            print(f"Could not save report: {e}")
    
    def run(self):
        print("IMPROVED DYSON SWARM FUZZY CONTROLLER")
        print("=" * 50)
        print("Enhanced with: Real control, Energy tracking, Clear metrics")
        print()
        
        sim_started = self.start_cpp_simulations()
        if not sim_started:
            print("Running in demo mode with simulated data...")
        
        print("Fuzzy control active. Press Ctrl+C to stop and save report...\n")
        
        iteration = 0
        
        while self.running:
            try:
                data = self.read_and_process_data()

                controls = self.compute_fuzzy_control(data)

                self.apply_control_to_simulation(controls)

                self.display_status(data, controls)

                self.control_history.append({
                    'iteration': iteration,
                    'time_hours': (time.time() - self.start_time) / 3600,
                    'metrics': self.display_metrics.copy(),
                    'controls': controls.copy()
                })
                
                if len(self.control_history) > 1000:
                    self.control_history = self.control_history[-500:]
                
                iteration += 1
                time.sleep(1.0)  
                
            except KeyboardInterrupt:
                break
            except Exception as e:
                print(f"\nControl loop error: {e}")
                time.sleep(0.5)
        
        self.save_final_report()
        self.cleanup()

def main():
    controller = ImprovedDysonController()
    controller.run()

if __name__ == "__main__":
    main()