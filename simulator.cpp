// g++ -o simulator simulator.cpp -lGL -lGLU -lglut -std=c++17
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <thread>
#include <signal.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <cstdarg>

namespace Constants {
    constexpr double G = 6.67430e-11;
    constexpr double c = 2.99792458e8;
    constexpr double AU = 1.495978707e11;
    constexpr double PI = 3.14159265359;
    constexpr double SOLAR_MASS = 1.989e30;
    constexpr double SOLAR_RADIUS = 6.96e8;
    constexpr double SOLAR_LUMINOSITY = 3.828e26;
}

struct Vec3 {
    double x, y, z;
    Vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
    double length() const { return sqrt(x*x + y*y + z*z); }
    Vec3 normalize() const { 
        double len = length();
        return len > 0 ? Vec3(x/len, y/len, z/len) : Vec3(0, 0, 0);
    }
    double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
};

struct SolarData {
    double flux_Wm2 = 1367.0;
    double radiation_pressure_Pa = 4.56e-6;
    double solar_flare_active = 0.0;
    double flux_multiplier = 1.0;
    double degradation_factor_per_year = 0.05;
    double equilibrium_temperature_K = 300.0;
};

struct DysonData {
    double total_power_absorbed_W = 1e12;
    double thermodynamic_efficiency = 0.3;
    double avg_exterior_temp_K = 350.0;
    int operational_components = 500;
    int num_components = 1000;
    double simulation_time_years = 0.0;
    double avg_degradation = 0.05;
};

class DysonPanel {
public:
    Vec3 position, velocity, normal;
    double area, efficiency, temperature;
    double energy_collected, total_energy_lifetime;
    bool active;
    double orbital_radius, orbital_angle, orbital_speed;
    double orbital_inclination, orbital_phase;
    
    DysonPanel(Vec3 pos, double panel_area = 1000.0) 
        : position(pos), area(panel_area), efficiency(0.25), temperature(300.0),
          energy_collected(0.0), total_energy_lifetime(0.0), active(true) {
        
        orbital_radius = position.length();
        orbital_angle = atan2(position.y, position.x);
        orbital_inclination = atan2(sqrt(position.x*position.x + position.y*position.y), position.z);
        orbital_phase = 0.0;
        orbital_speed = sqrt(Constants::G * Constants::SOLAR_MASS / orbital_radius);
        update_orientation();
    }
    
    void update_orientation() {
        Vec3 to_sun = (Vec3(0,0,0) - position).normalize();
        normal = to_sun;
    }
    
    void update_orbit(double dt, const DysonData& dyson_data) {
        double angular_velocity = orbital_speed / orbital_radius;
        orbital_phase += angular_velocity * dt;
        
        if (orbital_phase > 2.0 * Constants::PI) {
            orbital_phase -= 2.0 * Constants::PI;
        }
        
        position.x = orbital_radius * sin(orbital_inclination) * cos(orbital_angle + orbital_phase);
        position.y = orbital_radius * sin(orbital_inclination) * sin(orbital_angle + orbital_phase);
        position.z = orbital_radius * cos(orbital_inclination);
        
        velocity.x = -orbital_speed * sin(orbital_inclination) * sin(orbital_angle + orbital_phase);
        velocity.y = orbital_speed * sin(orbital_inclination) * cos(orbital_angle + orbital_phase);
        velocity.z = -orbital_speed * cos(orbital_inclination) * sin(orbital_angle + orbital_phase);
        
        update_orientation();
    }
    
    double collect_energy(double dt, const SolarData& solar_data) {
        if (!active) return 0.0;
        
        double distance = position.length();
        if (distance < 1e-6) return 0.0;
        
        double flux = solar_data.flux_Wm2;
        Vec3 to_sun = (Vec3(0,0,0) - position).normalize();
        double cos_angle = std::max(0.0, normal.dot(to_sun));
        
        double power_collected = flux * area * efficiency * cos_angle;
        double energy_this_step = power_collected * dt;
        
        energy_collected += energy_this_step;
        total_energy_lifetime += energy_this_step;
        update_temperature(flux * cos_angle, dt, solar_data);
        
        return energy_this_step;
    }
    
private:
    void update_temperature(double incident_flux, double dt, const SolarData& solar_data) {
        double absorbed_power = incident_flux * area * 0.9;
        double radiated_power = 5.67e-8 * area * pow(temperature, 4);
        double net_power = absorbed_power - radiated_power;
        
        double thermal_mass = area * 50.0;
        double heat_capacity = 900.0;
        double dT_dt = net_power / (thermal_mass * heat_capacity);
        
        temperature += dT_dt * dt;
        temperature = std::max(200.0, std::min(800.0, temperature));
        
        if (temperature > 350.0) {
            efficiency = 0.25 * (1.0 - 0.001 * (temperature - 350.0));
            efficiency = std::max(0.05, efficiency);
        }
        
        if (temperature > 700.0) active = false;
    }
};

struct EnergyParticle {
    Vec3 position, velocity;
    float color[3];
    double lifetime, energy;
    
    EnergyParticle(Vec3 pos, Vec3 vel, double e) 
        : position(pos), velocity(vel), lifetime(0.0), energy(e) {
        color[0] = 1.0f; color[1] = 0.8f; color[2] = 0.2f;
    }
};

class DysonSimulator {
private:
    std::vector<DysonPanel> panels;
    std::vector<EnergyParticle> particles;
    std::mt19937 rng;
    
    double sim_time, time_step, total_energy_collected;
    double power_generation_rate, completion_percentage;
    
    double camera_distance, camera_angle_x, camera_angle_y;
    bool show_orbits, show_energy_flow, auto_rotate_camera;
    int window_width, window_height;
    
    SolarData solar_data;
    DysonData dyson_data;
    
public:
    DysonSimulator() : rng(std::chrono::steady_clock::now().time_since_epoch().count()),
                      sim_time(0.0), time_step(0.1), total_energy_collected(0.0),
                      power_generation_rate(0.0), completion_percentage(0.0),
                      camera_distance(8.0), camera_angle_x(20.0), camera_angle_y(0.0),
                      show_orbits(true), show_energy_flow(true), auto_rotate_camera(true),
                      window_width(1200), window_height(900) {
        initialize_dyson_swarm_realistic();
    }

    void set_window_size(int w, int h) { window_width = w; window_height = h; }

    void initialize_dyson_swarm_realistic() {
        panels.clear();
        
        int num_panels = 150;
        std::vector<double> shell_radii = {2.5, 2.8, 3.0, 3.2, 3.5};
        
        std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * Constants::PI);
        std::uniform_real_distribution<double> inclination_dist(0.0, Constants::PI);
        std::uniform_real_distribution<double> area_dist(800.0, 1200.0);
        
        for (size_t shell_idx = 0; shell_idx < shell_radii.size(); shell_idx++) {
            double radius = shell_radii[shell_idx];
            int panels_in_shell = num_panels / shell_radii.size();
            
            for (int i = 0; i < panels_in_shell; i++) {
                double theta = angle_dist(rng);
                double phi = inclination_dist(rng);
                
                Vec3 pos(
                    radius * sin(phi) * cos(theta),
                    radius * sin(phi) * sin(theta),
                    radius * cos(phi)
                );
                
                double area = area_dist(rng);
                panels.emplace_back(pos, area);
            }
        }
    }
    
    void initialize_dyson_sphere() {
        panels.clear();
        
        int num_panels = 200;
        double sphere_radius = 3.0;
        
        std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * Constants::PI);
        std::uniform_real_distribution<double> inclination_dist(0.0, Constants::PI);
        std::uniform_real_distribution<double> area_dist(800.0, 1200.0);
        
        for (int i = 0; i < num_panels; i++) {
            double theta = angle_dist(rng);
            double phi = inclination_dist(rng);
            
            Vec3 pos(
                sin(phi) * cos(theta) * sphere_radius,
                sin(phi) * sin(theta) * sphere_radius,
                cos(phi) * sphere_radius
            );
            
            double area = area_dist(rng);
            panels.emplace_back(pos, area);
        }
    }
    
    void load_physics_data() {
        std::ifstream solar_file("solar_realtime.json");
        if (solar_file.is_open()) {
            std::string line;
            while (std::getline(solar_file, line)) {
                if (line.find("flux_Wm2") != std::string::npos) {
                    size_t pos = line.find(": ");
                    if (pos != std::string::npos) {
                        solar_data.flux_Wm2 = std::stod(line.substr(pos + 2));
                    }
                }
                else if (line.find("solar_flare_active") != std::string::npos) {
                    size_t pos = line.find(": ");
                    if (pos != std::string::npos) {
                        solar_data.solar_flare_active = std::stod(line.substr(pos + 2));
                    }
                }
                else if (line.find("flux_multiplier") != std::string::npos) {
                    size_t pos = line.find(": ");
                    if (pos != std::string::npos) {
                        solar_data.flux_multiplier = std::stod(line.substr(pos + 2));
                    }
                }
            }
            solar_file.close();
        }
        
        std::ifstream dyson_file("dyson_status.json");
        if (dyson_file.is_open()) {
            std::string line;
            while (std::getline(dyson_file, line)) {
                if (line.find("total_power_absorbed_W") != std::string::npos) {
                    size_t pos = line.find(": ");
                    if (pos != std::string::npos) {
                        dyson_data.total_power_absorbed_W = std::stod(line.substr(pos + 2));
                    }
                }
                else if (line.find("thermodynamic_efficiency") != std::string::npos) {
                    size_t pos = line.find(": ");
                    if (pos != std::string::npos) {
                        dyson_data.thermodynamic_efficiency = std::stod(line.substr(pos + 2));
                    }
                }
                else if (line.find("avg_exterior_temp_K") != std::string::npos) {
                    size_t pos = line.find(": ");
                    if (pos != std::string::npos) {
                        dyson_data.avg_exterior_temp_K = std::stod(line.substr(pos + 2));
                    }
                }
                else if (line.find("operational_components") != std::string::npos) {
                    size_t pos = line.find(": ");
                    if (pos != std::string::npos) {
                        dyson_data.operational_components = std::stoi(line.substr(pos + 2));
                    }
                }
            }
            dyson_file.close();
        }
    }
    
    void update_simulation(double dt) {
        sim_time += dt;
        load_physics_data();
        
        double energy_this_step = 0.0;
        
        for (auto& panel : panels) {
            panel.update_orbit(dt, dyson_data);
            energy_this_step += panel.collect_energy(dt, solar_data);
        }
        
        total_energy_collected += energy_this_step;
        power_generation_rate = energy_this_step / dt;
        
        double theoretical_max = Constants::SOLAR_LUMINOSITY * dt;
        completion_percentage = std::min(100.0, (energy_this_step / theoretical_max) * 100.0);
        
        if (show_energy_flow) {
            update_energy_particles(dt);
        }
    }
    
    void update_energy_particles(double dt) {
        particles.erase(
            std::remove_if(particles.begin(), particles.end(),
                [](const EnergyParticle& p) { return p.lifetime > 5.0 || p.position.length() > 12.0; }),
            particles.end()
        );
        
        if (particles.size() < 200 && (rng() % 100) < 40) {
            double theta = 2.0 * Constants::PI * (rng() % 1000) / 1000.0;
            double phi = Constants::PI * (rng() % 1000) / 1000.0;
            
            Vec3 direction(sin(phi) * cos(theta), sin(phi) * sin(theta), cos(phi));
            Vec3 start_pos = direction * 0.2;
            Vec3 velocity = direction * 2.5;
            
            particles.emplace_back(start_pos, velocity, solar_data.flux_Wm2);
        }
        
        for (auto& particle : particles) {
            particle.position = particle.position + particle.velocity * dt;
            particle.lifetime += dt;
            
            float fade = std::max(0.0f, 1.0f - (float)particle.lifetime / 5.0f);
            particle.color[0] = fade;
            particle.color[1] = fade * 0.8f;
            particle.color[2] = fade * 0.2f;
        }
    }
    
    void render() {
        glViewport(0, 0, window_width, window_height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        render_ui();
        
        int ui_width = window_width * 0.3;
        glViewport(ui_width, 0, window_width - ui_width, window_height);
        
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(45.0, (double)(window_width - ui_width) / window_height, 0.1, 100.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_LIGHTING);
        
        glTranslatef(0, 0, -camera_distance);
        glRotatef(camera_angle_x, 1, 0, 0);
        glRotatef(camera_angle_y, 0, 1, 0);
        
        if (auto_rotate_camera) {
            camera_angle_y += 0.3;
            if (camera_angle_y > 360) camera_angle_y -= 360;
        }
        
        render_sun();
        render_dyson_panels();
        
        if (show_orbits) render_orbital_paths();
        if (show_energy_flow) render_energy_particles();
        
        glutSwapBuffers();
    }
    
    void render_sun() {
        glPushMatrix();
        
        if (solar_data.solar_flare_active > 0.5) {
            glColor3f(1.0f, 0.4f, 0.1f);
        } else {
            glColor3f(1.0f, 0.9f, 0.2f);
        }
        
        glutSolidSphere(0.15, 20, 20);
        
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(1.0f, 0.7f, 0.1f, 0.4f);
        glutSolidSphere(0.25, 15, 15);
        glDisable(GL_BLEND);
        
        glPopMatrix();
    }
    
    void render_dyson_panels() {
        for (const auto& panel : panels) {
            glPushMatrix();
            glTranslatef(panel.position.x, panel.position.y, panel.position.z);

            if (panel.active) {
                double efficiency = std::min(1.0, panel.energy_collected / 1e6);
                if (efficiency > 0.7) {
                    glColor3f(0.1f, 0.7f, 1.0f);
                } else if (efficiency > 0.4) {
                    glColor3f(0.2f, 1.0f, 0.3f);
                } else if (efficiency > 0.1) {
                    glColor3f(1.0f, 0.8f, 0.1f);
                } else {
                    glColor3f(1.0f, 0.4f, 0.1f);
                }
            } else {
                glColor3f(0.6f, 0.1f, 0.1f);
            }

            glutSolidCube(0.05);
            glPopMatrix();
        }
    }
    
    void render_orbital_paths() {
        std::vector<double> radii = {2.5, 2.8, 3.0, 3.2, 3.5};
        
        glColor3f(0.2f, 0.3f, 0.6f);
        glLineWidth(1.0f);
        
        for (double radius : radii) {
            glBegin(GL_LINE_LOOP);
            for (int i = 0; i < 64; i++) {
                double angle = 2.0 * Constants::PI * i / 64.0;
                glVertex3f(radius * cos(angle), radius * sin(angle), 0);
            }
            glEnd();
        }
    }

    void render_energy_particles() {
        glPointSize(4.0f);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE);

        glBegin(GL_POINTS);
        for (const auto& particle : particles) {
            float alpha = std::max(0.0f, 1.0f - (float)particle.lifetime / 5.0f);
            glColor4f(1.0f, 0.6f, 0.2f, alpha);
            glVertex3f(particle.position.x, particle.position.y, particle.position.z);
        }
        glEnd();

        glDisable(GL_BLEND);
    }

    void render_text_small(int x, int y, const std::string& text, float r = 1.0f, float g = 1.0f, float b = 1.0f) {
        glColor3f(r, g, b);
        glRasterPos2i(x, y);
        for (char c : text) {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, c);
        }
    }

    void render_text_medium(int x, int y, const std::string& text, float r = 1.0f, float g = 1.0f, float b = 1.0f) {
        glColor3f(r, g, b);
        glRasterPos2i(x, y);
        for (char c : text) {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, c);
        }
    }

    void render_ui() {
        int ui_width = window_width * 0.3;
        glViewport(0, 0, ui_width, window_height);
        
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, ui_width, window_height, 0, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(0.0f, 0.0f, 0.0f, 0.9f);
        glBegin(GL_QUADS);
        glVertex2i(0, 0);
        glVertex2i(ui_width, 0);
        glVertex2i(ui_width, window_height);
        glVertex2i(0, window_height);
        glEnd();
        glDisable(GL_BLEND);
        
        int line_y = 30;
        int line_height = 20;
        int margin = 10;
        
        render_text_medium(margin, line_y, "DYSON CONTROL", 0.2f, 1.0f, 0.4f);
        line_y += 40;
        
        render_text_small(margin, line_y, format_string("Runtime: %.1f min", sim_time / 60.0), 0.9f, 0.9f, 0.9f);
        line_y += line_height;
        
        double power_tw = power_generation_rate / 1e12;
        render_text_small(margin, line_y, format_string("Power: %.2f TW", power_tw), 0.2f, 1.0f, 0.2f);
        line_y += line_height;
        
        render_text_small(margin, line_y, format_string("Energy: %.1f PJ", total_energy_collected / 1e15), 0.9f, 0.9f, 0.9f);
        line_y += line_height;
        
        render_text_small(margin, line_y, format_string("Efficiency: %.1f%%", completion_percentage), 0.9f, 0.9f, 0.9f);
        line_y += line_height + 10;
        
        int active = count_active_panels();
        int total = panels.size();
        render_text_small(margin, line_y, format_string("Panels: %d / %d", active, total), 0.9f, 0.9f, 0.9f);
        line_y += line_height;
        
        double avg_temp = calculate_average_temperature();
        render_text_small(margin, line_y, format_string("Temp: %.0f K", avg_temp), 0.9f, 0.9f, 0.9f);
        line_y += line_height + 10;
        
        if (solar_data.solar_flare_active > 0.5) {
            render_text_small(margin, line_y, "SOLAR FLARE!", 1.0f, 0.2f, 0.2f);
        } else {
            render_text_small(margin, line_y, "Solar: Normal", 0.2f, 1.0f, 0.2f);
        }
        line_y += line_height + 15;
        
        render_text_small(margin, line_y, "CONTROLS:", 1.0f, 0.8f, 0.2f);
        line_y += 25;
        render_text_small(margin, line_y, "SPACE - Reset", 0.7f, 0.7f, 0.7f);
        line_y += line_height;
        render_text_small(margin, line_y, "S - Full Sphere", 0.7f, 0.7f, 0.7f);
        line_y += line_height;
        render_text_small(margin, line_y, "W - Swarm Mode", 0.7f, 0.7f, 0.7f);
        line_y += line_height;
        render_text_small(margin, line_y, "O - Orbits", 0.7f, 0.7f, 0.7f);
        line_y += line_height;
        render_text_small(margin, line_y, "E - Particles", 0.7f, 0.7f, 0.7f);
        line_y += line_height;
        render_text_small(margin, line_y, "A - Auto Camera", 0.7f, 0.7f, 0.7f);
    }
    
    std::string format_string(const char* format, ...) {
        char buffer[256];
        va_list args;
        va_start(args, format);
        vsnprintf(buffer, sizeof(buffer), format, args);
        va_end(args);
        return std::string(buffer);
    }
    
    int count_active_panels() const {
        return std::count_if(panels.begin(), panels.end(), 
                           [](const DysonPanel& p) { return p.active; });
    }
    
    double calculate_average_temperature() const {
        if (panels.empty()) return 0.0;
        double sum = 0.0;
        for (const auto& panel : panels) {
            sum += panel.temperature;
        }
        return sum / panels.size();
    }
    
    void reset_simulation() {
        sim_time = 0.0;
        total_energy_collected = 0.0;
        power_generation_rate = 0.0;
        completion_percentage = 0.0;
        particles.clear();
        
        for (auto& panel : panels) {
            panel.energy_collected = 0.0;
            panel.total_energy_lifetime = 0.0;
            panel.temperature = 300.0;
            panel.active = true;
            panel.efficiency = 0.25;
        }
    }
    
    void toggle_orbits() { show_orbits = !show_orbits; }
    void toggle_energy_flow() { show_energy_flow = !show_energy_flow; }
    void toggle_auto_camera() { auto_rotate_camera = !auto_rotate_camera; }
    
    void set_camera(double dist, double angle_x, double angle_y) {
        camera_distance = dist;
        camera_angle_x = angle_x;
        camera_angle_y = angle_y;
    }
};

DysonSimulator* sim = nullptr;

void display() {
    if (sim) sim->render();
}

void idle() {
    if (sim) {
        sim->update_simulation(0.016);
    }
    glutPostRedisplay();
    std::this_thread::sleep_for(std::chrono::milliseconds(16));
}

void keyboard(unsigned char key, int x, int y) {
    if (sim) {
        switch (key) {
            case ' ': sim->reset_simulation(); break;
            case 's': case 'S': sim->initialize_dyson_sphere(); break;
            case 'w': case 'W': sim->initialize_dyson_swarm_realistic(); break;
            case 'o': case 'O': sim->toggle_orbits(); break;
            case 'e': case 'E': sim->toggle_energy_flow(); break;
            case 'a': case 'A': sim->toggle_auto_camera(); break;
            case 27: exit(0); break;
        }
    }
}

int mouse_x = 0, mouse_y = 0;
bool mouse_pressed = false;

void mouse(int button, int state, int x, int y) {
    mouse_pressed = (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
    mouse_x = x; mouse_y = y;
}

void motion(int x, int y) {
    if (mouse_pressed && sim) {
        static double angle_x = 20, angle_y = 0;
        angle_x += (y - mouse_y) * 0.5;
        angle_y += (x - mouse_x) * 0.5;
        
        angle_x = std::max(-90.0, std::min(90.0, angle_x));
        sim->set_camera(8.0, angle_x, angle_y);
        
        mouse_x = x; mouse_y = y;
    }
}

void reshape(int width, int height) {
    if (sim) {
        sim->set_window_size(width, height);
    }
}

void init_opengl() {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    
    GLfloat light_pos[] = {0.0, 0.0, 0.0, 1.0};
    GLfloat light_ambient[] = {0.3, 0.3, 0.3, 1.0};
    GLfloat light_diffuse[] = {1.0, 1.0, 0.8, 1.0};
    
    glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    
    glClearColor(0.02, 0.02, 0.08, 1.0);
}

volatile sig_atomic_t keep_running = 1;
void signal_handler(int signal) {
    keep_running = 0;
    std::cout << "\nShutting down Optimized Dyson Simulator...\n";
    exit(0);
}

int main(int argc, char** argv) {
    signal(SIGINT, signal_handler);
    
    std::cout << "OPTIMIZED DYSON SPHERE SIMULATOR WITH PHYSICS INTEGRATION\n";
    std::cout << "========================================================\n";
    std::cout << "Features:\n";
    std::cout << "  • Reads data from solar_calculations.cpp output\n";
    std::cout << "  • Integrates with dyson_sphere.cpp physics\n";
    std::cout << "  • Compact UI (30% of screen)\n";
    std::cout << "  • Real-time physics data integration\n\n";
    
    std::cout << "Controls:\n";
    std::cout << "  SPACE:  Reset simulation\n";
    std::cout << "  S:      Full Dyson Sphere\n";
    std::cout << "  W:      Dyson Swarm\n";
    std::cout << "  O:      Toggle orbital paths\n";
    std::cout << "  E:      Toggle energy particles\n";
    std::cout << "  A:      Toggle auto camera\n";
    std::cout << "  Mouse:  Manual camera control\n";
    std::cout << "  ESC:    Exit\n\n";
    
    std::cout << "To run with full physics integration:\n";
    std::cout << "1. Compile and run: ./solar_calculator\n";
    std::cout << "2. Compile and run: ./dyson_simulation\n";
    std::cout << "3. Run this simulator to visualize integrated data\n\n";
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1200, 900);
    glutCreateWindow("Optimized Dyson Sphere Simulator - Physics Integrated");
    
    init_opengl();
    
    sim = new DysonSimulator();
    
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutReshapeFunc(reshape);
    
    std::cout << "Starting optimized simulation...\n";
    std::cout << "UI occupies 30% of screen width\n";
    std::cout << "Reading physics data from JSON files\n\n";
    
    glutMainLoop();
    
    delete sim;
    return 0;
}