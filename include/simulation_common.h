#ifndef SIMULATION_COMMON_H
#define SIMULATION_COMMON_H

#include <vector>

constexpr double G = 6.67430e-11;

// Body structure
struct Body {
    double x, y, z;       // Position (meters)
    double vx, vy, vz;    // Velocity (meters/second)
    double mass;          // Mass (kilograms)
    double diameter;      // Diameter (meters)
};

// Simulation result structure
struct SimulationResult {
    int timestep;
    int bodyIndex;
    double x, y, z;
    double vx, vy, vz;
};

#endif // SIMULATION_COMMON_H
