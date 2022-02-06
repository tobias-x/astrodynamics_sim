#pragma once
#include <vector>
#include <string>

struct Body {
    std::string name;
    double x,y,z;
    double vx,vy,vz;
    double mass;
    double diameter;
};

struct GridBody {
    double initial_x, initial_y, initial_z;
    double x, y, z;
    int i, j; // <<< NEW
};


struct SimulationResult {
    int timestep;
    int bodyIndex;
    double x, y, z, vx, vy, vz;
};
