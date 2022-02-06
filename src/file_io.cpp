#include "file_io.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

static constexpr double KM_TO_M = 1000.0;
static constexpr double SOLAR_MASS_IN_KG = 1.9891e30;

std::vector<Body> parseConfig(const std::string& filename) {
    std::vector<Body> bodies;
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Cannot open " << filename << "\n";
        return bodies;
    }

    std::string line;

    // Skip the header line
    if (!std::getline(in, line)) {
        std::cerr << "Empty config file!\n";
        return bodies;
    }

    while (std::getline(in, line)) {
        if (line.empty()) continue; // Skip empty lines

        std::istringstream iss(line);
        Body b;
        char comma;

        // Read name
        if (!std::getline(iss, b.name, ',')) {
            std::cerr << "Warning: Could not parse name from line:\n" << line << "\n";
            continue;
        }

        // Read all other fields
        if (!(iss >> b.x >> comma >> b.y >> comma >> b.z >> comma 
                  >> b.vx >> comma >> b.vy >> comma >> b.vz >> comma 
                  >> b.mass >> comma >> b.diameter)) {
            std::cerr << "Warning: Skipping badly formatted line:\n" << line << "\n";
            continue;
        }

        // Unit conversions
        b.x *= KM_TO_M;
        b.y *= KM_TO_M;
        b.z *= KM_TO_M;
        b.vx *= KM_TO_M;
        b.vy *= KM_TO_M;
        b.vz *= KM_TO_M;
        b.mass *= SOLAR_MASS_IN_KG;
        b.diameter *= KM_TO_M; // optional if you want meters for diameter

        bodies.push_back(b);
    }

    std::cout << "Loaded " << bodies.size() << " bodies from " << filename << "\n";
    return bodies;
}

void writeCSV(const std::vector<SimulationResult>& results, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Cannot open output file " << filename << "\n";
        return;
    }

    out << "timestep,bodyIndex,x,y,z,vx,vy,vz\n";
    for (const auto& r : results) {
        out << r.timestep << "," << r.bodyIndex << ","
            << r.x << "," << r.y << "," << r.z << ","
            << r.vx << "," << r.vy << "," << r.vz << "\n";
    }
    std::cout << "Wrote " << results.size() << " results to " << filename << "\n";
}
