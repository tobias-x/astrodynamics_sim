#include "file_io.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

// Constants for unit conversion
constexpr double KM_TO_M = 1000.0;             // Kilometers to meters
constexpr double SOLAR_MASS_IN_KG = 1.9891e30; // Solar masses to kilograms

// Parse configuration file to initialize bodies
std::vector<Body> parseConfig(const std::string& filename, int& referenceBody) {
    std::vector<Body> bodies;
    std::ifstream infile(filename);
    std::string line;

    if (!infile.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return bodies;
    }

    // Read the first line for the reference body index
    if (std::getline(infile, line)) {
        size_t commaPos = line.find(',');
        if (commaPos != std::string::npos) {
            try {
                referenceBody = std::stoi(line.substr(commaPos + 1));
            } catch (...) {
                std::cerr << "Error parsing reference body index.\n";
                return bodies;
            }
        } else {
            std::cerr << "Error: Missing reference body index.\n";
            return bodies;
        }
    }

    // Skip the header line
    if (!std::getline(infile, line)) {
        std::cerr << "Error: Missing header line in configuration file.\n";
        return bodies;
    }

    // Parse the bodies
    while (std::getline(infile, line)) {
        Body body;
        std::istringstream iss(line);

        std::string name;
        char delimiter;

        // Parse line in the format: name,x,y,z,vx,vy,vz,mass,diameter
        std::getline(iss, name, ',');
        if (!(iss >> body.x >> delimiter
                  >> body.y >> delimiter
                  >> body.z >> delimiter
                  >> body.vx >> delimiter
                  >> body.vy >> delimiter
                  >> body.vz >> delimiter
                  >> body.mass >> delimiter
                  >> body.diameter)) {
            std::cerr << "Error parsing line: " << line << std::endl;
            continue;
        }

        // Convert units
        body.x *= KM_TO_M;           // Position from km to m
        body.y *= KM_TO_M;
        body.z *= KM_TO_M;

        body.vx *= KM_TO_M;          // Velocity from km/s to m/s
        body.vy *= KM_TO_M;
        body.vz *= KM_TO_M;

        body.mass *= SOLAR_MASS_IN_KG; // Mass from solar masses to kg

        body.diameter *= KM_TO_M;    // Diameter from km to m

        bodies.push_back(body);
    }

    return bodies;
}

// Write simulation results to a CSV file
void writeCSV(const std::vector<SimulationResult>& results, const std::string& filename) {
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << " for writing." << std::endl;
        return;
    }

    // Write header
    outfile << "timestep,bodyIndex,x,y,z,vx,vy,vz\n";

    // Write results
    for (const auto& result : results) {
        outfile << result.timestep << ","
                << result.bodyIndex << ","
                << result.x << ","
                << result.y << ","
                << result.z << ","
                << result.vx << ","
                << result.vy << ","
                << result.vz << "\n";
    }

    std::cout << "Successfully wrote simulation results to " << filename << std::endl;
}
