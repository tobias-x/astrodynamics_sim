#include "cpu_sim.h"  // Now includes the header instead of cpp
#include "file_io.h"
#include <iostream>
#include <vector>
#include <string>

void runTimestepSim(const std::string& configFile, double timeStep, int totalSteps, bool useCuda, const std::string& outputFile) {
    int referenceBody = 0;

    // Parse the configuration file to initialize bodies
    auto bodies = parseConfig(configFile, referenceBody);
    if (bodies.empty()) {
        std::cerr << "Error: No bodies found in the configuration file.\n";
        return;
    }

    std::cout << "Starting timestep simulation with " << bodies.size() << " bodies...\n";

    // Use CPU-based simulation
    CpuSim simulator;

    // Store results
    std::vector<SimulationResult> results;

    // Run the simulation
    simulator.runSim(bodies, timeStep, totalSteps, results);

    // Write results to the output file
    writeCSV(results, outputFile);

    std::cout << "Timestep simulation completed. Results written to " << outputFile << "\n";
}
