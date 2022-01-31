#include <iostream>
#include <string>
#include <vector>
#include "runtime_config.h"
#include "file_io.h"

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: ./sim <use_cuda> <config_file> <time_step> <total_steps> <output_file>\n";
        return 1;
    }

    bool useCuda = std::stoi(argv[1]);
    std::string configFile = argv[2];
    double timeStep = std::stod(argv[3]);
    int totalSteps = std::stoi(argv[4]);
    std::string outputFile = argv[5];

    int referenceBody = 0;
    auto bodies = parseConfig(configFile, referenceBody);
    std::vector<SimulationResult> results;

    try {
        auto simulator = getSimImplementation(useCuda);
        simulator->runSim(bodies, timeStep, totalSteps, results);
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    writeCSV(results, outputFile);
    return 0;
}
