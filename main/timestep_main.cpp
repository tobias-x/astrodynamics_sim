#include "run_timestep_sim.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: ./timestep_sim[(_gpu)] <config_file> <time_step> <total_steps> <output_file>\n";
        return 1;
    }

    std::string configFile = argv[1];
    double timeStep = std::stod(argv[2]);
    int totalSteps = std::stoi(argv[3]);
    std::string outputFile = argv[4];

#ifdef USE_CUDA
    bool useCuda = true;
#else
    bool useCuda = false;
#endif

    runTimestepSim(configFile, timeStep, totalSteps, useCuda, outputFile);
    return 0;
}
