#include "run_realtime_sim.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: ./runtime_sim[(_gpu)] <config_file> <time_step> <duration> <output_file>\n";
        return 1;
    }

    std::string configFile = argv[1];
    double timeStep = std::stod(argv[2]);
    double duration = std::stod(argv[3]);
    std::string outputFile = argv[4];

#ifdef USE_CUDA
    bool useCuda = true;
#else
    bool useCuda = false;
#endif

    runRealtimeSim(configFile, timeStep, duration, useCuda, outputFile);
    return 0;
}
