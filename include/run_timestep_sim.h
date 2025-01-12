#ifndef RUN_TIMESTEP_SIM_H
#define RUN_TIMESTEP_SIM_H

#include <string>

void runTimestepSim(const std::string& configFile, double timeStep, int totalSteps, bool useCuda, const std::string& outputFile);

#endif // RUN_TIMESTEP_SIM_H
