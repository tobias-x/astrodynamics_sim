#ifndef RUN_REALTIME_SIM_H
#define RUN_REALTIME_SIM_H

#include <string>

void runRealtimeSim(const std::string& configFile, double timeStep, double duration, bool useCuda, const std::string& outputFile);

#endif // RUN_REALTIME_SIM_H
