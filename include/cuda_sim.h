#pragma once
#include <vector>
#include "file_io.h"

#ifdef USE_CUDA
class CudaSim {
public:
    void step(std::vector<Body>& bodies, double dt) const;
    void runSim(const std::vector<Body>& bodies,
                double timeStep,
                int totalSteps,
                std::vector<SimulationResult>& results) const;
};
#else
class CudaSim {
public:
    void step(std::vector<Body>&, double) const {}
    void runSim(const std::vector<Body>&, double, int, std::vector<SimulationResult>&) const {}
};
#endif
