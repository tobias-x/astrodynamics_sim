#pragma once
#include <vector>
#include "file_io.h"
#include "types.h"

class CpuSim {
public:
    void step(std::vector<Body>& bodies, double dt) const;
    void calcGrid(std::vector<GridBody>& gridBodies, const std::vector<Body>& bodies) const;
    void runSim(const std::vector<Body>& bodies,
                double timeStep,
                int totalSteps,
                std::vector<SimulationResult>& results) const;
};
