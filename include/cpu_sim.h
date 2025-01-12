#ifndef CPU_SIM_H
#define CPU_SIM_H

#include "simulation_common.h"
#include <vector>

// CPU-based simulation class
class CpuSim {
public:
    // Runs the simulation
    void runSim(const std::vector<Body>& bodies, double timeStep, int totalSteps, std::vector<SimulationResult>& results) const;
};

#endif // CPU_SIM_H
