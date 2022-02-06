#pragma once
#include "types.h"
#include <vector>
#include <atomic>

class SimulationEngine {
public:
    SimulationEngine(std::vector<Body> initialBodies, std::vector<GridBody> initialGridBodies, double dt);
    void run();
    void stop();
    const std::vector<Body>& latest() const;
    const std::vector<GridBody>& latestGrid() const;

private:
    std::vector<Body> buffers[2];
    std::vector<GridBody> gridBuffers[2];
    std::atomic<int> currentRead{0};
    std::atomic<bool> running{true};
    double dt;
};
