#include "simulation_engine.h"
#include "cpu_sim.h"
#include <iostream>

#include <thread>
#include <chrono>

SimulationEngine::SimulationEngine(std::vector<Body> initialBodies, std::vector<GridBody> initialGridBodies, double dt)
    : buffers{ std::move(initialBodies), {} },
      gridBuffers{ std::move(initialGridBodies), {} },
      dt(dt)
{
    buffers[1] = buffers[0];

    for (const auto& body : buffers[0]) {
        std::cout << body.name << std::endl;
    }

    gridBuffers[1] = gridBuffers[0];
}

void SimulationEngine::run() {
    CpuSim sim;
    while (running.load(std::memory_order_acquire)) {
        int r = currentRead.load(std::memory_order_acquire);
        int w = 1 - r;
        buffers[w] = buffers[r];
        gridBuffers[w] = gridBuffers[r];
        sim.step(buffers[w], dt);
        sim.calcGrid(gridBuffers[w], buffers[w]);
        currentRead.store(w, std::memory_order_release);
    }
}

void SimulationEngine::stop() {
    running.store(false, std::memory_order_release);
}

const std::vector<Body>& SimulationEngine::latest() const {
    return buffers[currentRead.load(std::memory_order_acquire)];
}

const std::vector<GridBody>& SimulationEngine::latestGrid() const {
    return gridBuffers[currentRead.load(std::memory_order_acquire)];
}
