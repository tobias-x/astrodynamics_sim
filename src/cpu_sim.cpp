#include "cpu_sim.h"
#include <cmath>
#include <iostream>
#include <vector>

void CpuSim::runSim(const std::vector<Body>& bodies, double timeStep, int totalSteps, std::vector<SimulationResult>& results) const {
    auto localBodies = bodies;

    for (int step = 0; step < totalSteps; ++step) {
        std::vector<Body> updatedBodies = localBodies;

        for (size_t i = 0; i < localBodies.size(); ++i) {
            double ax = 0.0, ay = 0.0, az = 0.0;

            for (size_t j = 0; j < localBodies.size(); ++j) {
                if (i == j) continue;

                // Gravitational force calculation
                double dx = localBodies[j].x - localBodies[i].x;
                double dy = localBodies[j].y - localBodies[i].y;
                double dz = localBodies[j].z - localBodies[i].z;

                double distSquared = dx * dx + dy * dy + dz * dz + 1e-5; // Softening factor
                double dist = std::sqrt(distSquared);

                double force = G * localBodies[j].mass / distSquared;
                ax += force * dx / dist;
                ay += force * dy / dist;
                az += force * dz / dist;
            }

            // Update velocities
            updatedBodies[i].vx += ax * timeStep;
            updatedBodies[i].vy += ay * timeStep;
            updatedBodies[i].vz += az * timeStep;

            // Update positions
            updatedBodies[i].x += updatedBodies[i].vx * timeStep;
            updatedBodies[i].y += updatedBodies[i].vy * timeStep;
            updatedBodies[i].z += updatedBodies[i].vz * timeStep;

            results.push_back({
                step,
                static_cast<int>(i),
                updatedBodies[i].x,
                updatedBodies[i].y,
                updatedBodies[i].z,
                updatedBodies[i].vx,
                updatedBodies[i].vy,
                updatedBodies[i].vz
            });
        }

        localBodies = std::move(updatedBodies);
    }
}
