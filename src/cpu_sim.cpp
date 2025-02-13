#include "cpu_sim.h"
#include <cmath>
#include <iostream>
#include <vector>

void CpuSim::runSim(const std::vector<Body>& bodies, double timeStep, int totalSteps, std::vector<SimulationResult>& results) const {
    // Copy initial state
    std::vector<Body> localBodies = bodies;
    const double softening = 1e-5;

    // For each simulation step
    for (int step = 0; step < totalSteps; ++step) {
        // Temporary containers for RK4 derivatives
        std::vector<Body> newBodies(localBodies.size());
        std::vector<Body> k1(localBodies.size());
        std::vector<Body> k2(localBodies.size());
        std::vector<Body> k3(localBodies.size());
        std::vector<Body> k4(localBodies.size());
        std::vector<Body> tempBodies(localBodies.size());

        // --------------------------------------------
        // k1: Evaluate derivatives at t
        // For each body, dx/dt = vx and dv/dt = acceleration from localBodies.
        for (size_t i = 0; i < localBodies.size(); ++i) {
            k1[i].x = localBodies[i].vx;
            k1[i].y = localBodies[i].vy;
            k1[i].z = localBodies[i].vz;
            
            double ax = 0.0, ay = 0.0, az = 0.0;
            for (size_t j = 0; j < localBodies.size(); ++j) {
                if (i == j) continue;
                double dx = localBodies[j].x - localBodies[i].x;
                double dy = localBodies[j].y - localBodies[i].y;
                double dz = localBodies[j].z - localBodies[i].z;
                double distSq = dx*dx + dy*dy + dz*dz + softening;
                double dist = std::sqrt(distSq);
                double force = G * localBodies[j].mass / distSq;
                ax += force * dx / dist;
                ay += force * dy / dist;
                az += force * dz / dist;
            }
            k1[i].vx = ax;
            k1[i].vy = ay;
            k1[i].vz = az;
            k1[i].mass = localBodies[i].mass;
        }

        // --------------------------------------------
        // k2: Evaluate derivatives at t + dt/2 using state = localBodies + (dt/2)*k1
        for (size_t i = 0; i < localBodies.size(); ++i) {
            tempBodies[i].x  = localBodies[i].x  + (timeStep / 2.0) * k1[i].x;
            tempBodies[i].y  = localBodies[i].y  + (timeStep / 2.0) * k1[i].y;
            tempBodies[i].z  = localBodies[i].z  + (timeStep / 2.0) * k1[i].z;
            tempBodies[i].vx = localBodies[i].vx + (timeStep / 2.0) * k1[i].vx;
            tempBodies[i].vy = localBodies[i].vy + (timeStep / 2.0) * k1[i].vy;
            tempBodies[i].vz = localBodies[i].vz + (timeStep / 2.0) * k1[i].vz;
            tempBodies[i].mass = localBodies[i].mass;
        }
        for (size_t i = 0; i < localBodies.size(); ++i) {
            k2[i].x = tempBodies[i].vx;
            k2[i].y = tempBodies[i].vy;
            k2[i].z = tempBodies[i].vz;
            
            double ax = 0.0, ay = 0.0, az = 0.0;
            for (size_t j = 0; j < localBodies.size(); ++j) {
                if (i == j) continue;
                double dx = tempBodies[j].x - tempBodies[i].x;
                double dy = tempBodies[j].y - tempBodies[i].y;
                double dz = tempBodies[j].z - tempBodies[i].z;
                double distSq = dx*dx + dy*dy + dz*dz + softening;
                double dist = std::sqrt(distSq);
                double force = G * tempBodies[j].mass / distSq;
                ax += force * dx / dist;
                ay += force * dy / dist;
                az += force * dz / dist;
            }
            k2[i].vx = ax;
            k2[i].vy = ay;
            k2[i].vz = az;
            k2[i].mass = tempBodies[i].mass;
        }

        // --------------------------------------------
        // k3: Evaluate derivatives at t + dt/2 using state = localBodies + (dt/2)*k2
        for (size_t i = 0; i < localBodies.size(); ++i) {
            tempBodies[i].x  = localBodies[i].x  + (timeStep / 2.0) * k2[i].x;
            tempBodies[i].y  = localBodies[i].y  + (timeStep / 2.0) * k2[i].y;
            tempBodies[i].z  = localBodies[i].z  + (timeStep / 2.0) * k2[i].z;
            tempBodies[i].vx = localBodies[i].vx + (timeStep / 2.0) * k2[i].vx;
            tempBodies[i].vy = localBodies[i].vy + (timeStep / 2.0) * k2[i].vy;
            tempBodies[i].vz = localBodies[i].vz + (timeStep / 2.0) * k2[i].vz;
            tempBodies[i].mass = localBodies[i].mass;
        }
        for (size_t i = 0; i < localBodies.size(); ++i) {
            k3[i].x = tempBodies[i].vx;
            k3[i].y = tempBodies[i].vy;
            k3[i].z = tempBodies[i].vz;
            
            double ax = 0.0, ay = 0.0, az = 0.0;
            for (size_t j = 0; j < localBodies.size(); ++j) {
                if (i == j) continue;
                double dx = tempBodies[j].x - tempBodies[i].x;
                double dy = tempBodies[j].y - tempBodies[i].y;
                double dz = tempBodies[j].z - tempBodies[i].z;
                double distSq = dx*dx + dy*dy + dz*dz + softening;
                double dist = std::sqrt(distSq);
                double force = G * tempBodies[j].mass / distSq;
                ax += force * dx / dist;
                ay += force * dy / dist;
                az += force * dz / dist;
            }
            k3[i].vx = ax;
            k3[i].vy = ay;
            k3[i].vz = az;
            k3[i].mass = tempBodies[i].mass;
        }

        // --------------------------------------------
        // k4: Evaluate derivatives at t + dt using state = localBodies + dt*k3
        for (size_t i = 0; i < localBodies.size(); ++i) {
            tempBodies[i].x  = localBodies[i].x  + timeStep * k3[i].x;
            tempBodies[i].y  = localBodies[i].y  + timeStep * k3[i].y;
            tempBodies[i].z  = localBodies[i].z  + timeStep * k3[i].z;
            tempBodies[i].vx = localBodies[i].vx + timeStep * k3[i].vx;
            tempBodies[i].vy = localBodies[i].vy + timeStep * k3[i].vy;
            tempBodies[i].vz = localBodies[i].vz + timeStep * k3[i].vz;
            tempBodies[i].mass = localBodies[i].mass;
        }
        for (size_t i = 0; i < localBodies.size(); ++i) {
            k4[i].x = tempBodies[i].vx;
            k4[i].y = tempBodies[i].vy;
            k4[i].z = tempBodies[i].vz;
            
            double ax = 0.0, ay = 0.0, az = 0.0;
            for (size_t j = 0; j < localBodies.size(); ++j) {
                if (i == j) continue;
                double dx = tempBodies[j].x - tempBodies[i].x;
                double dy = tempBodies[j].y - tempBodies[i].y;
                double dz = tempBodies[j].z - tempBodies[i].z;
                double distSq = dx*dx + dy*dy + dz*dz + softening;
                double dist = std::sqrt(distSq);
                double force = G * tempBodies[j].mass / distSq;
                ax += force * dx / dist;
                ay += force * dy / dist;
                az += force * dz / dist;
            }
            k4[i].vx = ax;
            k4[i].vy = ay;
            k4[i].vz = az;
            k4[i].mass = tempBodies[i].mass;
        }

        // --------------------------------------------
        // Combine increments to update the state:
        // new_state = old_state + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
        for (size_t i = 0; i < localBodies.size(); ++i) {
            newBodies[i].x  = localBodies[i].x  + (timeStep / 6.0) * (k1[i].x  + 2.0 * k2[i].x  + 2.0 * k3[i].x  + k4[i].x);
            newBodies[i].y  = localBodies[i].y  + (timeStep / 6.0) * (k1[i].y  + 2.0 * k2[i].y  + 2.0 * k3[i].y  + k4[i].y);
            newBodies[i].z  = localBodies[i].z  + (timeStep / 6.0) * (k1[i].z  + 2.0 * k2[i].z  + 2.0 * k3[i].z  + k4[i].z);
            newBodies[i].vx = localBodies[i].vx + (timeStep / 6.0) * (k1[i].vx + 2.0 * k2[i].vx + 2.0 * k3[i].vx + k4[i].vx);
            newBodies[i].vy = localBodies[i].vy + (timeStep / 6.0) * (k1[i].vy + 2.0 * k2[i].vy + 2.0 * k3[i].vy + k4[i].vy);
            newBodies[i].vz = localBodies[i].vz + (timeStep / 6.0) * (k1[i].vz + 2.0 * k2[i].vz + 2.0 * k3[i].vz + k4[i].vz);
            newBodies[i].mass = localBodies[i].mass;

            // Save the simulation result for this body and step.
            results.push_back({
                step,
                static_cast<int>(i),
                newBodies[i].x,
                newBodies[i].y,
                newBodies[i].z,
                newBodies[i].vx,
                newBodies[i].vy,
                newBodies[i].vz
            });
        }

        localBodies = std::move(newBodies);
    }
}
