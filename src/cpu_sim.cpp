#include "cpu_sim.h"
#include <cmath>
#include <iostream>

static constexpr double G_CONST = 6.6743e-11;
static constexpr double SOFTENING = 1e-5;

void CpuSim::step(std::vector<Body>& bodies, double dt) const {
    int n = bodies.size();
    if (n == 0) return;

    struct Deriv { double vx, vy, vz, ax, ay, az; };

    std::vector<Deriv> k1(n), k2(n), k3(n), k4(n);
    std::vector<Body> temp(n), next(n);

    for (int i = 0; i < n; ++i) {
        k1[i].vx = bodies[i].vx;
        k1[i].vy = bodies[i].vy;
        k1[i].vz = bodies[i].vz;
        double ax = 0, ay = 0, az = 0;
        for (int j = 0; j < n; ++j) if (i != j) {
            double dx = bodies[j].x - bodies[i].x;
            double dy = bodies[j].y - bodies[i].y;
            double dz = bodies[j].z - bodies[i].z;
            double distSq = dx * dx + dy * dy + dz * dz + SOFTENING;
            double invDist = 1.0 / std::sqrt(distSq);
            double f = G_CONST * bodies[j].mass * invDist * invDist;
            ax += f * dx * invDist;
            ay += f * dy * invDist;
            az += f * dz * invDist;
        }
        k1[i].ax = ax; k1[i].ay = ay; k1[i].az = az;
    }

    auto compute = [&](const std::vector<Deriv>& kPrev, std::vector<Deriv>& kOut, double scale) {
        for (int i = 0; i < n; ++i) {
            temp[i].x = bodies[i].x + kPrev[i].vx * dt * scale;
            temp[i].y = bodies[i].y + kPrev[i].vy * dt * scale;
            temp[i].z = bodies[i].z + kPrev[i].vz * dt * scale;
            temp[i].vx = bodies[i].vx + kPrev[i].ax * dt * scale;
            temp[i].vy = bodies[i].vy + kPrev[i].ay * dt * scale;
            temp[i].vz = bodies[i].vz + kPrev[i].az * dt * scale;
            temp[i].mass = bodies[i].mass;
        }
        for (int i = 0; i < n; ++i) {
            kOut[i].vx = temp[i].vx;
            kOut[i].vy = temp[i].vy;
            kOut[i].vz = temp[i].vz;
            double ax = 0, ay = 0, az = 0;
            for (int j = 0; j < n; ++j) if (i != j) {
                double dx = temp[j].x - temp[i].x;
                double dy = temp[j].y - temp[i].y;
                double dz = temp[j].z - temp[i].z;
                double distSq = dx * dx + dy * dy + dz * dz + SOFTENING;
                double invDist = 1.0 / std::sqrt(distSq);
                double f = G_CONST * temp[j].mass * invDist * invDist;
                ax += f * dx * invDist;
                ay += f * dy * invDist;
                az += f * dz * invDist;
            }
            kOut[i].ax = ax; kOut[i].ay = ay; kOut[i].az = az;
        }
    };

    compute(k1, k2, 0.5);
    compute(k2, k3, 0.5);
    compute(k3, k4, 1.0);

    for (int i = 0; i < n; ++i) {
        next[i].name = bodies[i].name;
        next[i].mass = bodies[i].mass;
        next[i].x = bodies[i].x + dt / 6.0 * (k1[i].vx + 2 * k2[i].vx + 2 * k3[i].vx + k4[i].vx);
        next[i].y = bodies[i].y + dt / 6.0 * (k1[i].vy + 2 * k2[i].vy + 2 * k3[i].vy + k4[i].vy);
        next[i].z = bodies[i].z + dt / 6.0 * (k1[i].vz + 2 * k2[i].vz + 2 * k3[i].vz + k4[i].vz);
        next[i].vx = bodies[i].vx + dt / 6.0 * (k1[i].ax + 2 * k2[i].ax + 2 * k3[i].ax + k4[i].ax);
        next[i].vy = bodies[i].vy + dt / 6.0 * (k1[i].ay + 2 * k2[i].ay + 2 * k3[i].ay + k4[i].ay);
        next[i].vz = bodies[i].vz + dt / 6.0 * (k1[i].az + 2 * k2[i].az + 2 * k3[i].az + k4[i].az);
    }

    bodies = std::move(next);
}

void CpuSim::runSim(const std::vector<Body>& bodies,
                    double timeStep,
                    int totalSteps,
                    std::vector<SimulationResult>& results) const
{
    std::vector<Body> state = bodies;
    for (int stepIndex = 0; stepIndex < totalSteps; ++stepIndex) {
        step(state, timeStep);
        for (int i = 0; i < (int)state.size(); ++i) {
            results.push_back({
                stepIndex, i,
                state[i].x, state[i].y, state[i].z,
                state[i].vx, state[i].vy, state[i].vz
            });
        }
    }
}

void CpuSim::calcGrid(std::vector<GridBody>& gridBodies, const std::vector<Body>& bodies) const {
    if (bodies.empty() || gridBodies.empty()) return;

    const double soften = 1e-5; // small number to avoid divide-by-zero
    const double baseMass = 1e28; // grid point mass

    for (auto& gb : gridBodies) {
        double weightedX = baseMass * gb.initial_x;
        double weightedY = baseMass * gb.initial_y;
        double weightedZ = baseMass * gb.initial_z;
        double totalWeight = baseMass;

        for (const auto& body : bodies) {
            double dx = body.x - gb.initial_x;
            double dy = body.y - gb.initial_y;
            double dz = body.z - gb.initial_z;

            double distance = std::abs(dx) + std::abs(dy) + std::abs(dz) + soften; // <<< NO square, NO sqrt, purely linear

            double weight = body.mass / distance; // <<< LINEAR

            weightedX += weight * body.x;
            weightedY += weight * body.y;
            weightedZ += weight * body.z;
            totalWeight += weight;
        }

        gb.x = weightedX / totalWeight;
        gb.y = weightedY / totalWeight;
        gb.z = weightedZ / totalWeight;
    }
}
