#ifdef USE_CUDA
#include "cuda_sim.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include <iostream>

static constexpr double G_CONST = 6.6743e-11;
static constexpr double SOFTENING = 1e-5;

__global__ void newtonianKernel(Body* bodies, int n, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    extern __shared__ Body shared[];
    for (int j=threadIdx.x; j<n; j+=blockDim.x)
        shared[j] = bodies[j];
    __syncthreads();

    Body me = shared[i];
    bodies[i] = me;
}

void CudaSim::step(std::vector<Body>& bodies, double dt) const {
    int n = bodies.size();
    Body* d_bodies;
    cudaMalloc(&d_bodies, n*sizeof(Body));
    cudaMemcpy(d_bodies, bodies.data(), n*sizeof(Body), cudaMemcpyHostToDevice);

    int block = 128, grid = (n+block-1)/block;
    newtonianKernel<<<grid,block,n*sizeof(Body)>>>(d_bodies, n, dt);
    cudaDeviceSynchronize();

    cudaMemcpy(bodies.data(), d_bodies, n*sizeof(Body), cudaMemcpyDeviceToHost);
    cudaFree(d_bodies);
}

void CudaSim::runSim(const std::vector<Body>& bodies,
                     double timeStep,
                     int totalSteps,
                     std::vector<SimulationResult>& results) const
{
    std::vector<Body> state = bodies;
    for (int step=0; step<totalSteps; ++step) {
        step(state, timeStep);
        for (int i=0; i<(int)state.size(); ++i) {
            results.push_back({
                step, i,
                state[i].x, state[i].y, state[i].z,
                state[i].vx, state[i].vy, state[i].vz
            });
        }
    }
}
#endif
