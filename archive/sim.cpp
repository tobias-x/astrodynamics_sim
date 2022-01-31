#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

// cmake -S . -B build -G Ninja -DUSE_CUDA=ON

#define G 6.67430e-11   // Gravitational constant
#define C 299792458.0   // Speed of light as double to avoid overflow

// Structure to hold body data
struct Body {
    double x, y, z;
    double vx, vy, vz;
    double mass;
    double diameter;
};


// Structure to hold simulation results for each timestep
struct SimulationResult {
    int timestep;
    int bodyIndex;
    double x, y, z, vx, vy, vz; // Newtonian fields
};

const double SOLAR_MASS_IN_KG = 1.9891e30; // Conversion factor from solar masses to kilograms

std::vector<Body> parseConfig(const std::string& filename, int& referenceBody) {
    std::vector<Body> bodies;
    std::ifstream infile(filename);
    std::string line;

    if (!infile.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return bodies;
    }

    std::cout << "Starting to parse the configuration file: " << filename << std::endl;

    // Read the first line for the reference body index
    if (std::getline(infile, line)) {
        std::cout << "First line read for reference body index: " << line << std::endl;
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        // Parse line in the format "label,number"
        std::size_t commaPos = line.find(',');
        if (commaPos != std::string::npos) {
            std::string label = line.substr(0, commaPos);
            std::string value = line.substr(commaPos + 1);
            try {
                referenceBody = std::stoi(value);
                std::cout << "Reference body index parsed: " << referenceBody << std::endl;
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error: Invalid number format in the reference body index: " << value << std::endl;
                return bodies;
            }
        } else {
            std::cerr << "Error: Incorrect format for the reference body index line: " << line << std::endl;
            return bodies;
        }
    } else {
        std::cerr << "Failed to read the first line (reference body index)." << std::endl;
        return bodies;
    }

    // Skip the header line
    if (std::getline(infile, line)) {
        std::cout << "Header line skipped: " << line << std::endl;
    } else {
        std::cerr << "Failed to skip the header line or file is empty after the reference line." << std::endl;
        return bodies;
    }

    // Read each subsequent line and parse the data
    while (std::getline(infile, line)) {
        std::cout << "Reading line: " << line << std::endl;
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

        std::istringstream iss(line);
        std::string name;
        Body body;
        char delimiter;

        std::getline(iss, name, ',');
        if (!(iss >> body.x >> delimiter
                  >> body.y >> delimiter
                  >> body.z >> delimiter
                  >> body.vx >> delimiter
                  >> body.vy >> delimiter
                  >> body.vz >> delimiter
                  >> body.mass >> delimiter
                  >> body.diameter)) {
            std::cerr << "Failed to parse line: " << line << std::endl;
            continue;
        }

        // Convert units to ensure compatibility with gravitational constant
        body.x *= 1000; // Convert from km to meters
        body.y *= 1000;
        body.z *= 1000;
        body.vx *= 1000; // Convert from km/s to m/s
        body.vy *= 1000;
        body.vz *= 1000;
        
        // Convert mass from solar masses to kilograms
        body.mass *= SOLAR_MASS_IN_KG;
        bodies.push_back(body);
        
        std::cout << "Parsed body: " << name << " - x (m): " << body.x << ", y (m): " << body.y 
                  << ", z (m): " << body.z << ", vx (m/s): " << body.vx << ", vy (m/s): " << body.vy 
                  << ", vz (m/s): " << body.vz << ", mass (kg): " << body.mass << ", diameter: " << body.diameter
                  << std::endl;
    }

    std::cout << "Finished parsing configuration. Total bodies parsed: " << bodies.size() << std::endl;
    return bodies;
}

#ifdef USE_CUDA
// CUDA kernel for Newtonian simulation
__global__ void newtonianKernel(Body* d_bodies, int numBodies, double timeStep, int totalSteps, SimulationResult* d_results) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numBodies) return;

    extern __shared__ Body sharedBodies[];
    for (int i = threadIdx.x; i < numBodies; i += blockDim.x) {
        sharedBodies[i] = d_bodies[i];
    }
    __syncthreads();

    Body body = sharedBodies[idx];

    for (int step = 0; step < totalSteps; ++step) {
        double ax = 0.0, ay = 0.0, az = 0.0;
        for (int j = 0; j < numBodies; ++j) {
            if (j != idx) {
                double dx = sharedBodies[j].x - body.x;
                double dy = sharedBodies[j].y - body.y;
                double dz = sharedBodies[j].z - body.z;
                double r_squared = dx*dx + dy*dy + dz*dz;
                double r = sqrt(r_squared);
                double mass_j = sharedBodies[j].mass;
                double force = G * mass_j / r_squared;
                ax += force * dx / r;
                ay += force * dy / r;
                az += force * dz / r;
            }
        }

        body.vx += ax * timeStep;
        body.vy += ay * timeStep;
        body.vz += az * timeStep;
        body.x += body.vx * timeStep;
        body.y += body.vy * timeStep;
        body.z += body.vz * timeStep;

        __syncthreads();

        int resultIdx = step * numBodies + idx;
        d_results[resultIdx] = {step, idx, body.x, body.y, body.z, body.vx, body.vy, body.vz};

        sharedBodies[idx] = body;
        __syncthreads();
    }

    d_bodies[idx] = body;
}

#endif // USE_CUDA

void rungeKutta4CPU(std::vector<Body>& bodies, double timeStep, int totalSteps, std::vector<SimulationResult>& results) {
    int numBodies = bodies.size();
    std::vector<std::tuple<double, double, double>> initialPositions;

    // Store the initial positions of the bodies
    for (const auto& body : bodies) {
        initialPositions.emplace_back(body.x, body.y, body.z);
    }

    std::cout << "Starting Runge-Kutta 4 CPU simulation with " << numBodies << " bodies and " << totalSteps << " total steps." << std::endl;

    for (int step = 0; step < totalSteps; ++step) {
        std::vector<Body> k1(numBodies), k2(numBodies), k3(numBodies), k4(numBodies);
        
        // Calculate k1
        for (int i = 0; i < numBodies; ++i) {
            double ax1 = 0.0, ay1 = 0.0, az1 = 0.0;

            for (int j = 0; j < numBodies; ++j) {
                if (i != j) {
                    double dx = bodies[j].x - bodies[i].x;
                    double dy = bodies[j].y - bodies[i].y;
                    double dz = bodies[j].z - bodies[i].z;
                    double r_squared = dx*dx + dy*dy + dz*dz;
                    double r = sqrt(r_squared);
                    double mass_j = bodies[j].mass;
                    double force = G * mass_j / r_squared;
                    ax1 += force * dx / r;
                    ay1 += force * dy / r;
                    az1 += force * dz / r;
                }
            }

            k1[i] = {bodies[i].vx, bodies[i].vy, bodies[i].vz, ax1, ay1, az1};
        }

        // Calculate k2
        for (int i = 0; i < numBodies; ++i) {
            double ax2 = 0.0, ay2 = 0.0, az2 = 0.0;

            for (int j = 0; j < numBodies; ++j) {
                if (i != j) {
                    double dx = (bodies[j].x + k1[j].x * timeStep / 2) - (bodies[i].x + k1[i].x * timeStep / 2);
                    double dy = (bodies[j].y + k1[j].y * timeStep / 2) - (bodies[i].y + k1[i].y * timeStep / 2);
                    double dz = (bodies[j].z + k1[j].z * timeStep / 2) - (bodies[i].z + k1[i].z * timeStep / 2);
                    double r_squared = dx*dx + dy*dy + dz*dz;
                    double r = sqrt(r_squared);
                    double mass_j = bodies[j].mass;
                    double force = G * mass_j / r_squared;
                    ax2 += force * dx / r;
                    ay2 += force * dy / r;
                    az2 += force * dz / r;
                }
            }

            k2[i] = {bodies[i].vx + k1[i].vx * timeStep / 2, bodies[i].vy + k1[i].vy * timeStep / 2, bodies[i].vz + k1[i].vz * timeStep / 2, ax2, ay2, az2};
        }

        // Calculate k3
        for (int i = 0; i < numBodies; ++i) {
            double ax3 = 0.0, ay3 = 0.0, az3 = 0.0;

            for (int j = 0; j < numBodies; ++j) {
                if (i != j) {
                    double dx = (bodies[j].x + k2[j].x * timeStep / 2) - (bodies[i].x + k2[i].x * timeStep / 2);
                    double dy = (bodies[j].y + k2[j].y * timeStep / 2) - (bodies[i].y + k2[i].y * timeStep / 2);
                    double dz = (bodies[j].z + k2[j].z * timeStep / 2) - (bodies[i].z + k2[i].z * timeStep / 2);
                    double r_squared = dx*dx + dy*dy + dz*dz;
                    double r = sqrt(r_squared);
                    double mass_j = bodies[j].mass;
                    double force = G * mass_j / r_squared;
                    ax3 += force * dx / r;
                    ay3 += force * dy / r;
                    az3 += force * dz / r;
                }
            }

            k3[i] = {bodies[i].vx + k2[i].vx * timeStep / 2, bodies[i].vy + k2[i].vy * timeStep / 2, bodies[i].vz + k2[i].vz * timeStep / 2, ax3, ay3, az3};
        }

        // Calculate k4
        for (int i = 0; i < numBodies; ++i) {
            double ax4 = 0.0, ay4 = 0.0, az4 = 0.0;

            for (int j = 0; j < numBodies; ++j) {
                if (i != j) {
                    double dx = (bodies[j].x + k3[j].x * timeStep) - (bodies[i].x + k3[i].x * timeStep);
                    double dy = (bodies[j].y + k3[j].y * timeStep) - (bodies[i].y + k3[i].y * timeStep);
                    double dz = (bodies[j].z + k3[j].z * timeStep) - (bodies[i].z + k3[i].z * timeStep);
                    double r_squared = dx*dx + dy*dy + dz*dz;
                    double r = sqrt(r_squared);
                    double mass_j = bodies[j].mass;
                    double force = G * mass_j / r_squared;
                    ax4 += force * dx / r;
                    ay4 += force * dy / r;
                    az4 += force * dz / r;
                }
            }

            k4[i] = {bodies[i].vx + k3[i].vx * timeStep, bodies[i].vy + k3[i].vy * timeStep, bodies[i].vz + k3[i].vz * timeStep, ax4, ay4, az4};
        }

        // Update positions and velocities based on RK4
        for (int i = 0; i < numBodies; ++i) {
            bodies[i].x += (k1[i].x + 2 * k2[i].x + 2 * k3[i].x + k4[i].x) * timeStep / 6;
            bodies[i].y += (k1[i].y + 2 * k2[i].y + 2 * k3[i].y + k4[i].y) * timeStep / 6;
            bodies[i].z += (k1[i].z + 2 * k2[i].z + 2 * k3[i].z + k4[i].z) * timeStep / 6;
            bodies[i].vx += (k1[i].vx + 2 * k2[i].vx + 2 * k3[i].vx + k4[i].vx) * timeStep / 6;
            bodies[i].vy += (k1[i].vy + 2 * k2[i].vy + 2 * k3[i].vy + k4[i].vy) * timeStep / 6;
            bodies[i].vz += (k1[i].vz + 2 * k2[i].vz + 2 * k3[i].vz + k4[i].vz) * timeStep / 6;
        }

        for (int i = 0; i < numBodies; ++i) {
            results.push_back({step, i, bodies[i].x, bodies[i].y, bodies[i].z, bodies[i].vx, bodies[i].vy, bodies[i].vz});
        }
    }

    std::cout << "Finished Runge-Kutta 4 CPU simulation." << std::endl;
}

// Function to write output to Parquet file
void writeCSV(const std::vector<SimulationResult>& results, const std::string& filename) {
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << " for writing." << std::endl;
        return;
    }

    // Write the CSV header
    outfile << "timestep,bodyIndex,x,y,z,vx,vy,vz" << std::endl;

    // Write the data for each result
    for (const auto& result : results) {
        outfile << result.timestep << ","
                << result.bodyIndex << ","
                << result.x << ","
                << result.y << ","
                << result.z << ","
                << result.vx << ","
                << result.vy << ","
                << result.vz << std::endl;
    }

    std::cout << "Successfully wrote results to CSV file: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    std::cout << "Starting the simulation..." << std::endl;

    if (argc < 7) {
        std::cerr << "Usage: ./simulation <config_file> <time_steps> <step_size> <output_file> <use_gpu>\n";
        return -1;
    }

    std::string configFile = argv[1];
    int totalSteps = std::stoi(argv[2]);
    double timeStep = std::stod(argv[3]);
    std::string outputFile = argv[4];
    bool isNewtonian = std::stoi(argv[5]);
    bool useGPU = std::stoi(argv[6]);

    std::cout << "Configuration file: " << configFile << std::endl;
    std::cout << "Total steps: " << totalSteps << ", Time step size: " << timeStep << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;
    std::cout << ", Use GPU: " << useGPU << std::endl;

    int referenceBody = 0;
    std::vector<Body> bodies = parseConfig(configFile, referenceBody);
    int numBodies = bodies.size();
    std::vector<SimulationResult> results;

    if (useGPU) {
#ifdef USE_CUDA
        std::cout << "Running Newtonian simulation on GPU..." << std::endl;
        Body* d_bodies;
        SimulationResult* d_results;
        cudaMalloc((void**)&d_bodies, numBodies * sizeof(Body));
        cudaMalloc((void**)&d_results, numBodies * totalSteps * sizeof(SimulationResult));
        cudaMemcpy(d_bodies, bodies.data(), numBodies * sizeof(Body), cudaMemcpyHostToDevice);

        int threadsPerBlock = 256;
        int blocksPerGrid = (numBodies + threadsPerBlock - 1) / threadsPerBlock;
        size_t sharedMemSize = numBodies * sizeof(Body);
        newtonianKernel<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(d_bodies, numBodies, timeStep, totalSteps, d_results);
        cudaDeviceSynchronize();

        results.resize(numBodies * totalSteps);
        cudaMemcpy(results.data(), d_results, numBodies * totalSteps * sizeof(SimulationResult), cudaMemcpyDeviceToHost);

        std::cout << "GPU Newtonian simulation complete." << std::endl;

        cudaFree(d_bodies);
        cudaFree(d_results);
#else
        std::cerr << "CUDA is not supported on this build. Running CPU-only version.\n";
        rungeKutta4CPU(bodies, timeStep, totalSteps, results);
#endif
    } else {
        rungeKutta4CPU(bodies, timeStep, totalSteps, results);
    }

    // Call the new CSV writer function
    writeCSV(results, outputFile);

    std::cout << "Simulation completed successfully." << std::endl;

    return 0;
}
