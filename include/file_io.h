#ifndef FILE_IO_H
#define FILE_IO_H

#include <string>
#include <vector>
#include "simulation_common.h"

std::vector<Body> parseConfig(const std::string& filename, int& referenceBody);
void writeCSV(const std::vector<SimulationResult>& results, const std::string& filename);

#endif // FILE_IO_H
