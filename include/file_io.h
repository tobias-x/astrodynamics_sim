#pragma once
#include <vector>
#include <string>
#include "types.h"
std::vector<Body> parseConfig(const std::string& filename);
void writeCSV(const std::vector<SimulationResult>& results, const std::string& filename);
