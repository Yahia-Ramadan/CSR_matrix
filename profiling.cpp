#include "CSR_matrix.hpp"

#include <chrono>
#include <iostream>
#include <random>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <valgrind/callgrind.h>


float generate_random_float(float min, float max, std::mt19937& gen) {
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

int main() {
    
    CSR_matrix a("MM_matrix_files/rajat31.mtx");
    std::vector<double> vec(a.getNRows(), 0.0);
    std::vector<double> res(a.getNCols(), 0.0);
    
    a.SpMV(res, vec, 2.5, -0.5, 1);

    
    CALLGRIND_START_INSTRUMENTATION;
    
    a.SpMV(res, vec, 2.5, -0.5, 8);
    
    CALLGRIND_STOP_INSTRUMENTATION;

    return 0;
}