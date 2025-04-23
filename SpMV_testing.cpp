#include "CSR_matrix.hpp"
#include "mkl_spblas.h"

#include <chrono>
#include <iostream>
#include <random>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>

namespace fs = std::filesystem;

// Helper to get all .mtx files from a directory
std::vector<std::string> get_matrix_files(const std::string& folder_path) {
    std::vector<std::string> files;
    for (const auto& entry : fs::directory_iterator(folder_path)) {
        if (entry.path().extension() == ".mtx") {
            files.push_back(entry.path().string());
        }
    }
    return files;
}

// Random float generator
float generate_random_float(float min, float max, std::mt19937& gen) {
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

int main() {
    std::string folder_path = "MM_matrix_files/";
    std::ofstream log("benchmark_results.csv");
    log << "Matrix,Rows,Cols,NNZ,CustomAvgTimeMS,MKLAvgTimeMS\n";

    std::vector<std::string> matrix_files = get_matrix_files(folder_path);

    for (const auto& file : matrix_files) {
        std::cout << "Processing: " << file << std::endl;

        // Load matrix once
        CSR_matrix a(file);
        int rows = a.getNCols(); 
        int cols = a.getNRows();
        int nnz = a.getNNZ();

        std::vector<double> vec(cols, 0.0);
        std::mt19937 gen(std::random_device{}());
        for (auto& v : vec) v = generate_random_float(-100.0f, 100.0f, gen);

        double total_time_mine = 0.0;
        double total_time_mkl  = 0.0;

        // MKL setup
        sparse_matrix_t A;
        std::vector<double> values = a.getValues();
        std::vector<int> col_index = a.getCol();
        std::vector<int> row_pointer = a.getRow();

        sparse_status_t status = mkl_sparse_d_create_csr(
            &A, SPARSE_INDEX_BASE_ZERO,
            rows, cols,
            row_pointer.data(),
            row_pointer.data() + 1,
            col_index.data(),
            values.data()
        );

        if (status != SPARSE_STATUS_SUCCESS) {
            std::cerr << "MKL failed to create CSR matrix for: " << file << "\n";
            continue;
        }

        struct matrix_descr descr;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        mkl_sparse_optimize(A);

        // Run both methods 250 times
        for (int i = 0; i < 250; ++i) {
            std::vector<double> res(rows, 0.0);

            auto start1 = std::chrono::high_resolution_clock::now();
            std::vector<double> out = a.SpMV(res, vec, 2.5, -0.5, 1);
            auto end1 = std::chrono::high_resolution_clock::now();
            total_time_mine += std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1).count() / 1e6;

            std::vector<double> y(rows, 0.0);
            auto start2 = std::chrono::high_resolution_clock::now();
            status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 2.5, A, descr, vec.data(), -0.5, y.data());
            auto end2 = std::chrono::high_resolution_clock::now();
            total_time_mkl += std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2).count() / 1e6;

            if (status != SPARSE_STATUS_SUCCESS) {
                std::cerr << "MKL multiplication failed on iteration " << i << " for file: " << file << "\n";
                break;
            }
        }

        double avg_custom = total_time_mine / 250.0;
        double avg_mkl = total_time_mkl / 250.0;

        std::cout << "Avg custom SpMV: " << avg_custom << " ms\n";
        std::cout << "Avg MKL    SpMV: " << avg_mkl    << " ms\n\n";

        log << file << "," << rows << "," << cols << "," << nnz << ","
            << avg_custom << "," << avg_mkl << "\n";

        mkl_sparse_destroy(A);
    }

    log.close();
    std::cout << "Benchmarking complete. Results saved to 'benchmark_results.csv'\n";
    return 0;
}