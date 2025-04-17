#include "CSR_matrix.hpp"
#include "mkl_spblas.h"
#include <mkl.h>
#include <mkl_spblas.h>

#include <chrono>
#include <iostream>
#include <random>
#include <limits>

float generate_random_float(float min, float max) {
    // Seed the random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Define the distribution
    std::uniform_real_distribution<> dis(min, max);

    // Generate and return the random float
    return dis(gen);
}

int main(){

    /************************************************
     * 
     * 
     *      First, execute this library's version
     *      of simple SpMV
     * 
     * 
    * *************************************************/


    //Load matrix
    auto start = std::chrono::high_resolution_clock::now();
    CSR_matrix a("MM_matrix_files/memplus.mtx");
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "\nThis library's time to load matrix: " << duration.count() << " ms" << std::endl;

    //Generate random dense vector with size 17758
    float min_value = -100.0;
    float max_value = 100.0;

    std::vector<double> vec(17758, 0.0);
    std::vector<double> res(17758, 0.0);
    for (size_t i = 0; i < 17758; i++)
    {
        vec.at(i) = generate_random_float(min_value, max_value);
    }
    
    //Multiply
 
    auto start1 = std::chrono::high_resolution_clock::now();
    std::vector<double> out;
    out = a.SpMV(res, vec, 2.5, -0.5);

    auto end1 = std::chrono::high_resolution_clock::now();
    auto duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1);
    double elapsed_ms = duration_ns.count() / 1'000'000.0;
    std::cout << "This library's time to perform SpMV: " << elapsed_ms << " ms\n" << std::endl;



    /************************************************
     * 
     * 
     *      Now, execute this Intel OpenMKL's verision
     *      of SpMV
     * 
     * 
    * *************************************************/
    
    
    std::vector<double> values = a.getValues();
    std::vector<int> col_index   = a.getCol();
    std::vector<int> row_pointer = a.getRow();

    int m = 17758;  // rows
    int n = 17758;  // columns
 
    // Create sparse matrix handle
    sparse_matrix_t A;
    sparse_status_t status;


    auto s_mkl_csr = std::chrono::high_resolution_clock::now();

    status = mkl_sparse_d_create_csr(
        &A, SPARSE_INDEX_BASE_ZERO,
        m, n,
        row_pointer.data(),
        row_pointer.data() + 1,
        col_index.data(),
        values.data()
    );

    if (status != SPARSE_STATUS_SUCCESS) {
        std::cerr << "Failed to create CSR matrix.\n";
        return 1;
    }

    auto e_mkl_csr = std::chrono::high_resolution_clock::now();

    auto mkl_csr = std::chrono::duration_cast<std::chrono::milliseconds>(e_mkl_csr - s_mkl_csr);
    std::cout << "Intel MKL's time to load matrix: " << mkl_csr.count() << " ms" << std::endl;

    // Set matrix description
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;

    // Optional optimization
    mkl_sparse_optimize(A);

    // Input vector
    std::vector<double> x = vec; // size = n
    std::vector<double> y(m, 0.0);  // size = m

    // Perform y = A * x

    auto s_mkl_mult = std::chrono::high_resolution_clock::now();

    status = mkl_sparse_d_mv(
        SPARSE_OPERATION_NON_TRANSPOSE,
        1.0, A,
        descr,
        x.data(),
        0.0, y.data()
    );

    auto e_mkl_mult = std::chrono::high_resolution_clock::now();
    if (status != SPARSE_STATUS_SUCCESS) {
        std::cerr << "Matrix-vector multiplication failed.\n";
        return 1;
    }

    auto duration_ns_mkl = std::chrono::duration_cast<std::chrono::nanoseconds>(e_mkl_mult - s_mkl_mult);
    double elapsed_ms_mkl = duration_ns_mkl.count() / 1'000'000.0;
    std::cout << "Intel MKL's time to perform SpMV: " << elapsed_ms_mkl << " ms\n" << std::endl;


    // Clean up
    mkl_sparse_destroy(A);

    return 0;

}