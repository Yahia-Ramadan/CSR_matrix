#include "CSR_matrix.hpp"
#include <mkl.h>

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

bool compare_vectors(const std::vector<double>& a, const std::vector<double>& b, double eps = 1e-9) {
    if (a.size() != b.size()) {
    std::cerr << "Vectors have different sizes!\n";
    return false;
    }

    for (size_t i = 0; i < a.size(); ++i) {
        double diff = std::abs(a[i] - b[i]);
        if (diff > eps) {
            std::cerr << "Mismatch at index " << i << ": a = " << a[i] << ", b = " << b[i] << ", diff = " << diff << "\n";
        return false;
        }
    }

    return true;
}

int main(){
    mkl_set_num_threads(1); // VM default only has 2 cores

    /************************************************
     * 
     * 
     *      First, execute this library's version
     *      of simple SpMV
     * 
     * 
    * *************************************************/
    int m = 17758;  // rows
    int n = 17758;  // columns

    //Load matrix
    auto start = std::chrono::high_resolution_clock::now();
    CSR_matrix a("MM_matrix_files/memplus.mtx");
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "\nThis library's time to load matrix: " << duration.count() << " ms" << std::endl;

    //Generate random dense vector with size 17758
    float min_value = -100.0;
    float max_value = 100.0;

    std::vector<double> vec(m, 0.0);

    for (size_t i = 0; i < 17758; i++)
    {
        vec.at(i) = generate_random_float(min_value, max_value);
    }
    
    //Multiply
 
    auto start1 = std::chrono::high_resolution_clock::now();
    std::vector<double> res(m, 0.0);
    // a.SpMV(res, vec, 2.5, -0.5);
    std::vector<double> result = a.vectorMultiply(vec); // using traditional vectorMultiply
    
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

    // Optional optimization
    mkl_sparse_optimize(A); 

    auto e_mkl_csr = std::chrono::high_resolution_clock::now();

    auto mkl_csr = std::chrono::duration_cast<std::chrono::milliseconds>(e_mkl_csr - s_mkl_csr);
    std::cout << "Intel MKL's time to load matrix: " << mkl_csr.count() << " ms" << std::endl;

    // Set matrix description
    struct matrix_descr descr;
    descr.type = SPARSE_MATRIX_TYPE_GENERAL;
    descr.diag = SPARSE_DIAG_NON_UNIT;  // diagonal's not all default 1.0

    std::vector<double> y(m, 0.0);  // size = m

    // Perform y = A * x

    auto s_mkl_mult = std::chrono::high_resolution_clock::now();

    // status = mkl_sparse_d_mv(
    //     SPARSE_OPERATION_NON_TRANSPOSE,
    //     1.0, 
    //     A,
    //     descr,
    //     vec.data(),
    //     0.0, 
    //     y.data()
    // );

    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, A, descr, vec.data(), 0.0, y.data());

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

    if (compare_vectors(result, y)) {
        std::cout << "Implementation matches MKL" << std::endl;
    } else {
        std::cout << "Mismatch detected." << std::endl;
    }

    return 0;

}