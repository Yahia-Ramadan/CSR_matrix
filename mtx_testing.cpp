#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include "mmio.h"

void coo_to_csr(int M, int N, int nnz,
                std::vector<int> &row_coo,
                std::vector<int> &col_coo,
                std::vector<double> &val_coo,
                std::vector<int> &row_ptr,
                std::vector<int> &col_ind,
                std::vector<double> &csr_val) {
    row_ptr.assign(M + 1, 0);
    col_ind.resize(nnz);
    csr_val.resize(nnz);

    // Count non-zeros per row
    for (int i = 0; i < nnz; i++) {
        row_ptr[row_coo[i] + 1]++;
    }

    // Prefix sum
    for (int i = 0; i < M; i++) {
        row_ptr[i + 1] += row_ptr[i];
    }

    std::vector<int> temp_row_ptr = row_ptr;
    for (int i = 0; i < nnz; i++) {
        int row = row_coo[i];
        int dest = temp_row_ptr[row]++;
        col_ind[dest] = col_coo[i];
        csr_val[dest] = val_coo[i];
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " matrix_file.mtx\n";
        return 1;
    }

    FILE *f;
    MM_typecode matcode;
    int M, N, nnz;

    if ((f = fopen(argv[1], "r")) == NULL) {
        std::cerr << "Cannot open file.\n";
        return 1;
    }

    if (mm_read_banner(f, &matcode) != 0) {
        std::cerr << "Could not process Matrix Market banner.\n";
        return 1;
    }

    if (!mm_is_matrix(matcode) || !mm_is_coordinate(matcode) || !mm_is_real(matcode)) {
        std::cerr << "Unsupported Matrix Market format.\n";
        return 1;
    }

    mm_read_mtx_crd_size(f, &M, &N, &nnz);

    std::vector<int> I(nnz), J(nnz);
    std::vector<double> val(nnz);

    for (int i = 0; i < nnz; i++) {
        fscanf(f, "%d %d %lf", &I[i], &J[i], &val[i]);
        I[i]--; // Convert to 0-based
        J[i]--;
    }

    fclose(f);

    std::vector<int> row_ptr, col_ind;
    std::vector<double> csr_val;

    coo_to_csr(M, N, nnz, I, J, val, row_ptr, col_ind, csr_val);

    // Print CSR matrix
    std::cout << "CSR Matrix:\n";
    std::cout << "row_ptr: ";
    for (auto r : row_ptr) std::cout << r << " ";
    std::cout << "\ncol_ind: ";
    for (auto c : col_ind) std::cout << c << " ";
    std::cout << "\nval:     ";
    for (auto v : csr_val) std::cout << v << " ";
    std::cout << std::endl;

    return 0;
}
