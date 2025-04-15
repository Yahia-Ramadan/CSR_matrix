#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include "mmio.h"

int main(int argc, char *argv[]) {
    FILE *f;
    MM_typecode matcode;
    int M, N, nz;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " matrix_file.mtx" << std::endl;
        return 1;
    }

    if ((f = fopen(argv[1], "r")) == NULL) {
        std::cerr << "Cannot open file: " << argv[1] << std::endl;
        return 1;
    }

    if (mm_read_banner(f, &matcode) != 0) {
        std::cerr << "Could not process Matrix Market banner.\n";
        return 1;
    }

    if (!mm_is_matrix(matcode) || !mm_is_coordinate(matcode) || !mm_is_real(matcode)) {
        std::cerr << "Only real coordinate matrices are supported.\n";
        return 1;
    }

    mm_read_mtx_crd_size(f, &M, &N, &nz);

    std::vector<int> I(nz), J(nz);
    std::vector<double> val(nz);

    for (int i = 0; i < nz; i++) {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        I[i]--; J[i]--; // convert from 1-based to 0-based
    }

    fclose(f);

    // Print the matrix
    std::cout << "Matrix size: " << M << "x" << N << ", " << nz << " non-zeros.\n";
    for (int i = 0; i < nz; i++) {
        std::cout << "(" << I[i] << ", " << J[i] << ") = " << val[i] << "\n";
    }

    return 0;
}
