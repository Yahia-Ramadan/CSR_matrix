#include "CSR_matrix.h"

int main() {
    std::vector<std::vector<double>> matrix1 = {
        {0, 5, 0},
        {2, 0, 0},
        {1, 0, 3}
    };

    std::vector<std::vector<double>> matrix2 = {
        {0, 0, 1},
        {2, 0, 0},
        {0, 0, 3}
    };

    CSR_matrix a(matrix1);
    CSR_matrix b(matrix2);

    CSR_matrix CSC = a.matrixMultiply(b);
    CSC.print();

    return 0;
}