#include "CSR_matrix.hpp"

int main() {
    
    // std::vector<std::vector<double>> matrix1 = {
    //     {0, 5, 0},
    //     {2, 0, 0},
    //     {1, 0, 3},
    //     {0, 0, 2}
    // };

    // std::vector<std::vector<double>> matrix2 = {
    //     {0, 0, 1, 0},
    //     {2, 0, 0, 0},
    //     {0, 0, 3, 1}
    // };

    // CSR_matrix a(matrix1);
    // CSR_matrix b(matrix2);

    // CSR_matrix CSC = a.matrixMultiply(b);
    // CSC.print();


    CSR_matrix a("MM_matrix_files/jgl009.mtx");
    CSR_matrix b("MM_matrix_files/jgl009.mtx");

    a.print();

    // CSR_matrix out = a.matrixMultiply(b);
    // out.print();

    return 0;
}