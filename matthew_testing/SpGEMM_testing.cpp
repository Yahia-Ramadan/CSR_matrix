#include <iostream>
#include <eigen3/Eigen/Sparse>
#include <cstdlib>
#include <ctime>
#include <cmath>
// #include "mergeKSorted.cpp"
#include "../CSR_matrix.hpp"
#include <vector>

using namespace Eigen;
using namespace std;

// Generates Sparse Matrix
SparseMatrix<double> createRandomSparseMatrix(int rows, int cols, double density) {
    SparseMatrix<double> mat(rows, cols);
    std::vector<Triplet<double>> triplets;
    int maxNonZeros = floor(rows * cols * density);

    for (int i = 0; i < maxNonZeros; ++i) {
        int row = rand() % rows;
        int col = rand() % cols;
        double value = static_cast<double>(rand()) / RAND_MAX; // value between [0,1]
        triplets.emplace_back(row, col, value);
    }

    mat.setFromTriplets(triplets.begin(), triplets.end());
    mat.makeCompressed(); // CSR format
    return mat;
}

CSR_matrix SpGEMM(CSR_matrix matA, CSR_matrix matB) {
    vector<double> result_values;
    vector<int> result_col;
    vector<int> result_row = {0};

    vector<int> rowA = matA.getRow();
    vector<int> rowB = matB.getRow();

    vector<int> colA = matA.getCol();
    vector<int> colB = matB.getCol();

    vector<double> valuesA = matA.getValues();
    vector<double> valuesB = matB.getValues();

    // for each row in A
    for(int i=0; i<rowA.size()-1; i++) {
        unordered_map<int, double> mergedRow; 
        // cout << "rowAStart: " << rowA[i] << " rowAEnd: " << rowA[i+1] << endl;
        // for each nonzero index j in the current row in A
        for (int j = rowA[i]; j < rowA[i + 1]; j++) {
            double currentValueA = valuesA[j];
            // find row in B corresponding to column of value in A
            int rowB_start = rowB[colA[j]];
            int rowB_end = rowB[colA[j]+1];

            // cout << "rowBStart: " << rowB_start << " rowBEnd: " << rowB_end << endl;
            // for each nonzero index k for corresponding row in B
            for (int k = rowB_start; k<rowB_end; k++) {
                int key = colB[k];
                double value = currentValueA * valuesB[k];
                // cout << "currentValueA: " << currentValueA << " valuesB[k]: " << valuesB[k] << " ";
                // cout << "key, value: " << key << " " << value << endl;
                mergedRow[key] += value;
            }

        }
        
        // mergeRow contains properly merged values, sorted makes the keys sorted
        map<int, double> sorted(mergedRow.begin(), mergedRow.end());

        for(auto item: sorted) {
            int key = item.first;
            double value = item.second;
            if(value != 0.0) { // even though mergedRow contains nonzero elements, adding multiple elements in same column may result in 0
                
                // is there any way to do it without push_back?
                result_values.push_back(value);
                result_col.push_back(key);

            }
        }
        result_row.push_back(result_values.size());
    }

    CSR_matrix matC(result_values, result_col, result_row);
    return matC;
}

CSR_matrix makeCSR(SparseMatrix<double, RowMajor> matrix) {
    vector<double> values(matrix.nonZeros());
    vector<int> col(matrix.nonZeros());
    vector<int> row(matrix.rows() + 1);

    for (int i = 0; i < matrix.nonZeros(); i++) {
        values[i] = matrix.valuePtr()[i];      
        col[i] = matrix.innerIndexPtr()[i];
    }

    for (int i = 0; i <= matrix.rows(); i++) {
        row[i] = matrix.outerIndexPtr()[i];
    }
    CSR_matrix res(values, col, row);
    return res;
}

// // HashMap implementation
// // Note that pair<int, double> is column_index, value
// vector<pair<int,double>> mergeKListsHashMap(vector<vector<pair<int, double>>> lists) {
//     if(lists.empty()) {
//         return vector<pair<int,double>>{};
//     }

//     unordered_map<int, double> result_map; // key = column index, key value = value

//     for(int i=0; i<lists.size(); i++) {
//         for(int j=0; j<lists[i].size(); j++) {
//             // lists[i][col] gives pair<int, double>
//             int column_index = lists[i][j].first;
//             int value = lists[i][j].second;
//             if (result_map.find(column_index) == result_map.end()) { // does not exist in map
//                 result_map[column_index] = value;
//             } else {
//                 result_map[column_index] += value;
//             }
//         }
//     }
// }


int main() {
    srand(42); // seed

    int rows = 100, cols = 100;
    double density = 0.5;
    SparseMatrix<double, RowMajor> sm1 = createRandomSparseMatrix(rows, cols, density);
    SparseMatrix<double, RowMajor> sm2 = createRandomSparseMatrix(rows, cols, density);

    SparseMatrix<double, RowMajor> Eigen_Result = sm1 * sm2;
    CSR_matrix A = makeCSR(sm1);
    CSR_matrix B = makeCSR(sm2);

    // cout << "Eigen Results: " << endl;
    // for (int k = 0; k < Eigen_Result.outerSize(); ++k)
    //     for (SparseMatrix<double, RowMajor>::InnerIterator it(Eigen_Result, k); it; ++it)
    //         cout << "(" << it.row() << "," << it.col() << ") = " << it.value() << endl;
    // cout << endl << endl;

    CSR_matrix SpGEMM_Result = SpGEMM(A, B);

    // vector<double> values = SpGEMM_Result.getValues();
    // vector<int> col = SpGEMM_Result.getCol();
    // vector<int> row = SpGEMM_Result.getRow();
    
    // Get non-zero count and dimensions
    int nnz = Eigen_Result.nonZeros();
    int numRows = Eigen_Result.rows();

    // Extract pointers
    const double* valPtr = Eigen_Result.valuePtr();       // values array
    const int*    colPtr = Eigen_Result.innerIndexPtr();  // column indices (for RowMajor: these are columns)
    const int*    rowPtr = Eigen_Result.outerIndexPtr();  // row pointers (CSR-style)

    // Convert to vectors
    std::vector<double> values(valPtr, valPtr + nnz);
    std::vector<int>    columns(colPtr, colPtr + nnz);
    std::vector<int>    rowss(rowPtr, rowPtr + numRows + 1);  // +1 because CSR row array is size (n_rows + 1)

    if(values == SpGEMM_Result.getValues() && columns == SpGEMM_Result.getCol() && rowss == SpGEMM_Result.getRow()) {
        cout << "Results match!" << endl;
    }
    else {
        cout << "Results do not match..." << endl;
    }
    // cout << "values.size(): " << values.size() << endl;
    // cout << "col.size(): " << col.size() << endl;
    // cout << "row.size(): " << row.size() << endl;

    // for(int i=0; i<values.size(); i++) {
    //     cout << values[i] << " ";
    // }
    // cout << endl;

    // for(int i=0; i<col.size(); i++) {
    //     cout << col[i] << " ";
    // }
    // cout << endl;

    // for(int i=0; i<row.size(); i++) {
    //     cout << row[i] << " ";
    // }
    // cout << endl;

    // for(int i=0; i<row.size()-1; i++) {
    //     for(int k=row[i]; k<row[i+1]; k++) {
    //         cout << "(" << i << "," << col[k] << ") = " << values[k] << endl;
    //     }
    // }

    return 0;
}
