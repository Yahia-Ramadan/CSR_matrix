#include <iostream>
#include <vector>

using namespace std;

void CSR(const vector<vector<int>>& matrix, 
    vector<int>& values, vector<int>& col_index, vector<int>& row_ptr) {

    int rows = matrix.size();
    int cols = matrix[0].size();
    
    values = {};
    col_index = {};
    row_ptr = {0};

    // loop through matrix by row (left to right)
    for(int i=0; i<rows; i++) {
        int num_elements = 0;
        for(int j=0; j<cols; j++) {
            if(matrix[i][j] != 0) {
                values.push_back(matrix[i][j]);
                col_index.push_back(j);
                num_elements++;
            }
        }
        row_ptr.push_back(row_ptr.back() + num_elements);
    }
}

// vector multiplication (A*x), returns resulting vector, x=vector
vector<int> VectorMult(vector<int>& x, vector<int>& values,
    vector<int>& col_index, vector<int>& row_ptr) {

    int rows = row_ptr.size() - 1; // total rows = size of row_ptr-1
    vector<int> result;
    
    for(int i=0; i<rows; i++) { // loop through each row
        int sum = 0;
        for(int j=row_ptr[i]; j<row_ptr[i+1]; j++) {    // contains indices in values corresponding to this row
            sum += values[j] * x[col_index[j]];
        }
        result.push_back(sum);
    }
    return result;    
}

vector<vector<int>> MatrixMultDense(vector<vector<int>>& matrix_B, vector<int>& values,
    vector<int>& col_index, vector<int>& row_ptr) {
    // matrix A in CSR form
    int rowsA = row_ptr.size() - 1; // (rowA x colA) * (rowB x colB) = (rowA x colB)
    int colsB = matrix_B[0].size();
    
    // cout << "rowsA: " << rowsA << " colsB: " << colsB << endl;
    vector<vector<int>> result(rowsA, vector<int>(colsB, 0)); // initialize result matrix
        
    for(int i=0; i<rowsA; i++) {    // each row in A
        for(int j=row_ptr[i]; j<row_ptr[i+1]; j++) { // select only non-zero values
            int val = values[j];
            int col = col_index[j];

            for(int k=0; k<colsB; k++) {
                result[i][k] += val * matrix_B[col][k];
            }
        }
    }

    return result;
}

void print(vector<int> vector) {
    cout << "[ ";
    for(const auto& element: vector) {
        cout << element << " ";
    }
    cout << "]" << endl;
}

int main() {


    // 3 arrays: values, col_index, and row_ptr
    vector<int> values;
    vector<int> col_index;
    vector<int> row_ptr;

    // Example matrix for testing
    vector<vector<int>> matrix = {
        {0, 8, 0, 0, 0}, 
        {0, 0, 0, 0, 0}, 
        {5, 0, 0, 0, 9},
        {0, 0, 3, 0, 0}
    };

    CSR(matrix, values, col_index, row_ptr);

    //// TESTING ARRAY IMPLEMENTATION ////
    // cout << "values: ";
    // print(values);

    // cout << "col_index: ";
    // print(col_index);
    
    // cout << "row_ptr: ";
    // print(row_ptr);
    /////////////////////////////////////


    //// TESTING VECTOR MULTIPLICATION (DENSE) ////
    // vector<int> x = {1,2,3,4,5};
    // vector<int> result = VectorMult(x, values, col_index, row_ptr);

    // cout << "result vector: ";
    // print(result);
    ///////////////////////////////////////////////


    //// TESTING MATRIX MULTIPLICATION (DENSE) ////
    vector<vector<int>> matrix_B = {
        {1,2,3},
        {4,5,6},
        {7,8,9},
        {0,0,1},
        {2,1,0},
    };
    vector<vector<int>> result_matrix = MatrixMultDense(matrix_B, values, col_index, row_ptr);

    cout << "result matrix:" << endl;
    for (int i=0; i<result_matrix.size(); i++) {
        for (int j=0; j<result_matrix[i].size(); j++) {
          std::cout << result_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    return 0;
    ///////////////////////////////////////////////



    //// TESTING MATRIX MULTIPLICATION (SPARSE) ////




    ///////////////////////////////////////////////




    return 0;
}

