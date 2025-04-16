#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdexcept>

class CSR_matrix {
private:
    std::vector<double> values;
    std::vector<int> col;
    std::vector<int> row;

    //TODO: update constructors to support the new values;
    // int n_cols;
    // int n_rows;
    // int n_nonzero;

public:

    CSR_matrix();
    CSR_matrix(std::vector<std::vector<double>>& matrix);
    CSR_matrix(CSR_matrix& cpy);
    CSR_matrix(std::vector<double> v, std::vector<int> c, std::vector<int> r);
    CSR_matrix(const std::string& filename);

    ~CSR_matrix();

    std::vector<std::vector<double>> toMatrix();
    std::vector<double> vectorMultiply(std::vector<double>& vec);
    CSR_matrix matrixMultiply(CSR_matrix& other);
    void print();
    
    double dotHelper( CSR_matrix& other, int idx1, int idx);
    CSR_matrix toCSC();

    // Getters
    const std::vector<double>& getValues() const { return values; }
    const std::vector<int>& getCol() const { return col; }
    const std::vector<int>& getRow() const { return row; }
};

CSR_matrix::CSR_matrix() = default;


CSR_matrix::CSR_matrix(std::vector<std::vector<double>>& matrix) {

    int len = matrix.size();
    int len2 = matrix[0].size();

    row.resize(len+1, 0);

    for (int i = 0; i < len; i++){
        row[i] = values.size();
        for (int j = 0; j < len2; j++){
            if (matrix[i][j] != 0) {
                values.push_back(matrix[i][j]);
                col.push_back(j);
            }
        }   
    }

    row[len] = values.size();
}

CSR_matrix::CSR_matrix(CSR_matrix& cpy): values(cpy.values), col(cpy.col), row(cpy.row) {}

CSR_matrix::CSR_matrix(std::vector<double> v, std::vector<int> c, std::vector<int> r): values(v), col(c), row(r) {}


CSR_matrix::CSR_matrix(const std::string& filename){

    //load file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << "\n";
        return;
    }

    //first parse header line

    //can add support for pattern and complex later, as well as symmetries other than general
    //must be matrix, coordinate, real, general

    std::string line;
    std::getline(file, line);
    std::stringstream header(line);
    try {
        std::string curr;
        
        //skip %%MatrixMarket portion
        header >> curr;

        //parse rest of line
        header >> curr;
        if (curr != "matrix") {
            std::cerr << "Error: First line of matrix file must be: /%/%MatrixMarket matrix coordinate real general\n";
            return;
        }
        header >> curr;
        if (curr != "coordinate") {
            std::cerr << "Error: First line of matrix file must be: /%/%MatrixMarket matrix coordinate real general\n";
            return;
        }
        header >> curr;
        if (curr != "real") {
            std::cerr << "Error: First line of matrix file must be: /%/%MatrixMarket matrix coordinate real general\n";
            return;
        }
        header >> curr;
        if (curr != "general") {
            std::cerr << "Error: First line of matrix file must be: /%/%MatrixMarket matrix coordinate real general\n";
            return;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error while parsing header: " << e.what() << '\n';
    }


    //skip comments
    std::getline(file, line);
    while (line[0] == '%') {
        std::getline(file, line);
    }

    //parse dimensions
    int rows, cols, nnz;
    std::istringstream metadata(line);
    metadata >> rows >> cols >> nnz;

    //temp for COO
    std::vector<int> row_indices(nnz), col_indices(nnz);
    std::vector<double> values(nnz);

    //read in COO
    for (int i = 0; i < nnz; ++i) {
        int row, col;
        double val;
        file >> row >> col >> val;

        //converts to 0 based index
        row_indices[i] = row - 1;
        col_indices[i] = col - 1;
        values[i] = val;
    }


    //initialze csr vecs
    std::vector<int> row_ptr(rows + 1, 0);
    std::vector<int> csr_col(nnz);
    std::vector<double> csr_val(nnz);

    //entries per row
    for (int i = 0; i < nnz; ++i) {
        row_ptr[row_indices[i] + 1]++;
    }

    //sum to get row
    for (int i = 0; i < rows; ++i) {
        row_ptr[i + 1] += row_ptr[i];
    }

    std::vector<int> offset = row_ptr;

    //fill in col/val arrays
    for (int i = 0; i < nnz; ++i) {
        int row = row_indices[i];
        int dest = offset[row]++;
        csr_col[dest] = col_indices[i];
        csr_val[dest] = values[i];
    }
    
    this->values = csr_val;
    this->col = csr_col;
    this->row = row_ptr;    
}

CSR_matrix::~CSR_matrix() = default;

//Converts CSR format to normal matrix. If final n cols empty, then they will not be put to the matrix
std::vector<std::vector<double>> CSR_matrix::toMatrix(){

    int num_rows = row.size()-1;
    int num_col = 0;
    for (int c: col){if (c+1 > num_col){num_col = c+1;}}
    
    std::vector<std::vector<double>> out(num_rows, std::vector<double>(num_col, 0.0));
    
    for (int i = 1; i < num_rows+1; i++){
        for (int j = row[i-1]; j < row[i]; j++){
            out[i-1][col[j]] = values[j];
        }
    }
    
    return out;
}

//perform SpMV
std::vector<double> CSR_matrix::vectorMultiply(std::vector<double>& vec) {

    std::vector<double> out(row.size() - 1, 0.0);
    for (int i = 1; i < row.size(); i++){
        for (int j = row[i-1]; j < row[i]; j++){
            out[i-1] += (values[j] * vec[col[j]]);
        }
    }
    
    return out;
}

CSR_matrix CSR_matrix::matrixMultiply(CSR_matrix& other) {

    CSR_matrix csc = other.toCSC();

    std::vector<double> out_v;
    std::vector<int> out_c;
    std::vector<int> out_r(this->row.size(), 0);

    for (int i = 0; i < this->row.size()-1; i++){
        for (int j = 0; j < csc.col.size()-1; j++){
            double res = dotHelper(csc, i, j);
            if (res != 0.0){
                out_v.push_back(res);
                out_c.push_back(j);
                out_r[i+1]++;
            }
        }
    }

    for (int i = 0; i < out_r.size()-1; i++){
        out_r[i+1] += out_r[i];
    }
    
    return CSR_matrix(out_v, out_c, out_r);
}


//creates a csc representation from a csr representin, with col indicies in col, row values in row
CSR_matrix CSR_matrix::toCSC(){

    //get max number of cols
    int num_col = 0;
    for (int c: this->col) {if (c + 1 > num_col) {num_col = c + 1;}}

    //initialze vectors for csc
    std::vector<double> v(this->values.size(), 0.0);
    std::vector<int> r(this->values.size(), 0);
    std::vector<int> c(num_col+1, 0);

    //first get counts by col and then do sum
    for (int i = 0; i < this->col.size(); i++){c[this->col[i]+1]++;}
    for (int i = 0; i < c.size()-1; i++){c[i+1] += c[i];}

    std::vector<int> counts = c;
    for (int i = 0; i < this->row.size()-1; i++){
        for (int j = this->row[i]; j < this->row[i+1]; j++){

            int col = this->col[j];
            int ptr = counts[col]++;

            v[ptr] = this->values[j];
            r[ptr] = i;
        }
    }
    
    return CSR_matrix(v, c, r);
}


double CSR_matrix::dotHelper(CSR_matrix& other, int idx1, int idx2) {

    int i = this->row[idx1];
    int j = other.col[idx2];
    double res = 0.0;    

    while (i < this->row[idx1+1] && j < other.col[idx2+1]){
        if (this->col[i] == other.row[j]){
            res += this->values[i] * other.values[j];
            i++; j++;
        }
        else if (this->col[i] < other.row[j]){i++;}
        else{j++;}
    }

    return res;
}


void CSR_matrix::print() {
    std::cout << "Values: (";
    for (size_t i = 0; i < values.size(); i++) {
        std::cout << values[i];
        if (i != values.size()-1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;

    std::cout << "Columns: (";
    for (size_t i = 0; i < col.size(); i++) {
        std::cout << col[i];
        if (i != col.size()-1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;

    std::cout << "Rows: (";
    for (size_t i = 0; i < row.size(); i++) {
        std::cout << row[i];
        if (i != row.size()-1) std::cout << ", ";
    }
    std::cout << ")" << std::endl;
    
}

