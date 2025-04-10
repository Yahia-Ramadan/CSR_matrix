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


//this constructor is still buggy
CSR_matrix::CSR_matrix(const std::string& filename) {

    //parse the comment lines
    std::ifstream file(filename);
    std::string line;

    std::getline(file, line);
    while (line[0] == '%'){
        std::getline(file, line);
    }

    //now string should contain dimension lines
    int num_row;
    int num_col;
    int num_entries;

    std::istringstream iss(line);
    iss >> num_row >> num_col >> num_entries;

    // std::cout << num_row << " " << num_col << "\n";
    
    //make the matrix with zeros
    std::vector<std::vector<double>> out(num_row, std::vector<double>(num_col, 0.0));

    //popualte
    int in_v, in_c, in_r;
    while (file >> in_r >> in_c >> in_v){
        // std::cout << "\nitem\n";
        out[in_r-1][in_c-1] = in_v;
    }

    // for (const auto& row : out) {
    //     for (const auto& val : row) {
    //         std::cout << val << " ";
    //     }
    //     std::cout << std::endl;
    // }

    *this = CSR_matrix(out);    
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

//multiply a CSR matrix with a column vector
std::vector<double> CSR_matrix::vectorMultiply(std::vector<double>& vec) {

    std::vector<double> out(vec.size(), 0.0);
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

