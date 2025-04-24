#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <omp.h>
#include <queue>


class CSR_matrix {
private:
    std::vector<double> values;
    std::vector<int> col;
    std::vector<int> row;

    //TODO: update constructors to support the new values;
    int n_cols;
    int n_rows;
    int nnz;

    //private structs used for merge list step of SpGEMM
    struct SpGEMM_mergelists_in{
        int* col;         // col.data()
        double* val;      // val.data()
        unsigned int idx; // 0
        int len;          // col.size()
        double scalar;

        //operator overload for pq formation
        bool operator>(const SpGEMM_mergelists_in& other) const{
            return col[idx] > other.col[other.idx];
        }
    };

    //vector of these pairs is col and val vectors of the new row
    struct SpGEMM_mergelists_out{
        int col;
        double val;
    };

public:    

    CSR_matrix();
    CSR_matrix(std::vector<std::vector<double>>& matrix);
    CSR_matrix(CSR_matrix& cpy);
    CSR_matrix(std::vector<double> v, std::vector<int> c, std::vector<int> r);
    CSR_matrix(const std::string& filename);

    ~CSR_matrix();

    std::vector<std::vector<double>> toMatrix();
    // std::vector<double> vectorMultiply(std::vector<double>& vec);
    std::vector<double> SpMV(std::vector<double>& y, std::vector<double>& x, double a, double b, int numthreads) const;
    // CSR_matrix matrixMultiply(CSR_matrix& other);
    std::vector<std::pair<int, double>> tuple_mergeSorted(std::vector<std::vector<std::pair<int, double>>>);
    std::vector<SpGEMM_mergelists_out> addr_mergeSorted(std::vector<SpGEMM_mergelists_in> vecs);
    CSR_matrix SpGEMM(CSR_matrix& B);

    void print();
    
    // double dotHelper( CSR_matrix& other, int idx1, int idx)
    // CSR_matrix toCSC();

    // Getters
    const std::vector<double>& getValues() const { return values; }
    const std::vector<int>& getCol() const { return col; }
    const std::vector<int>& getRow() const { return row; }
    int getNCols() const { return n_cols; }
    int getNRows() const { return n_rows; }
    int getNNZ() const { return nnz; }
    
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


    //TODO: update these values properly
    n_cols = 0;
    n_rows = 0;
    nnz = 0;
    
}

CSR_matrix::CSR_matrix(CSR_matrix& cpy): values(cpy.values), col(cpy.col), row(cpy.row), n_cols(cpy.n_cols), n_rows(cpy.n_rows), nnz(cpy.nnz) {}

CSR_matrix::CSR_matrix(std::vector<double> v, std::vector<int> c, std::vector<int> r): values(v), col(c), row(r) {
    
    //TODO: update these values properly
    n_cols = 0;
    n_rows = 0;
    nnz = 0;
}


CSR_matrix::CSR_matrix(const std::string& filename){

    //support special filename, ends in .mtxtxt
    //first row is: num_col num_row nnz
    //then the next nnz rows represent the val vector, then col vector, then the row+1 represnt row ptr

    if (filename.substr(filename.find_last_of(".") + 1) == "mtxtxt")
    {
        //load file
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open the file " << filename << "\n";
            return;
        }

        std::string line;
        std::getline(file, line);
        std::stringstream str(line);
        str >> n_cols >> n_rows >> nnz;

        double d;
        int in;

        values.resize(nnz);
        col.resize(nnz);
        row.resize(n_rows+1);

        for (size_t i = 0; i < nnz; ++i)
        {
            file >> d;
            values.at(i) = d;
        }
        
        for (size_t i = 0; i < nnz; ++i)
        {
            file >> in;
            col.at(i) = in;
        }

        for (size_t i = 0; i < n_rows+1; ++i)
        {
            file >> in;
            row.at(i) = in;
        }
        file.close();
    }
    else /* .mtx format */ {

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
                std::cerr << "Error in file " << filename << ": First line of matrix file must be: /%/%MatrixMarket matrix coordinate real general\n";
                std::cerr << "Actual: " << curr << "\n";
                return;
            }
            header >> curr;
            if (curr != "coordinate") {
                std::cerr << "Error in file " << filename << ": First line of matrix file must be: /%/%MatrixMarket matrix coordinate real general\n";
                std::cerr << "Actual: " << curr << "\n";
                return;
            }
            header >> curr;
            if (curr != "real") {
                std::cerr << "Error in file " << filename << ": First line of matrix file must be: /%/%MatrixMarket matrix coordinate real general\n";
                std::cerr << "Actual: " << curr << "\n";
                return;
            }
            header >> curr;
            if (curr != "general") {
                std::cerr << "Error in file " << filename << ": First line of matrix file must be: /%/%MatrixMarket matrix coordinate real general\n";
                std::cerr << "Actual: " << curr << "\n";
                return;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error in file " << filename << " while parsing header: " << e.what() << '\n';
        }


        //skip comments
        std::getline(file, line);
        while (line[0] == '%') {
            std::getline(file, line);
        }

        //parse dimensions
        int rows, cols, nonzero;
        std::istringstream metadata(line);
        metadata >> rows >> cols >> nonzero;

        n_cols = cols;
        n_rows = rows;
        nnz = nonzero;

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
        file.close();
    }
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

//perform SpMV of form y = α · Ax + β · y
//where y is the output vector, α and β are scalars, x is the dense vector, and A is the sparse matrix 
std::vector<double> CSR_matrix::SpMV(std::vector<double>& y, std::vector<double>& x, double a, double b, int numthreads) const {
    
    std::vector<double> out(y.size(), 0.0);

    #pragma omp parallel for num_threads(numthreads)
    for (int i = 0; i < n_rows; i++){
        double rowSum = 0.0;
        //calculate Ax
        #pragma omp simd reduction(+:rowSum)
        for (int j = row[i]; j < row[i + 1]; j++) {
            rowSum += values[j] * x[col[j]];
        }
        
        //add a * Ax + b * y to the output vector
        out[i] = a * rowSum + b * y[i];
    }

    return out;
}

//could try to use struct instead of tuple
std::vector<std::pair<int, double>> CSR_matrix::tuple_mergeSorted(std::vector<std::vector<std::pair<int, double>>> vecs){

    using triple = std::tuple<std::pair<int, double>, int, int>;

    auto comparator = [](const triple& a, const triple& b){
        return std::get<0>(a).first > std::get<0>(b).first;
    };

    std::vector<triple> temp(vecs.size());

    for (size_t i = 0; i < vecs.size(); i++){
        temp.at(i) = {vecs[i][0], i, 0};
    }

    std::priority_queue<triple, std::vector<triple>, decltype(comparator)> pq(comparator, std::move(temp));

    std::vector<std::pair<int, double>> res(vecs.size() * vecs.size());

    while (!pq.empty())
    {
        triple element = pq.top();
        pq.pop();

        std::pair<int, double> val = std::get<0>(element);
        int global_i = std::get<1>(element);
        int vector_i = std::get<2>(element);

        vector_i++;

        if (vector_i < vecs[global_i].size()){
            pq.push({vecs[global_i][vector_i], global_i, vector_i});
        }

        if (res.empty() || res.at(res.size()-1).first != val.first){
            res.emplace_back(val);
        }
        else{
            res[res.size()-1].second += val.second;
        }
    }

    return res;
}


//try to use memory addresses with structs
//vecs is an array of pointer pairs, each indicating the start of a vector
std::vector<CSR_matrix::SpGEMM_mergelists_out> CSR_matrix::addr_mergeSorted(std::vector<SpGEMM_mergelists_in> vecs){

    std::priority_queue<SpGEMM_mergelists_in, std::vector<SpGEMM_mergelists_in>, std::greater<SpGEMM_mergelists_in>> pq(std::greater<SpGEMM_mergelists_in>(), std::move(vecs));
    
    std::vector<SpGEMM_mergelists_out> out;
    out.reserve(vecs.size() * vecs.size());

    while (!pq.empty())
    {
        SpGEMM_mergelists_in r = pq.top();
        pq.pop();

        //if the output vector is empty or if the column of current item is not already present
        if (out.empty() || out[out.size() - 1].col != r.col[r.idx])
        {
            out.emplace_back(SpGEMM_mergelists_out{r.col[r.idx], r.val[r.idx] * r.scalar});
        }
        //if we need to add to the value of the greatest column in the output vector
        else{
            out[out.size() - 1].val += r.val[r.idx] * r.scalar;
        }

        ++r.idx;
        //emplace if the size isnt more than the length of the vector
        if (r.idx < r.len)
        {
            pq.emplace(SpGEMM_mergelists_in{r.col, r.val, r.idx, r.len, r.scalar});
        }
    }

    return out;
    
}

CSR_matrix CSR_matrix::SpGEMM(CSR_matrix& B){

    std::vector<double> new_v;
    std::vector<int> new_c;
    std::vector<int> new_r(n_rows + 1);
    new_v.reserve(nnz);
    new_c.reserve(nnz);

    std::vector<SpGEMM_mergelists_in> vec;
    std::vector<SpGEMM_mergelists_out> out;
    vec.reserve(n_cols);
    out.reserve(n_cols);

    for (size_t i = 0; i < n_rows; i++)
    {
        vec.clear();
        out.clear();
        for (size_t j = row[i]; j < row[i+1]; j++)
        {
            vec.emplace_back(SpGEMM_mergelists_in{&B.col[B.row[col[j]]], &B.values[B.row[col[j]]], 0, B.row[col[j] + 1 ] - B.row[col[j]], values[j]});
        }
        out = addr_mergeSorted(vec);

        for (size_t k = 0; k < out.size(); k++)
        {
            SpGEMM_mergelists_out item = out[k];
            new_c.emplace_back(item.col);
            new_v.emplace_back(item.val);
            ++new_r[i];
        }
    }

    for (int i = 0; i < n_rows; ++i) {
        new_r[i + 1] += new_r[i];
    }

    return CSR_matrix(new_v, new_c, new_r);
}


// std::vector<double> CSR_matrix::vectorMultiply(std::vector<double>& vec) {

//     std::vector<double> out(row.size() - 1, 0.0);
//     for (int i = 1; i < row.size(); i++){
//         for (int j = row[i-1]; j < row[i]; j++){
//             out[i-1] += (values[j] * vec[col[j]]);
//         }
//     }
    
//     return out;
// }




// CSR_matrix CSR_matrix::matrixMultiply(CSR_matrix& other) {

//     CSR_matrix csc = other.toCSC();

//     std::vector<double> out_v;
//     std::vector<int> out_c;
//     std::vector<int> out_r(this->row.size(), 0);

//     for (int i = 0; i < this->row.size()-1; i++){
//         for (int j = 0; j < csc.col.size()-1; j++){
//             double res = dotHelper(csc, i, j);
//             if (res != 0.0){
//                 out_v.push_back(res);
//                 out_c.push_back(j);
//                 out_r[i+1]++;
//             }
//         }
//     }

//     for (int i = 0; i < out_r.size()-1; i++){
//         out_r[i+1] += out_r[i];
//     }
    
//     return CSR_matrix(out_v, out_c, out_r);
// }


// //creates a csc representation from a csr representin, with col indicies in col, row values in row
// CSR_matrix CSR_matrix::toCSC(){

//     //get max number of cols
//     int num_col = 0;
//     for (int c: this->col) {if (c + 1 > num_col) {num_col = c + 1;}}

//     //initialze vectors for csc
//     std::vector<double> v(this->values.size(), 0.0);
//     std::vector<int> r(this->values.size(), 0);
//     std::vector<int> c(num_col+1, 0);

//     //first get counts by col and then do sum
//     for (int i = 0; i < this->col.size(); i++){c[this->col[i]+1]++;}
//     for (int i = 0; i < c.size()-1; i++){c[i+1] += c[i];}

//     std::vector<int> counts = c;
//     for (int i = 0; i < this->row.size()-1; i++){
//         for (int j = this->row[i]; j < this->row[i+1]; j++){

//             int col = this->col[j];
//             int ptr = counts[col]++;

//             v[ptr] = this->values[j];
//             r[ptr] = i;
//         }
//     }
    
//     return CSR_matrix(v, c, r);
// }


// double CSR_matrix::dotHelper(CSR_matrix& other, int idx1, int idx2) {

//     int i = this->row[idx1];
//     int j = other.col[idx2];
//     double res = 0.0;    

//     while (i < this->row[idx1+1] && j < other.col[idx2+1]){
//         if (this->col[i] == other.row[j]){
//             res += this->values[i] * other.values[j];
//             i++; j++;
//         }
//         else if (this->col[i] < other.row[j]){i++;}
//         else{j++;}
//     }

//     return res;
// }


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

