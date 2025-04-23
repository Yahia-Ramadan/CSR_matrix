#include "CSR_matrix.hpp"

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

int main() {

    //Load matrix
    CSR_matrix a("MM_matrix_files/west0067.mtx");

    //Generate random dense vector with size 67
    float min_value = -100.0;
    float max_value = 100.0;

    std::vector<double> x(67, 0.0);
    std::vector<double> y(67, 0.0);

    for (size_t i = 0; i < 67; i++)
    {
        x.at(i) = generate_random_float(min_value, max_value);
    }

    //Multiply

    std::vector<double> out;
    out = a.SpMV(y, x, 2.5, 0.0);

    // Write the matrix, input vector, and output vector to text files
    std::ofstream input_vector_file("input_vector.txt");
    std::ofstream output_vector_file("output_vector.txt");
    std::ofstream row_file("matrix_rows.txt");
    std::ofstream col_file("matrix_cols.txt");
    std::ofstream value_file("matrix_values.txt");

    const auto& rows = a.getRow();
    const auto& cols = a.getCol();
    const auto& values = a.getValues();


    if (row_file.is_open()){
        for (const auto& val: a.getRow())
        {
            row_file << val << "\n";
        }
        row_file.close();
    }

    if (col_file.is_open()){
        for (const auto& val: a.getCol())
        {
            col_file << val << "\n";
        }
        col_file.close();
    }

    if (value_file.is_open()){
        for (const auto& val: a.getValues())
        {
            value_file << val << "\n";
        }
        value_file.close();
    }


    if (input_vector_file.is_open()) {
        for (const auto& val : x) {
            input_vector_file << val << "\n";
        }
        input_vector_file.close();
    } else {
        std::cerr << "Unable to open file for writing input vector." << std::endl;
    }

    if (output_vector_file.is_open()) {
        for (const auto& val : out) {
            output_vector_file << val << "\n";
        }
        output_vector_file.close();
    } else {
        std::cerr << "Unable to open file for writing output vector." << std::endl;
    }

}