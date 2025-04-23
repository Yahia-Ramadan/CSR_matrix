import numpy as np
from scipy.sparse import csr_matrix

# Example CSR matrix

# Read data, indices, and indptr from files
v = np.loadtxt('matrix_values.txt', dtype=np.float128)
c = np.loadtxt('matrix_cols.txt', dtype=int)
r = np.loadtxt('matrix_rows.txt', dtype=int)

sparse_matrix = csr_matrix((v, c, r), shape=(67, 67))

# Dense vector
x = np.loadtxt('input_vector.txt')

# Perform Sparse Matrix-Vector Multiplication (SpMV)
result = sparse_matrix.dot(x)

# Multiply the result by a scalar
scalar = 2.5
result = result * scalar

# Load the comparison vector from a text file
comparison_vector = np.loadtxt('output_vector.txt')

# Compare each value in the result to the corresponding value in the comparison vector
# Allow for a tolerance in floating point comparison
tolerance = 1e-4
comparison = np.abs(result - comparison_vector) <= tolerance

# Print mismatched indices and values
mismatched_indices = np.where(comparison == False)[0]
for idx in mismatched_indices:
    print(f"Mismatch at index {idx}: result={result[idx]}, comparison_vector={comparison_vector[idx]}")
