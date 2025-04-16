The files contained in MM_matrix_files are from NIST available at https://math.nist.gov/MatrixMarket/data/.

To run, first insure you have the linux version of Intel's OpenMKL (available here: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) installed and do the following:

```bash
# Set up environment
source /opt/intel/oneapi/setvars.sh

# Build the binary
make

# Run the example
./SpMV