Some of the files contained in MM_matrix_files are from NIST and available at https://math.nist.gov/MatrixMarket/data/. Others are from the SuiteSparse Matrix Collection and available at https://sparse.tamu.edu/.

To run, first change the SRC in the makefile to be the file that you'd like to run and ensure you have the Linux version of Intel's OpenMKL (available here: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) installed and do the following:

To Install Eigen: sudo apt-get install libeigen3-dev

```bash
# Set up environment
source /opt/intel/oneapi/setvars.sh

# Build the binary
make

# Run the example
./SpMV

