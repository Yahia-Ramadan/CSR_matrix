# Path to Intel MKL's mkl_link_tool
MKL_LINK_TOOL = /opt/intel/oneapi/mkl/latest/bin/mkl_link_tool

# C++ compiler
CXX = g++
CXXFLAGS = -O2 -std=c++17 -fopenmp -ftree-vectorize -march=native

# Output binary name
TARGET = SpMV

# Source file
SRC = SpMV_testing_threaded.cpp

# Use mkl_link_tool to get MKL linker flags
MKL_FLAGS = $(shell $(MKL_LINK_TOOL) --compiler=gnu_c --mkl=threading -libs)

# Additional system libs needed by MKL
SYS_LIBS = -lpthread -lm -ldl

# Default rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET) $(MKL_FLAGS) $(SYS_LIBS)

clean:
	rm -f $(TARGET)
