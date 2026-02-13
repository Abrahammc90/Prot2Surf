# Compiler and flags
FC = gfortran
NVCC = nvcc
#FLAGS = -O0 -Wall -fno-fast-math
FLAGS = -O0 -Wall -llapack -lblas -fopenmp -march=native -funroll-loops -fno-fast-math
CUDA_FLAGS = -O2 -arch=sm_86 -std=c++11

# Source files
MOD_PDB_SRC = mod_pdb.f90
MOD_ASSOC_SRC = mod_assoc.f90
READ_INPUT_SRC = read_input.f90
MATHS_SRC = maths.f90
MOD_CUDA_SRC = mod_cuda.f90
CUDA_KERNELS_SRC = cuda_kernels.cu
MOD_MATRIX_SRC = mod_matrix.f90
MOD_CLUST_ALGORITHM_SRC = mod_clust_algorithm.f90
MOD_THRESHOLD_SRC = mod_threshold.f90
MAKE_MATRIX_SRC = make_matrix.f90
ANALYZE_RESIDUES_SRC = analyze_residues.f90
CLUST_SRC = clust.f90
CLUST_ALL_SRC = clust_all.f90
THRESH_SRC = threshold.f90

# Object files
MOD_PDB_OBJ = mod_pdb.o
MOD_ASSOC_OBJ = mod_assoc.o
READ_INPUT_OBJ = read_input.o
MATHS_OBJ = maths.o
MOD_CUDA_OBJ = mod_cuda.o
CUDA_KERNELS_OBJ = cuda_kernels.o
MOD_MATRIX_OBJ = mod_matrix.o
MOD_CLUST_ALGORITHM_OBJ = mod_clust_algorithm.o
MOD_THRESHOLD_OBJ = mod_threshold.o
MAKE_MATRIX_OBJ = make_matrix.o
ANALYZE_RESIDUES_OBJ = analyze_residues.o
CLUST_OBJ = clust.o
CLUST_ALL_OBJ = clust_all.o
THRESH_OBJ = threshold.o

# Executable
BIN_DIR = ../bin
EXE1 = $(BIN_DIR)/make_matrix
EXE2 = $(BIN_DIR)/clust
EXE3 = $(BIN_DIR)/analyze_residues
EXE4 = $(BIN_DIR)/threshold
EXE5 = $(BIN_DIR)/clust_all

# Default target: build the executable and documentation
all: $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) docs

# Rule to compile mod_pdb
$(MOD_PDB_OBJ): $(MOD_PDB_SRC)
	$(FC) $(FLAGS) -c $(MOD_PDB_SRC)

# Rule to compile mod_assoc
$(MOD_ASSOC_OBJ): $(MOD_ASSOC_SRC)
	$(FC) $(FLAGS) -c $(MOD_ASSOC_SRC)

# Rule to compile CUDA kernels
$(CUDA_KERNELS_OBJ): $(CUDA_KERNELS_SRC)
	$(NVCC) $(CUDA_FLAGS) -c $(CUDA_KERNELS_SRC) -o $(CUDA_KERNELS_OBJ)

# Rule to compile mod_cuda (depends on CUDA kernels)
$(MOD_CUDA_OBJ): $(MOD_CUDA_SRC) $(CUDA_KERNELS_OBJ)
	$(FC) $(FLAGS) -c $(MOD_CUDA_SRC)

# Rule to compile read_input (depends on mod_pdb and mod_assoc)
$(READ_INPUT_OBJ): $(READ_INPUT_SRC) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -c $(READ_INPUT_SRC)

# Rule to compile mod_clust_algorithm (depends on read_input and mod_cuda)
$(MOD_CLUST_ALGORITHM_OBJ): $(MOD_CLUST_ALGORITHM_SRC) $(READ_INPUT_OBJ) $(MOD_CUDA_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -c $(MOD_CLUST_ALGORITHM_SRC)

# Rule to compile maths
$(MATHS_OBJ): $(MATHS_SRC)
	$(FC) $(FLAGS) -c $(MATHS_SRC)



# Rule to compile mod_matrix (depends on maths)
$(MOD_MATRIX_OBJ): $(MOD_MATRIX_SRC) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(MOD_MATRIX_SRC)

# Rule to compile mod_threshold (depends on maths)
$(MOD_THRESHOLD_OBJ): $(MOD_THRESHOLD_SRC) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(MOD_THRESHOLD_SRC)



# Rule to compile make_matrix (depends on everything else)
$(MAKE_MATRIX_OBJ): $(MAKE_MATRIX_SRC) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_MATRIX_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(MAKE_MATRIX_SRC)

# Rule to compile clust (depends on mod_matrix and mod_clust_algorithm)
$(CLUST_OBJ): $(CLUST_SRC) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_MATRIX_OBJ) $(MOD_CUDA_OBJ) $(CUDA_KERNELS_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(CLUST_SRC)

# Rule to compile clust_all (depends on mod_matrix and mod_clust_algorithm)
$(CLUST_ALL_OBJ): $(CLUST_ALL_SRC) $(MOD_MATRIX_OBJ) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_CUDA_OBJ) $(CUDA_KERNELS_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(CLUST_ALL_SRC)

# Rule to compile analyze_residues (depends on read_input and maths)
$(ANALYZE_RESIDUES_OBJ): $(ANALYZE_RESIDUES_SRC) $(MATHS_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -c $(ANALYZE_RESIDUES_SRC)

# Rule to compile threshold (depends on everything else)
$(THRESH_OBJ): $(THRESH_SRC) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_THRESHOLD_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(THRESH_SRC)

# Rule to link the make_matrix executable (depends on all object files)
$(EXE1): $(MAKE_MATRIX_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_MATRIX_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -o $(EXE1) $(MAKE_MATRIX_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_MATRIX_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)

# Rule to link the clust executable (depends on mod_matrix and mod_clust_algorithm)
$(EXE2): $(CLUST_OBJ) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_MATRIX_OBJ) $(MOD_CUDA_OBJ) $(CUDA_KERNELS_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -o $(EXE2) $(CLUST_OBJ) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_MATRIX_OBJ) $(MATHS_OBJ) $(MOD_CUDA_OBJ) $(CUDA_KERNELS_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) -L/usr/local/cuda/lib64 -lcudart -lstdc++

# Rule to link the analyze_residues executable (depends on mod_matrix and mod_clust_algorithm)
$(EXE3): $(ANALYZE_RESIDUES_OBJ) $(MATHS_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -o $(EXE3) $(ANALYZE_RESIDUES_OBJ) $(MATHS_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)

# Rule to link the threshold executable (depends on all object files)
$(EXE4): $(THRESH_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_THRESHOLD_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -o $(EXE4) $(THRESH_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_THRESHOLD_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)

# Rule to link the clust_all executable (depends on mod_matrix and mod_clust_algorithm)
$(EXE5): $(CLUST_ALL_OBJ) $(MOD_MATRIX_OBJ) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_CUDA_OBJ) $(CUDA_KERNELS_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -o $(EXE5) $(CLUST_ALL_OBJ) $(MOD_MATRIX_OBJ) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MATHS_OBJ) $(MOD_CUDA_OBJ) $(CUDA_KERNELS_OBJ) -L/usr/local/cuda/lib64 -lcudart -lstdc++

# Unit test targets - Fortran tests
UNIT_TESTS_DIR = unit_test
TEST_MATHS = $(UNIT_TESTS_DIR)/test_maths
TEST_MATRIX = $(UNIT_TESTS_DIR)/test_matrix
TEST_THRESHOLD = $(UNIT_TESTS_DIR)/test_threshold
TEST_PDB = $(UNIT_TESTS_DIR)/test_pdb
TEST_CLUSTERING = $(UNIT_TESTS_DIR)/test_clustering
TEST_ASSOC = $(UNIT_TESTS_DIR)/test_assoc
TEST_READ_INPUT = $(UNIT_TESTS_DIR)/test_read_input

# Build and run all unit tests (Fortran only in src)
test: fortran_tests
	@echo ""
	@echo "========================================="
	@echo "All Fortran unit tests completed!"
	@echo "========================================="
	@echo "Note: Python tests are in ../tools/"

# Run only Fortran unit tests
fortran_tests: test_maths test_matrix test_threshold test_pdb test_clustering test_assoc test_read_input
	@echo ""
	@echo "All Fortran unit tests completed successfully!"
	@echo ""

# Test mathematical operations
test_maths: $(MOD_PDB_OBJ) $(UNIT_TESTS_DIR)/test_maths.f90
	$(FC) $(FLAGS) -o $(TEST_MATHS) $(UNIT_TESTS_DIR)/test_maths.f90 $(MATHS_OBJ) $(MOD_PDB_OBJ)
	$(TEST_MATHS)

# Test matrix operations
test_matrix: $(MOD_PDB_OBJ) $(MATHS_OBJ) $(UNIT_TESTS_DIR)/test_matrix.f90
	$(FC) $(FLAGS) -o $(TEST_MATRIX) $(UNIT_TESTS_DIR)/test_matrix.f90 $(MOD_MATRIX_OBJ) $(MATHS_OBJ) $(MOD_PDB_OBJ)
	$(TEST_MATRIX)

# Test threshold operations
test_threshold: $(UNIT_TESTS_DIR)/test_threshold.f90
	$(FC) $(FLAGS) -o $(TEST_THRESHOLD) $(UNIT_TESTS_DIR)/test_threshold.f90 $(MOD_THRESHOLD_OBJ) $(MATHS_OBJ)
	$(TEST_THRESHOLD)

# Test PDB operations
test_pdb: $(MOD_PDB_OBJ) $(UNIT_TESTS_DIR)/test_pdb.f90
	$(FC) $(FLAGS) -o $(TEST_PDB) $(UNIT_TESTS_DIR)/test_pdb.f90 $(MOD_PDB_OBJ)
	$(TEST_PDB)

# Test clustering operations
test_clustering: $(UNIT_TESTS_DIR)/test_clustering.f90
	$(FC) $(FLAGS) -o $(TEST_CLUSTERING) $(UNIT_TESTS_DIR)/test_clustering.f90 $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MATHS_OBJ)
	$(TEST_CLUSTERING)

# Test association file reading
test_assoc: $(MOD_ASSOC_OBJ) $(UNIT_TESTS_DIR)/test_assoc.f90
	$(FC) $(FLAGS) -o $(TEST_ASSOC) $(UNIT_TESTS_DIR)/test_assoc.f90 $(MOD_ASSOC_OBJ)
	$(TEST_ASSOC)

# Test read_input module
test_read_input: $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(UNIT_TESTS_DIR)/test_read_input.f90
	$(FC) $(FLAGS) -o $(TEST_READ_INPUT) $(UNIT_TESTS_DIR)/test_read_input.f90 $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(TEST_READ_INPUT)

# Documentation target: generate Doxygen documentation
docs:
	doxygen Doxyfile
	cd ../docs/latex && make pdf && cd ../../src

# Clean target: remove all compiled files (including test executables)
clean:
	rm -f $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) *.o *.mod
	rm -f $(TEST_MATHS) $(TEST_MATRIX) $(TEST_THRESHOLD) $(TEST_PDB) $(TEST_CLUSTERING) $(TEST_ASSOC) $(TEST_READ_INPUT)
	rm -rf ../docs/doxygen

# Clean only Fortran test files
clean_tests:
	rm -f $(TEST_MATHS) $(TEST_MATRIX) $(TEST_THRESHOLD) $(TEST_PDB) $(TEST_CLUSTERING) $(TEST_ASSOC) $(TEST_READ_INPUT)
