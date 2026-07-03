# Compiler and flags
FC = gfortran
#FLAGS = -O0 -Wall -fno-fast-math
FLAGS = -O0 -Wall -llapack -lblas -fopenmp -march=native -funroll-loops -fno-fast-math

# Source files
MOD_PDB_SRC = mod_pdb.f90
MOD_ASSOC_SRC = mod_assoc.f90
READ_INPUT_SRC = read_input.f90
MATHS_SRC = maths.f90
MOD_ARRAY_SRC = mod_array.f90
MOD_CLUST_ALGORITHM_SRC = mod_clust_algorithm.f90
MOD_THRESHOLD_SRC = mod_threshold.f90
MAKE_DATA_SRC = make_data.f90
CONTACT_MAP_SRC = contact_map.f90
CLUST_SRC = clust.f90
THRESH_SRC = threshold.f90

# Object files
MOD_PDB_OBJ = mod_pdb.o
MOD_ASSOC_OBJ = mod_assoc.o
READ_INPUT_OBJ = read_input.o
MATHS_OBJ = maths.o
MOD_ARRAY_OBJ = mod_array.o
MOD_CLUST_ALGORITHM_OBJ = mod_clust_algorithm.o
MOD_THRESHOLD_OBJ = mod_threshold.o
MAKE_DATA_OBJ = make_data.o
CONTACT_MAP_OBJ = contact_map.o
CLUST_OBJ = clust.o
THRESH_OBJ = threshold.o

# Executable
BIN_DIR = ../bin
EXE1 = $(BIN_DIR)/make_data
EXE2 = $(BIN_DIR)/clust
EXE3 = $(BIN_DIR)/contact_map
EXE4 = $(BIN_DIR)/threshold

# Default target: build executables, test binaries, and documentation
all: $(EXE1) $(EXE2) $(EXE3) $(EXE4) build_fortran_tests docs

# Rule to compile mod_pdb
$(MOD_PDB_OBJ): $(MOD_PDB_SRC)
	$(FC) $(FLAGS) -c $(MOD_PDB_SRC)

# Rule to compile mod_assoc
$(MOD_ASSOC_OBJ): $(MOD_ASSOC_SRC)
	$(FC) $(FLAGS) -c $(MOD_ASSOC_SRC)

# Rule to compile read_input (depends on mod_pdb and mod_assoc)
$(READ_INPUT_OBJ): $(READ_INPUT_SRC) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -c $(READ_INPUT_SRC)

# Rule to compile mod_clust_algorithm (depends on read_input)
$(MOD_CLUST_ALGORITHM_OBJ): $(MOD_CLUST_ALGORITHM_SRC) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -c $(MOD_CLUST_ALGORITHM_SRC)

# Rule to compile maths
$(MATHS_OBJ): $(MATHS_SRC)
	$(FC) $(FLAGS) -c $(MATHS_SRC)

# Rule to compile mod_array (depends on maths and mod_assoc)
$(MOD_ARRAY_OBJ): $(MOD_ARRAY_SRC) $(MATHS_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -c $(MOD_ARRAY_SRC)

# Rule to compile mod_threshold (depends on maths)
$(MOD_THRESHOLD_OBJ): $(MOD_THRESHOLD_SRC) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(MOD_THRESHOLD_SRC)



# Rule to compile make_data (depends on array/data input modules)
$(MAKE_DATA_OBJ): $(MAKE_DATA_SRC) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_ARRAY_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(MAKE_DATA_SRC)

# Rule to compile clust (depends on mod_array and mod_clust_algorithm)
$(CLUST_OBJ): $(CLUST_SRC) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_ARRAY_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(CLUST_SRC)

# Rule to compile contact_map (depends on read_input and maths)
$(CONTACT_MAP_OBJ): $(CONTACT_MAP_SRC) $(MATHS_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -c $(CONTACT_MAP_SRC)

# Rule to compile threshold
$(THRESH_OBJ): $(THRESH_SRC) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_THRESHOLD_OBJ) $(MOD_ARRAY_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(THRESH_SRC)

# Rule to link the make_data executable
$(EXE1): $(MAKE_DATA_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_ARRAY_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -o $(EXE1) $(MAKE_DATA_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_ARRAY_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)

# Rule to link the clust executable
$(EXE2): $(CLUST_OBJ) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_ARRAY_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -o $(EXE2) $(CLUST_OBJ) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_ARRAY_OBJ) $(MATHS_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)

# Rule to link the contact_map executable
$(EXE3): $(CONTACT_MAP_OBJ) $(MATHS_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -o $(EXE3) $(CONTACT_MAP_OBJ) $(MATHS_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)

# Rule to link the threshold executable (depends on all object files)
$(EXE4): $(THRESH_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_THRESHOLD_OBJ) $(MOD_ARRAY_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -o $(EXE4) $(THRESH_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_THRESHOLD_OBJ) $(MOD_ARRAY_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)

# Unit test targets - Fortran tests
UNIT_TESTS_DIR = unit_test
TEST_MATHS = $(UNIT_TESTS_DIR)/test_maths
TEST_ARRAY = $(UNIT_TESTS_DIR)/test_array
TEST_THRESHOLD = $(UNIT_TESTS_DIR)/test_threshold
TEST_PDB = $(UNIT_TESTS_DIR)/test_pdb
TEST_CLUSTERING = $(UNIT_TESTS_DIR)/test_clustering
TEST_ASSOC = $(UNIT_TESTS_DIR)/test_assoc
TEST_READ_INPUT = $(UNIT_TESTS_DIR)/test_read_input

# Build and run all unit tests (Fortran only in src)
test: fortran_tests
	@$(MAKE) cleanup_test_artifacts
	@echo ""
	@echo "========================================="
	@echo "All Fortran unit tests completed!"
	@echo "========================================="
	@echo "Note: Python tests are in ../tools/"

# Build Fortran test binaries only (no execution)
build_fortran_tests: $(TEST_MATHS) $(TEST_ARRAY) $(TEST_THRESHOLD) $(TEST_PDB) $(TEST_CLUSTERING) $(TEST_ASSOC) $(TEST_READ_INPUT)
	@echo ""
	@echo "Fortran test binaries built successfully!"
	@echo ""

# Run only Fortran unit tests
fortran_tests: test_maths test_array test_threshold test_pdb test_clustering test_assoc test_read_input
	@echo ""
	@echo "All Fortran unit tests completed successfully!"
	@echo ""

# Build test_maths binary
$(TEST_MATHS): $(MOD_PDB_OBJ) $(MATHS_OBJ) $(UNIT_TESTS_DIR)/test_maths.f90
	$(FC) $(FLAGS) -o $(TEST_MATHS) $(UNIT_TESTS_DIR)/test_maths.f90 $(MATHS_OBJ) $(MOD_PDB_OBJ)

# Test mathematical operations
test_maths: $(MOD_PDB_OBJ) $(MATHS_OBJ) $(UNIT_TESTS_DIR)/test_maths.f90
	@$(MAKE) $(TEST_MATHS)
	$(TEST_MATHS)

# Build test_array binary
$(TEST_ARRAY): $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MATHS_OBJ) $(MOD_ARRAY_OBJ) $(UNIT_TESTS_DIR)/test_array.f90
	$(FC) $(FLAGS) -o $(TEST_ARRAY) $(UNIT_TESTS_DIR)/test_array.f90 $(MOD_ARRAY_OBJ) $(MATHS_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)

# Test array operations
test_array: $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MATHS_OBJ) $(MOD_ARRAY_OBJ) $(UNIT_TESTS_DIR)/test_array.f90
	@$(MAKE) $(TEST_ARRAY)
	$(TEST_ARRAY)

# Build test_threshold binary
$(TEST_THRESHOLD): $(MOD_THRESHOLD_OBJ) $(MATHS_OBJ) $(UNIT_TESTS_DIR)/test_threshold.f90
	$(FC) $(FLAGS) -o $(TEST_THRESHOLD) $(UNIT_TESTS_DIR)/test_threshold.f90 $(MOD_THRESHOLD_OBJ) $(MATHS_OBJ)

# Test threshold operations
test_threshold: $(MOD_THRESHOLD_OBJ) $(MATHS_OBJ) $(UNIT_TESTS_DIR)/test_threshold.f90
	@$(MAKE) $(TEST_THRESHOLD)
	$(TEST_THRESHOLD)

# Build test_pdb binary
$(TEST_PDB): $(MOD_PDB_OBJ) $(UNIT_TESTS_DIR)/test_pdb.f90
	$(FC) $(FLAGS) -o $(TEST_PDB) $(UNIT_TESTS_DIR)/test_pdb.f90 $(MOD_PDB_OBJ)

# Test PDB operations
test_pdb: $(MOD_PDB_OBJ) $(UNIT_TESTS_DIR)/test_pdb.f90
	@$(MAKE) $(TEST_PDB)
	$(TEST_PDB)

# Build test_clustering binary
$(TEST_CLUSTERING): $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MATHS_OBJ) $(UNIT_TESTS_DIR)/test_clustering.f90
	$(FC) $(FLAGS) -o $(TEST_CLUSTERING) $(UNIT_TESTS_DIR)/test_clustering.f90 $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MATHS_OBJ)

# Test clustering operations
test_clustering: $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MATHS_OBJ) $(UNIT_TESTS_DIR)/test_clustering.f90
	@$(MAKE) $(TEST_CLUSTERING)
	$(TEST_CLUSTERING)

# Build test_assoc binary
$(TEST_ASSOC): $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(UNIT_TESTS_DIR)/test_assoc.f90
	$(FC) $(FLAGS) -o $(TEST_ASSOC) $(UNIT_TESTS_DIR)/test_assoc.f90 $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)

# Test association file reading
test_assoc: $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(UNIT_TESTS_DIR)/test_assoc.f90
	@$(MAKE) $(TEST_ASSOC)
	$(TEST_ASSOC)

# Build test_read_input binary
$(TEST_READ_INPUT): $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(UNIT_TESTS_DIR)/test_read_input.f90
	$(FC) $(FLAGS) -o $(TEST_READ_INPUT) $(UNIT_TESTS_DIR)/test_read_input.f90 $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)

# Test read_input module
test_read_input: $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(UNIT_TESTS_DIR)/test_read_input.f90
	@$(MAKE) $(TEST_READ_INPUT)
	$(TEST_READ_INPUT)

# Documentation target: generate Doxygen documentation
docs:
	@if command -v pdflatex >/dev/null 2>&1; then \
		echo "[docs] pdflatex found. Generating Doxygen and PDF docs..."; \
		doxygen Doxyfile; \
		cd ../docs/latex && make pdf; \
	else \
		echo "[docs] pdflatex not found. Skipping documentation generation."; \
	fi

# Remove files generated by unit tests during execution
cleanup_test_artifacts:
	rm -f $(UNIT_TESTS_DIR)/test_array_module.dat \
	      $(UNIT_TESTS_DIR)/test_pdb_temp.pdb \
	      $(UNIT_TESTS_DIR)/test_assoc.txt \
	      $(UNIT_TESTS_DIR)/test_count.txt \
	      $(UNIT_TESTS_DIR)/test_empty.txt \
	      $(UNIT_TESTS_DIR)/test_read_input_assoc.txt \
	      $(UNIT_TESTS_DIR)/test_exist_check.txt \
	      $(UNIT_TESTS_DIR)/test_header_only.txt \
	      $(UNIT_TESTS_DIR)/test_seq_1.txt \
	      $(UNIT_TESTS_DIR)/test_seq_2.txt \
	      $(UNIT_TESTS_DIR)/test_cluster_clusters.txt \
	      $(UNIT_TESTS_DIR)/test_mean_linkage_clust_info.txt \
	      $(UNIT_TESTS_DIR)/test_mean_linkage_clusters.txt \
	      $(UNIT_TESTS_DIR)/test_mean_linkage_cluster*_complexes \
	      test_array_module.dat \
	      test_pdb_temp.pdb \
	      test_assoc.txt \
	      test_count.txt \
	      test_empty.txt \
	      test_read_input_assoc.txt \
	      test_exist_check.txt \
	      test_header_only.txt \
	      test_seq_1.txt \
	      test_seq_2.txt \
	      test_cluster_clusters.txt \
	      test_mean_linkage_clust_info.txt \
	      test_mean_linkage_clusters.txt \
	      test_mean_linkage_cluster*_complexes

# Clean target: remove all compiled files (including test executables)
clean:
	rm -f $(EXE1) $(EXE2) $(EXE3) $(EXE4) *.o *.mod
	rm -f $(TEST_MATHS) $(TEST_ARRAY) $(TEST_THRESHOLD) $(TEST_PDB) $(TEST_CLUSTERING) $(TEST_ASSOC) $(TEST_READ_INPUT)
	@$(MAKE) cleanup_test_artifacts
	rm -rf ../docs/doxygen

# Clean only Fortran test files
clean_tests:
	rm -f $(TEST_MATHS) $(TEST_ARRAY) $(TEST_THRESHOLD) $(TEST_PDB) $(TEST_CLUSTERING) $(TEST_ASSOC) $(TEST_READ_INPUT)
	@$(MAKE) cleanup_test_artifacts
