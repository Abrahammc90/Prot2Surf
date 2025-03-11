# Compiler and flags
FC = gfortran
FLAGS = -O0 -Wall -fno-fast-math

# Source files
MOD_PDB_SRC = mod_pdb.f90
MOD_ASSOC_SRC = mod_assoc.f90
READ_INPUT_SRC = read_input.f90
MATHS_SRC = maths.f90
MOD_MATRIX_SRC = mod_matrix.f90
MOD_CLUST_ALGORITHM_SRC = mod_clust_algorithm.f90
MAKE_MATRIX_SRC = make_matrix.f90
ANALYZE_RESIDUES_SRC = analyze_residues.f90
CLUST_SRC = clust.f90

# Object files
MOD_PDB_OBJ = mod_pdb.o
MOD_ASSOC_OBJ = mod_assoc.o
READ_INPUT_OBJ = read_input.o
MATHS_OBJ = maths.o
MOD_MATRIX_OBJ = mod_matrix.o
MOD_CLUST_ALGORITHM_OBJ = mod_clust_algorithm.o
MAKE_MATRIX_OBJ = make_matrix.o
ANALYZE_RESIDUES_OBJ = analyze_residues.o
CLUST_OBJ = clust.o

# Executable
EXE1 = make_matrix
EXE2 = clust
EXE3 = analyze_residues

# Default target: build the executable
all: $(EXE1) $(EXE2) $(EXE3)

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



# Rule to compile mod_matrix (depends on maths)
$(MOD_MATRIX_OBJ): $(MOD_MATRIX_SRC) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(MOD_MATRIX_SRC)



# Rule to compile make_matrix (depends on everything else)
$(MAKE_MATRIX_OBJ): $(MAKE_MATRIX_SRC) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_MATRIX_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(MAKE_MATRIX_SRC)

# Rule to compile clust (depends on mod_matrix and mod_clust_algorithm)
$(CLUST_OBJ): $(CLUST_SRC) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_MATRIX_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -c $(CLUST_SRC)

# Rule to compile analyze_residues (depends on read_input and maths)
$(ANALYZE_RESIDUES_OBJ): $(ANALYZE_RESIDUES_SRC) $(MATHS_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -c $(ANALYZE_RESIDUES_SRC)


# Rule to link the make_matrix executable (depends on all object files)
$(EXE1): $(MAKE_MATRIX_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_MATRIX_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -o $(EXE1) $(MAKE_MATRIX_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_MATRIX_OBJ) $(READ_INPUT_OBJ) $(MATHS_OBJ)

# Rule to link the clust executable (depends on mod_matrix and mod_clust_algorithm)
$(EXE2): $(CLUST_OBJ) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ) $(MOD_MATRIX_OBJ) $(MATHS_OBJ)
	$(FC) $(FLAGS) -o $(EXE2) $(CLUST_OBJ) $(MOD_CLUST_ALGORITHM_OBJ) $(READ_INPUT_OBJ) $(MOD_MATRIX_OBJ) $(MATHS_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)

# Rule to link the analyze_residues executable (depends on mod_matrix and mod_clust_algorithm)
$(EXE3): $(ANALYZE_RESIDUES_OBJ) $(MATHS_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)
	$(FC) $(FLAGS) -o $(EXE3) $(ANALYZE_RESIDUES_OBJ) $(MATHS_OBJ) $(READ_INPUT_OBJ) $(MOD_PDB_OBJ) $(MOD_ASSOC_OBJ)


# Clean target: remove all compiled files
clean:
	rm -f $(EXE1) $(EXE2) $(EXE3) *.o *.mod
