# Define the Fortran compilers
MPIFORTRAN = mpif90

# Define the flags for each compiler
MPIFORTRANFLAGS = -O3 -fbounds-check -ffree-line-length-none -cpp

# Define the source files
SRCS = utils.f90 main.f90

# Define the object files
OBJS = $(SRCS:.f90=.o)

FC = $(MPIFORTRAN)
FFLAGS = $(MPIFORTRANFLAGS)
DEFINES = -D_MPI
EXE = gensdf_mpi

# Default target
all: $(EXE)

# Rule to build the executable
$(EXE): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $(EXE)

# Rule to compile Fortran files
%.o: %.f90
	$(FC) $(FFLAGS) $(DEFINES) -c $< -o $@

# Only remove module/object files
clean:
	rm -f $(OBJS) *.mod

# All Clean target to remove object files and the executable
allclean:
	rm -f $(OBJS) *.mod gensdf_mpi

# Phony targets
.PHONY: all clean
