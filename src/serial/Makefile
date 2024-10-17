# Define the Fortran compilers
GFORTRAN = gfortran

# Define the flags for each compiler
GFORTRANFLAGS = -O3 -fbounds-check -ffree-line-length-none -cpp

# Define the source files
SRCS = utils.f90 main.f90

# Define the object files
OBJS = $(SRCS:.f90=.o)

FC = $(GFORTRAN)
FFLAGS = $(GFORTRANFLAGS)
DEFINES =
EXE = gensdf_serial

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
	rm -f $(OBJS) *.mod gensdf_serial

# Phony targets
.PHONY: all clean
