# Define the Fortran compilers
NVFORTRAN = nvfortran
GFORTRAN = gfortran

# Define the flags for each compiler
GFORTRANFLAGS = -O3 -fbounds-check -ffree-line-length-none -cpp
NVFORTRANFLAGS = -O3 -cuda -fast -Mpreprocess

# Define the source files
SRCS = utils.f90 main.f90

# Define the object files
OBJS = $(SRCS:.f90=.o)

# Check if _ISCUDA is set and define the compiler and flags accordingly
ifeq ($(ISCUDA),True)
  FC = $(NVFORTRAN)
  FFLAGS = $(NVFORTRANFLAGS)
  DEFINES = -D_ISCUDA
  EXE = gensdf_cuda
else
  FC = $(GFORTRAN)
  FFLAGS = $(GFORTRANFLAGS)
  DEFINES =
  EXE = gensdf_cpu
endif

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
	rm -f $(OBJS) *.mod gensdf_cpu gensdf_cuda

# Phony targets
.PHONY: all clean
