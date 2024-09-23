# Define the Fortran compilers
NVFORTRAN = nvfortran
GFORTRAN = gfortran

# Define the flags for each compiler
GFORTRANFLAGS = -O3 -ffree-line-length-none -cpp
NVFORTRANFLAGS = -O3 -cuda -fast -Mpreprocess

# Define the source files
SRCS = utils.f90 main.f90

# Define the object files
OBJS = $(SRCS:.f90=.o)

# Output executable name
EXE = gensdf

# Check if _ISCUDA is set and define the compiler and flags accordingly
ifeq ($(ISCUDA),True)
  FC = $(NVFORTRAN)
  FFLAGS = $(NVFORTRANFLAGS)
  DEFINES = -D_ISCUDA
else
  FC = $(GFORTRAN)
  FFLAGS = $(GFORTRANFLAGS)
  DEFINES =
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
	rm -f $(OBJS)

# All Clean target to remove object files and the executable
allclean:
	rm -f $(OBJS) $(EXE)

# Phony targets
.PHONY: all clean