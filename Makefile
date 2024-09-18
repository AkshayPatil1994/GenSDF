FC = mpif90								# Name of the MPI compiler
EXECUTABLE = gensdf						# Name of the executable
SOURCES = $(wildcard *.f90)				# Sources used to compile the code
OBJECTS = $(SOURCES:.f90=.o)			# Object files generated

# Define how the executable is constructed using the object files
$(EXECUTABLE): $(OBJECTS)
	$(FC) -o $(EXECUTABLE) $(OBJECTS) -J mod/

%.o: %.f90
	$(FC) -ffree-line-length-none -c $< -o $@

# Clean method: Remove object files, module files, and the executable
clean:
	rm -f $(OBJECTS) mod/*.mod
allclean:
	rm -f $(OBJECTS) $(EXECUTABLE) mod/*.mod