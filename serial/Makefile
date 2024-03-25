##
# Fortran compiler
FC=gfortran

# Code location
VPATH=src/


# ~ LINKING ~
## make: Compiles the program.
all: clean MD.exe clean2
MD.exe:  pbc.o forces.o integrators.o  initialize.o  main.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) $^ -o $@

# ~ COMPILING ~

pbc.o: pbc/pbc.f90
	$(FC) -c $^

forces.o: forces/forces.f90
	$(FC) -c $^

integrators.o: integrators/integrators.f90
	$(FC) -c $^

initialize.o: initialization/initialize.f90
	$(FC) -c $^

main.o: main.f90
	$(FC) -c $^

# ~ CLEANING ~
## make clean: Removes all object files and the executable.
clean:
	del -f *.o *.mod MD.exe
clean2:
	del -f *.o *.mod

.PHONY:  clean all clean2