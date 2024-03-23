##
# Fortran compiler
FC=gfortran

# Code location
VPATH=src/


# ~ LINKING ~
## make: Compiles the program.
all: MD.x
MD.x:  pbc.o integrators.o forces.o initialize.o  main.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) $^ -o $@

# ~ COMPILING ~

pbc.o: pbc/pbc.f90
	$(FC) -c $^

integrators.o: integrators/integrators.f90
	$(FC) -c $^

forces.o: forces/forces.f90
	$(FC) -c $^

initialize.o: initialization/initialize.f90
	$(FC) -c $^

main.o: main.f90
	$(FC) -c $^

