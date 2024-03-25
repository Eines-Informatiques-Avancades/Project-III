
# Fortran compiler
FC=gfortran

# Code location
VPATH=src/


###################################################################
##                             HELP                              ##
###################################################################
##                                                                #
## - Commands:                                                                #
.PHONY. : help
help:
	@grep "^##" Makefile


# ~ LINKING ~
##  · make: Compiles the program and creates executable MD.exe    #
##                                                                #
all: MD.exe
MD.exe:  pbc.o forces.o integrators.o  initialize.o  main.o
	$(FC) $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) $^ -o $@

# ~ RUNNING THE PROGRAM ~
##  · make run: Runs the program. If the program is not compiled, # 
##            it compiles it first.                               #
##                                                                #

# ~ CLEANING ~
##  · make clean: Removes the object files and the executable.    #
##                                                                #

.PHONY.: cleanall
cleanall:
	rm -f *.o *.mod MD.exe 

.PHONY.: clean
clean:
	rm -f *.dat

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


# ~ PLOTS ~
##  · make plot: Plots the output data:                           #
##                                                                #
##      - Epot, Ekin, Etot vs time                                #
##      - Momentum vs time                                        #
##      - T ,  vs time                                            #
##      - Pressure vs time                                        #
##                                                                #

.PHONY.: plot
plot:
	@echo "Plotting files..."
	@python3 ./src/plots/plots.py
	@echo "Done! Check the plots in the 'plots' folder."


##                                                                #
##                                                                #
###################################################################