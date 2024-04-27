
# Fortran compiler
FC=mpif90

# Number of processors
N=5

# Code location
VPATH=src/


######################################################################
##               MOLECULAR DYNAMICS /    HELP                       ##
######################################################################
##                                                                   #
## - Contributors:                                                   #
##   - @anna-mr98 (Anna Monclús)                                     #
##   - @ainagaya (Aina Gaya)                                         #
##   - @Alplalo (Albert Plazas)                                      #
##   - @gluoon8 (Manel Serrano)                                      #
##                                                                   #
##                                                                   #
##    !IMPORTANT: Remember   to   choose    your                     #
##                parameters in namMD.nml file!!                     #
##                                                                   #
##                                                                   #
## - Commands:                                                       #
##                                                                   #
.PHONY. : help
help:
	@echo " "
	@grep "^##" Makefile

# ~ RUNNING THE PROGRAM ~
##  · make run: Compiles needed files and also runs the program.     #
##                                                                   #
.PHONY.: run
run: MD.exe
	mpirun -n ${N} ./MD.exe

# ~ PLOTS ~
##  · make plot: Plots the output data:                              #
##                                                                   #
##      - Epot, Ekin, Etot vs time                                   #
##      - Momentum vs time                                           #
##      - T vs time                                                  #
##      - Pressure vs time                                           #
##                                                                   #

.PHONY.: plot
plot:
	@echo "Plotting files..."
	@python3 ./src/plots/plots.py
	@cp ./src/plots/*.png ./
	@echo "Done! Check the plots in the main directory."

# ~ LINKING ~
##  · make all: Compiles the program and creates executable MD.exe   #
##                                                                   #
all: MD.exe
MD.exe:  pbc.o forces.o integrators.o  initialize.o  main.o
	$(FC) -O2 $(F_FLAGS) $(COMP_D_FLAGS) $(COMP_R_FLAGS) $^ -o $@


# ~ CLEANING ~
##  · make clean: Removes the modules, objects and executable        #
##                                                                   #
##  · make cleandata: Removes data files                             #
##                                                                   #
##  · make cleanplot: Removes plot files                             #
##                                                                   #
.PHONY.: clean
clean:
	@rm -rf *.o *.mod MD.exe hello.* HELLO
	@echo "Modules, objects and executable removed!" 

.PHONY.: cleandata
cleandata:
	@rm -f *.dat
	@echo "Data files removed!"

.PHONY.: cleanplot
cleanplot:
	@rm -f *.png
	@echo "Plot files removed!"

.PHONY.: qsub
qsub:
	@qsub -q iqtc07.q openmpi.sub
	@echo "Enviado a cola iqtc07"

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

##                                                                   #
##                                                                   #
######################################################################
