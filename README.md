# Project-III: Simulation of a Van der Waals gas using Molecular Dynamics.

## Brief description of the project

This project involves the development of a Molecular Dynamics code to simulate a Van der Waals gas, performing calculations in both serial and parallel modes.

Additionally, this project requires collaborative work to ensure the proper functioning of the code.


## Team Members and Responsibilities

Each team member has assigned tasks, which are indicated below:

4. **Anna Monclús (@anna-mr98)**:  Initialize the configuration and define boundary conditions. Also coordinates the GitHub repository.
1. **Aina Gaya (@ainagaya)**: Integration of Newton's equations.
2. **Albert Plazas (@Alplalo)**: Compute the forces for a Van der Waals interaction.
3. **Manel Serrano (@gluoon8)**:  Post-processing of data, statistics, and documentation of readme and makefile.



## Prerequisites
To execute the program, there are some pre-requisites:
- Make: to execute the program (https://www.gnu.org/software/make/#download)
- Gfortran: to run the MD simulation. (https://gcc.gnu.org/wiki/GFortran)
- Python 3.x : to generate the plots after simulation
  - Numpy (https://numpy.org/install/)
  - Matplotlib (https://matplotlib.org/stable/users/installing/index.html)
  - Scipy (https://scipy.org/)
  - Numba (https://numba.pydata.org/)


## How to

1. Clone repository to your local host
2. Choose *serial* or *parallel* folder with `cd serial` or `cd parallel` 
3. Use `make` or `make help` to see available commands.
4. Before starting a simulation, change your parameters in *namMD.nml* file  
5. To carry out the simulation, have a look to the ***quick quide***. 
6. Data is generated in *.dat* files. If you want to generate figures, use `make plot`


### Quick serial guide

To carry out a **serial** simulation after choosing parameters in *namMD.nml* file you can use:
```
make run (only works in serial mode)
make plot
```
And results will appear in your serial directory!


### Parallel run
To perform the simulation of the LJ gas, go to parallel dir with `cd parallel` and modify these 
ones for the ones that you desire. 
- Parameters in namMD.nml:
  - mass_real []
  - rho_real []
  - epsilon_real []
  - sigma_real []
  - Temp_real []
  - tfin_real []

There are two options to run the program in parallel mode: 
- **Run in local with OpenMPI**:
  - To run the job in a local device, OpenMPI is required. You can set the number of processors by changing the value of *N=Nprocessors* in Makefile file, inside the parallel directory. 
   `make run`

- **Run in a cluster**: 
  - To run the parallel code in iqtc07 of cerqt2, it is possible to set the number of processors (workers) by changing
    the variable NP=Nprocessors in openmpi.sub submit file.
  - Then, do `make qsub` and job will be sent to the queue manager. In parallel/ directory it is possible to identify
    errors in *hello.err* file that appear once the job is running. 
  - When the simulation is done, the output files can be found in parallel directory.
  - Do `make plot` to carry out statistics and plotting of relevant magnitudes along the simulation. 


In summary,

To carry out a **parallel** simulation in cerqt2 cluster after choosing parameters in *namMD.nml* file you can use:
(Make sure you choose the processors in *openmpi.sub* file)
```
make qsub
make plot
```
And results will appear in your parallel/results/ directory!

> [!TIP]
> There are some ways to clean generated files, have a look at `make clean`, `make cleandata` and `make cleanplot`.


## Help 
                          

- Commands:                                                       
  **Serial mode only**
  - `make run`: Compiles needed files and also runs the serial program.     
  
  **Parallel mode only**
  - `make multirun`: runs the job with openMPI in local device
  - `make qsub`: runs the job in iqtc07 queue in cerqt2

  **Works in both modes** 
  - `make plot`: Plots the output data:                              
     - Epot, Ekin, Etot vs time                                   
     - Momentum vs time                                           
     - T vs time                                                  
     - Pressure vs time
     - Radial Distribution Function
     - Averages and stdevs in Averages.dat
                                                
  - `make all`: Compiles the program and creates executable MD.exe   
 
  - `make clean`: Removes the modules, objects and executable        

  - `make cleandata`: Removes data files                             
 
  - `make cleanplot`: Removes plot files                             
 



## Contributors
|  Anna Monclús  |  Aina Gaya  |  Albert Plazas   |  Manel Serrano  |
| -------------- | ----------------- | ------------------ | ------------- |
| ![anna-mr98](./docs/anna-mr98.png "anna-mr98") | ![ainagaya](./docs/ainagaya.png "ainagaya") | ![Alplalo](./docs/Alplalo.png "Alplalo") | ![gluoon8](./docs/gluoon8.png "gluoon8") |
| [anna-mr98](https://github.com/anna-mr98)                                 | [ainagaya](https://github.com/ainagaya)| [Alplalo](https://github.com/Alplalo)                                  | [gluoon8](https://github.com/gluoon8)                                  |



Work developed in the Advanced Informatic Tools subject from [Master of Atomistic and Multiscale Computational Modelling in Physics, Chemistry and Biochemistry](http://www.ub.edu/computational_modelling/).

<table align="center">
  <tr>
    <td><img src="./docs/UB.png" alt="Logo UB"></td>
    <td><img src="./docs/UPC.png" alt="Logo UPC"></td>
  </tr>
  <tr>
    <td>Universitat de Barcelona</td>
    <td>Universitat Politècnica de Catalunya</td>
  </tr>
</table>
