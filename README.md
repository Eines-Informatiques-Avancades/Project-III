# Project-III: Simulation of a Van der Waals gas using Molecular Dynamics.

## Brief description of the project

This project involves the development of a Molecular Dynamics code to simulate a Van der Waals gas, performing calculations in both serial and parallel modes.

Additionally, this project requires collaborative work to ensure the proper functioning of the code.


## Team Members and Responsibilities

Each team member has assigned tasks, which are indicated below:

4. **Anna Monclús (@anna-mr98)**:  Initialize the configuration and define boundary conditions. Also coordinates the GitHub repository.
1. **Aina Gaya (@ainagaya)**: Integration Newton's equations.
2. **Albert Plazas (@Alplalo)**: Compute the forces for a Van der Waals interaction.
3. **Manel Serrano (@gluoon8)**:  Post-processing of data, statistics, and visualization.



## Prerequisites
To execute the program, there are some pre-requisites:
- Make: to execute the program (https://www.gnu.org/software/make/#download)
- Gfortran: to run the MD simulation. (https://gcc.gnu.org/wiki/GFortran)
- Python 3.x : to generate the plots after simulation
  - Numpy (https://numpy.org/install/)
  - Matplotlib (https://matplotlib.org/stable/users/installing/index.html)




## How to

1. Clone repository to your local host
2. Use `make` or `make help` to see the available commands.
3. Before starting a simulation, change your parameters in *namMD.nml* file  
4. To carry out the simulation, use `make run` and the program will be compiled and run. 
5. Data is generated in *.dat* files. If you want to generate figures, use `make plot`

### Quick guide

To carry out a simulation after choosing parameters in *namMD.nml* file you can use:
```
make run
make plot
```
And files will appear in your main directory!

> [!IMPORTANT]
> Current features are only available for serial code, parallel code is WIP!

> [!TIP]
> There are some ways to clean generated files, have a look at `make clean`, `make cleandata` and `make cleanplot`.


## Help 
                          

- Commands:                                                       
  
  - `make run`: Compiles needed files and also runs the program.     
 
  - `make plot`: Plots the output data:                              
     - Epot, Ekin, Etot vs time                                   
     - Momentum vs time                                           
     - T vs time                                                  
     - Pressure vs time                                           
  - `make all`: Compiles the program and creates executable MD.exe   
 
  - `make clean`: Removes the modules, objects and executable        

  - `make cleandata`: Removes data files                             
 
  - `make cleanplot`: Removes plot files                             
 



## Contributors
|  Anna Monclús  |  Aina Gaya  |  Albert Plazas   |  Manel Serrano  |
| -------------- | ----------------- | ------------------ | ------------- |
| ![anna-mr98](https://github.com/Eines-Informatiques-Avancades/Project-III/tree/master/docs/anna-mr98.png "anna-mr98") | ![ainagaya](https://github.com/Eines-Informatiques-Avancades/Project-III/tree/master/docs/ainagaya.png "ainagaya") | ![Alplalo](https://github.com/Eines-Informatiques-Avancades/Project-III/tree/master/docs/Alplalo.png "Alplalo") | ![gluoon8](https://github.com/Eines-Informatiques-Avancades/Project-III/tree/master/docs/gluoon8.png "gluoon8") |
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
