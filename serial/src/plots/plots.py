import numpy as np
import os
import matplotlib.pyplot as plt

############################################

#              LJ parameters

############################################

# Parameters for LJ potential of Kr found in:
'''
Rutkai et Al. (2016). How well does the Lennard-Jones potential 
         represent the thermodynamic properties of noble gases?. 

Molecular Physics. 115. 1-18. 10.1080/00268976.2016.1246760. 
'''

eps = 162.58        # [K]
sigma = 3.6274      # [Angstrom]

mass = 83.798       # [amu]


# Set path to the directory where the script is located
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# Load data from files
energy = np.loadtxt('../../energy_verlet.dat', skiprows=4, dtype=float)
temperature = np.loadtxt('../../Temperatures_verlet.dat', skiprows=1, dtype=float)
pressure = np.loadtxt('../../pressure_verlet.dat', skiprows=1, dtype=float)
momentum = energy[:, 4]

# Convert energy to kJ/mol
energy = energy*eps/1000

# Convert Temperature to K
temperature = temperature*eps

# Convert pressure to Pa
pressure = pressure*eps/(((sigma*1e-10)**3*6.022e23))/1e6

# convert time to ps
# WIP!!
#energy[:, 0] = energy[:, 0]*6.022e23*np.sqrt(mass/1000*sigma**2*1e-20/eps)
#temperature[:, 0] = temperature[:, 0]*6.022e23*np.sqrt(mass/1000*sigma**2*1e-20/eps)
#pressure[:, 0] = pressure[:, 0]*6.022e23*np.sqrt(mass/1000*sigma**2*1e-20/eps)



#
#       PLOT ENERGIES vs time
#

plt.figure(figsize=(10, 5))
plt.plot(energy[:, 0], energy[:, 1], label='Potential Energy', color='coral')
plt.plot(energy[:, 0], energy[:, 2], label='Kinetic Energy', color='darkorange')
plt.plot(energy[:, 0], energy[:, 3], label='Total Energy', color='crimson')
plt.xlabel('Timestep')
plt.ylabel('Energy (kJ/mol)')
plt.legend()
plt.title('Energy vs Step')
plt.savefig('Energies.png')
#plt.show()


#
#      PLOT MOMENTUM vs time
#

plt.figure()
plt.plot(energy[:, 0], momentum, label='Momentum', color='mediumaquamarine')
plt.ylim(np.mean(momentum)-1, np.mean(momentum)+1)
plt.xlabel('Timestep')
plt.ylabel('Momentum\'')
plt.legend()
plt.title('Momentum vs Time')
plt.savefig('Momentum.png')
#plt.show()


#
#      PLOT TEMPERATURE vs time
#

plt.figure()
plt.plot(temperature[:, 0], temperature[:, 1], label='Temperature', color='mediumvioletred')
plt.xlabel('Timestep')
plt.ylabel('Temperature (K)')
plt.legend()
plt.title('Temperature vs Time')
plt.savefig('Temperature.png')
#plt.show()


#
#      PLOT PRESSURE vs time
#

plt.figure()
plt.plot(pressure[:, 0], pressure[:, 1], label='Pressure', color='goldenrod')
plt.xlabel('Time (s)')
plt.ylabel('Pressure (MPa)')
plt.legend()
plt.title('Pressure vs Time')
plt.savefig('Pressure.png')
#plt.show()