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
momentum = np.loadtxt('../../momentum.dat', skiprows=1, dtype=float)

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

pos_fin = np.loadtxt('../../pos_out.dat', skiprows=1, dtype=float)


# Compute the radial distribution function
def rdf(pos, L, N, dr):
    # Number of particles
    n = len(pos)
    # Number of bins
    nbins = int(L/(2*dr))
    # Initialize the histogram
    hist = np.zeros(nbins)
    # Loop over all pairs of particles
    for i in range(n):
        for j in range(i+1, n):
            # Compute the distance between particles i and j
            rij = pos[i] - pos[j]
            # Apply periodic boundary conditions
            rij = rij - L * np.rint(rij/L)
            r = np.linalg.norm(rij)
            # Increment the histogram
            if r < L/2:
                k = int(r/dr)
                hist[k-1] += 2
    # Normalize the histogram
    rho = n / (L**3)
    for i in range(nbins):
        r = (i + 0.5) * dr
        hist[i] /= 4 * np.pi * r**2 * dr * rho * N
    return hist

#
#       PLOT RDF 
#
L = 10.77
N = 125
dr = 0.1
rdf = rdf(pos_fin, L, N, dr)
r = np.linspace(0.5*dr, L/2-0.5*dr, len(rdf))
plt.figure()
plt.plot(r, rdf, label='Radial Distribution Function', color='mediumseagreen')
plt.xlabel('r')
plt.ylabel('g(r)')
plt.legend()
plt.title('Radial Distribution Function')
plt.savefig('Radial_Distribution_function.png')
#plt.show()


#
#       PLOT ENERGIES vs time
#

plt.figure(figsize=(10, 5))
plt.plot(energy[:, 0], energy[:, 1], label='Potential Energy', color='#C75146')
plt.plot(energy[:, 0], energy[:, 2], label='Kinetic Energy', color='#AD2E24')
plt.plot(energy[:, 0], energy[:, 3], label='Total Energy', color='#EA8C55')
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
plt.plot(momentum[:, 0], momentum[:,1], label='Momentum', color='mediumaquamarine', linewidth=2)
# plt.ylim(np.mean(momentum)-1, np.mean(momentum)+1)
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
