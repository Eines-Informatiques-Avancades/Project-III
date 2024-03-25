import numpy as np
import os
import matplotlib.pyplot as plt

# Set path to the directory where the script is located
os.chdir(os.path.dirname(os.path.abspath(__file__)))



# Load data from files
energy = np.loadtxt('../../energy_verlet.dat', skiprows=4, dtype=float)
temperature = np.loadtxt('../../Temperatures_verlet.dat', skiprows=1, dtype=float)
pressure = np.loadtxt('../../pressure_verlet.dat', skiprows=1, dtype=float)
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
plt.plot(energy[:, 0], energy[:, 1], label='Potential Energy', color='coral')
plt.plot(energy[:, 0], energy[:, 2], label='Kinetic Energy', color='darkorange')
plt.plot(energy[:, 0], energy[:, 3], label='Total Energy', color='crimson')
plt.xlabel('Step')
plt.ylabel('Energy')
plt.legend()
plt.title('Energy vs Step')
plt.savefig('Energies.png')
#plt.show()


#
#      PLOT MOMENTUM vs time
#

plt.figure()
plt.plot(energy[:, 0], energy[:, 4], label='Momentum', color='mediumaquamarine')
plt.ylim(np.mean(energy[:, 4])-1, np.mean(energy[:, 4])+1)
plt.xlabel('Time')
plt.ylabel('Momentum')
plt.legend()
plt.title('Momentum vs Time')
plt.savefig('Momentum.png')
#plt.show()


#
#      PLOT TEMPERATURE vs time
#

plt.figure()
plt.plot(temperature[:, 0], temperature[:, 1], label='Temperature', color='mediumvioletred')
plt.xlabel('Time')
plt.ylabel('Temperature')
plt.legend()
plt.title('Temperature vs Time')
plt.savefig('Temperature.png')
#plt.show()


#
#      PLOT PRESSURE vs time
#

plt.figure()
plt.plot(pressure[:, 0], pressure[:, 1], label='Pressure', color='goldenrod')
plt.xlabel('Time')
plt.ylabel('Pressure')
plt.legend()
plt.title('Pressure vs Time')
plt.savefig('Pressure.png')
#plt.show()
