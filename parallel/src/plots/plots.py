import numpy as np
import os
import matplotlib.pyplot as plt
from numba import njit
import scipy as sp
import time
import matplotlib.animation as animation

############################################

#              LJ parameters

############################################

# Parameters for LJ potential of Kr found in:
"""
Rutkai et Al. (2016). How well does the Lennard-Jones potential 
         represent the thermodynamic properties of noble gases?. 

Molecular Physics. 115. 1-18. 10.1080/00268976.2016.1246760. 
"""

eps = 162.58        # [K]
sigma = 3.6274      # [Angstrom]
mass = 83.798       # [amu]

os.chdir(os.path.dirname(os.path.abspath(__file__)))    # Set path to the directory where the script is located


# Load data from files
energy = np.loadtxt("../../energy_verlet.dat", skiprows=4, dtype=float)
temperature = np.loadtxt("../../Temperatures_verlet.dat", skiprows=1, dtype=float)
pressure = np.loadtxt("../../pressure_verlet.dat", skiprows=1, dtype=float)
momentum = energy[:, 4]


energy[:,1:] = energy[:,1:] * eps / 1000                                         # Convert energy to kJ/mol
energy[:,0] = energy[:,0]* sigma * (mass/eps)**(0.5)                             # Convert time to ps
temperature[:,1:] = temperature[:,1:] * eps                                      # Convert Temperature to K
temperature[:,0] = temperature[:,0]* sigma * (mass/eps)**(0.5)                   # Convert time to ps
pressure[:,1:] = pressure[:,1:] * eps / ((sigma * 1e-10) ** 3 * 6.022e23) / 1e6  # Convert pressure to MPa
pressure[:,0] = pressure[:,0]* sigma * (mass/eps)**(0.5)                         # Convert time to ps

#-----------------------------------FUNCTIONS-----------------------------------|

@njit                                  # Decorator to compile the function with the JIT compiler
def rdf(pos, L, N, dr):
    
    n = len(pos)                                # Number of particles
    nbins = int(L/(2*dr))                       # Number of bins
    hist = np.zeros(nbins)                      # Initialize the histogram
    
    for i in range(n):                          # Loop over all pairs of particles
        for j in range(i+1, n):
            rij = pos[i] - pos[j]               # Distance between particles i and j
            rij = rij - L * np.rint(rij/L)      # Apply periodic boundary conditions
            r = np.linalg.norm(rij)
            
            if r < L/2:                                 # Increment the histogram
                k = int(r/dr)
                hist[k-1] += 2
    
    rho = N / (L**3)                                    # Density
    for i in range(nbins):
        r = (i + 0.5) * dr                              # Compute the radial distance
        hist[i] /= 4 * np.pi * r**2 * dr * rho * N      # Normalize the histogram
    return hist

def update(frame):          # Define the update function for the animation

    line.set_data(traj[frame][:, 0], traj[frame][:, 1])
    line.set_3d_properties(traj[frame][:, 2])
    return line,



#--------------------------------POSICIONS INICIALS--------------------------------|
pos_ini = np.loadtxt("../../pos_ini.dat", dtype=float)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(pos_ini[:, 0], pos_ini[:, 1], pos_ini[:, 2], color='#8F4B53')
ax.set_title('Initial Particle Positions')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.savefig('Initial_positions.png', dpi=300)
#plt.show()



#-----------------------------------TRAJECTORIES-----------------------------------|

pos_fin = np.loadtxt("../../traj.xyz", skiprows=1, usecols=(1,2,3),dtype=float)

N = 125
n_frames = int(len(pos_fin)/N)              # Number of frames
traj = np.array_split(pos_fin, n_frames)    # Split the array into frames

fig = plt.figure()                          # Create a figure for the animation
ax = fig.add_subplot(111, projection='3d')  # Add a 3D subplot to the figure
line, = ax.plot(traj[0][:, 0], traj[0][:, 1], traj[0][:, 2], 'o', color='#8F4B53', markeredgecolor='#662D34')   # Initialize the plot with the first frame

ani = animation.FuncAnimation(fig, update, frames=n_frames, interval=200)   # Create the animation

ax.set_title('Particle Trajectory')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ani.save('Trajectory.gif', writer='imagemagick', dpi=300, fps=10)   # Create a .gif file
#plt.show()      # Show the animation (disable if working on a cluster to avoid errors)



#----------------------------------RDF plot-------------------------------|
L = 8.55
N = 125
dr = 0.04
rdf_hist = np.zeros(int(L/(2*dr)))                  
count=0
for i in range(N,int(len(pos_fin)),N):
    rdf_hist_sum = rdf(pos_fin[i-N:i], L, N, dr)
    rdf_hist += rdf_hist_sum
    count+=1

rdf_hist /= (int(len(pos_fin)/N)-1)

#print('count=',count)
#print((int(len(pos_fin)/N)-1))

r = np.linspace(0.5*dr, L/2-0.5*dr, len(rdf_hist))
#print(len(rdf_hist))

plt.figure()
plt.axhline(y=1, color='black', linestyle='--', linewidth=0.5)
plt.plot(r, rdf_hist, label='g(r)', color='#5782FB')
plt.xlabel('r(Angstrom)')
plt.ylabel('g(r)')
plt.title('Radial Distribution Function')
plt.savefig('rdf.png', dpi=300)
#plt.show()

#
#       PLOT ENERGIES vs time
#

plt.figure(figsize=(10, 5))
plt.plot(energy[:, 0], energy[:, 1], label="Potential Energy", color="#042940")
plt.plot(energy[:, 0], energy[:, 2], label="Kinetic Energy", color="#005C53")
plt.plot(energy[:, 0], energy[:, 3], label="Total Energy", color="#9FC131")
plt.xlabel("Timestep (ps)")
plt.ylabel("Energy (kJ/mol)")
plt.legend(loc='lower right')
plt.title("Energy vs Step")
plt.savefig("Energies.png", dpi=300)
# plt.show()


#
#      PLOT MOMENTUM vs time
#

plt.figure()
plt.plot(
    energy[:, 0], momentum, label="Momentum", color="mediumaquamarine", linewidth=2
)
#plt.ylim(np.mean(momentum) - 3, np.mean(momentum) + 3)
plt.xlabel("Timestep (ps)")
plt.ylabel("Momentum'")
plt.legend()
plt.title("Momentum vs Time")
plt.savefig("Momentum.png", dpi=300)
# plt.show()


#
#      PLOT TEMPERATURE vs time
#

plt.figure()
plt.plot(
    temperature[:, 0], temperature[:, 1], label="Temperature", color="mediumvioletred"
)
plt.xlabel("Timestep (ps)")
plt.ylabel("Temperature (K)")
plt.legend()
plt.title("Temperature vs Time")
plt.savefig("Temperature.png", dpi=300)
# plt.show()


#
#      PLOT PRESSURE vs time
#

plt.figure()
plt.plot(pressure[:, 0], pressure[:, 1], label="Pressure", color="goldenrod")
plt.ylim(0,5)
plt.xlabel("Timestep (ps)")
plt.ylabel("Pressure (MPa)")
plt.legend()
plt.title("Pressure vs Time")
plt.savefig("Pressure.png", dpi=300)
# plt.show()
