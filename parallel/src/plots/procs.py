'''

Author: Manel Serrano Rodríguez

'''

# This script plots efficiency of parallelization as a function of the Cores used.
# Principal contributor: Manel Serrano
import numpy as np
import matplotlib.pyplot as plt
import os

# Set path to the directory where the script is located
os.chdir(os.path.dirname(os.path.abspath(__file__)))

data = np.loadtxt('./nprocs.dat', skiprows=1)

processors = data[:, 0]                 # Cores
time = data[:, 1]                       # Time taken by the simulation in seconds
speedup = data[:,2]                     # Speedup
efficiency = speedup / processors       # Efficiency

# Plot 1: Time vs Cores
plt.figure()
plt.plot(processors, time, 'o--', color='#395261', linewidth=1)
plt.xlabel('# of cores')
plt.ylabel('Time (s)')
plt.title('Time')
plt.grid()
plt.savefig('time.png', dpi=300)
#plt.show()

# Plot 2: Speedup vs Cores
plt.figure()
plt.plot(processors, processors, '--', label='Ideal', color='black', linewidth=1)
plt.plot(processors, speedup, 'o-', label='Measured', color='#543B8B')
for i in range(len(processors)):
    plt.text(processors[i], speedup[i]+0.5, str(round(speedup[i], 2)), fontsize=8, ha='center')
plt.xlabel('# of cores')
plt.ylabel('Speedup')
plt.title('Speedup')
plt.legend()
plt.grid()
plt.savefig('speedup.png', dpi=300)
#plt.show()

# Plot 3: Efficiency vs Cores
plt.figure()
plt.plot(processors, efficiency, 'o--', color='#39613A', linewidth=1)
plt.xlabel('# of cores')
plt.ylabel('Efficiency')
plt.title('Efficiency')
plt.grid()
plt.savefig('efficiency.png', dpi=300)
#plt.show()
