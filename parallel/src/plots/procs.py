# This script plots efficiency of parallelization as a function of the Cores used.

import numpy as np
import matplotlib.pyplot as plt
import os

# Set path to the directory where the script is located
os.chdir(os.path.dirname(os.path.abspath(__file__)))

data = np.loadtxt('../../nprocs.dat', skiprows=1)

processors = data[:, 0]                 # Cores
time = data[:, 1]                       # Time taken by the simulation in seconds
speedup = data[:,2]                     # Speedup
efficiency = speedup / processors       # Efficiency

# Plot time taken by the simulation as a function of the Cores
plt.figure(figsize=(12, 4))

# Subplot 1: Time vs Cores
plt.subplot(1, 3, 1)
plt.plot(processors, time, 'o-')
plt.xlabel('# of cores')
plt.ylabel('Time (s)')
plt.title('Time')
plt.grid()

# Subplot 2: Speedup vs Cores
plt.subplot(1, 3, 2)
plt.plot(processors, speedup, 'o-')
plt.xlabel('# of cores')
plt.ylabel('Speedup')
plt.title('Speedup')
plt.grid()

# Subplot 3: Efficiency vs Cores
plt.subplot(1, 3, 3)
plt.plot(processors, efficiency, 'o-')
plt.xlabel('# of cores')
plt.ylabel('Efficiency')
plt.title('Efficiency')
plt.grid()

plt.tight_layout()
plt.savefig('cores.png', dpi=300)
plt.show()