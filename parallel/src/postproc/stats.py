'''

Author: Manel Serrano Rodriguez

'''
import numpy as np
import os
import matplotlib.pyplot as plt
from numba import njit
import scipy as sp
import time

eps = 162.58        # [K]
sigma = 3.6274      # [Angstrom]
mass = 83.798       # [amu]

t0 = 0              # Set initial cpu time

#-----------------------------------Functions------------------------------------|

@njit
def calculate_stdev(E,n):
    '''
    Computes the standard deviation of the data using the unbiased estimator
    '''
    # Calculate mean
    sum = 0.
    for i in range(n):
        sum += E[i]/n
    mean = sum

    sum = 0.
    for i in range(n):
        sum += (E[i]-mean)**2/n         # Calculate standard deviation

    stdev = np.sqrt(sum/(n-1))          # Unbiased estimator of the standard deviation (n-1 in the denominator)
    return stdev

@njit
def perform_binning(pot, mbin=1):
    '''
    Perform binning of the data and return the average values and standard deviations for each block size
    '''
    midx = 0
    while ((len(pot)/mbin) > 20):       # Limit the number of bins to 25
        mbin = mbin*2                   # Increase the bin size by a factor of 2
        midx = midx + 1                 # Count the number of bins  

    mlist = np.empty(midx)              # List of bin sizes
    E_exp = np.empty(midx)              # List of average values
    E_stdev = np.empty(midx)            # List of standard deviations

    for i in range(midx):               # Compute bin sizes based on powers of 2
        mlist[i] = int(2**i)
        
    
    for i in range(midx):                           # Perform binning for each block size
        size = int(mlist[i])                        # Set the size of the bin
        E_binned = np.zeros(int(len(pot)/size))     # Initialize the array for the binned values    

        E_exp[i] = 0

        for j in range(int(len(pot)/size)):         # Loop over the number of bins
            E_binned[j] = 0

            for k in range(size):                       # Loop over the size of the bin
                idx = j * size + k                      # Calculate the index of the binned value
                E_binned[j] = E_binned[j] + pot[idx]    # Sum the values in the bin 

            
            E_binned[j] = E_binned[j]/size              # Average values for each bin
            E_exp[i] = E_exp[i] + E_binned[j]           # Sum the average values for each block size

        E_stdev[i] = calculate_stdev(E_binned,int(len(pot)/size))       # Calculate the standard deviation for each block size

        E_exp[i] = E_exp[i]/int((len(pot)/size))                # Calculate the average value for each block size

    return mlist, E_exp, E_stdev


def autocorr(x, a, b, tau):
    return a - b * np.exp(np.clip(-x/tau, -np.inf, np.log(np.finfo('d').max)))


def fitting(autocorr, mlist, stdev):
    parms, cov = sp.optimize.curve_fit(autocorr, mlist, stdev)
    return parms


def plot_stdev(mlist, E_stdev, parms, magnitude):
    '''
    Plot size of the bins vs the standard deviation and save the figure with the name of the plotted magnitude
    '''
    
    x = np.arange(mlist[0], mlist[-1])
    plt.figure(dpi=300)  
    plt.scatter(mlist, E_stdev, label='Standard deviation', color='#9E49B8')
    plt.plot(x, autocorr(x, *parms), label='Fitted autocorrelation', linestyle='--', color='black', linewidth=1)
    plt.legend()
    plt.xscale('log')
    plt.xlabel('Bin size')
    plt.ylabel('$\sigma$')
    plt.title(f'{magnitude} ($\sigma$) vs bin size')
    plt.grid()
    plt.savefig(f'{magnitude}_vs_binsize.png', dpi=300)  


#-------------------------------------------------------------------------|
#                               MAIN PROGRAM                              |
#-------------------------------------------------------------------------|

os.chdir(os.path.dirname(os.path.abspath(__file__)))    # Set path to the directory where the script is located

#---------------------------Read data from files---------------------------|

energy = np.loadtxt('../../energy_verlet.dat', skiprows=4, usecols= (1,2,3,4), dtype=float)
temperature = np.loadtxt('../../Temperatures_verlet.dat', skiprows=1, dtype=float)
pressure = np.loadtxt('../../pressure_verlet.dat', skiprows=1, dtype=float)

energy = energy[int(0.6*len(energy)):,:]                # Remove the first 60% of the data
temperature = temperature[int(0.6*len(temperature)):,1]        # Remove the first 60% of the data
pressure = pressure[int(0.6*len(pressure)):,1]              # Remove the first 60% of the data

#------------------------------Convert units--------------------------------|

energy = energy*eps/1000                                    # Convert energy to kJ/mol
temperature = temperature*eps                               # Convert Temperature to K
pressure = pressure*eps/(((sigma*1e-10)**3*6.022e23))/1e6   # Convert pressure to Pa


pot = energy[:,0]
kin = energy[:,1]
tot = energy[:,2]
momentum = energy[:,3]

#-----------------------------Perform binning for every magnitude-------------------------------|

start = time.time()     # Start time

mlist, pot_exp, pot_stdev = perform_binning(pot)    
mlist, kin_exp, kin_stdev = perform_binning(kin)
mlist, tot_exp, tot_stdev = perform_binning(tot)

mlist, mom_exp, mom_stdev = perform_binning(momentum)
mlist, temp_exp, temp_stdev = perform_binning(temperature)
mlist, pres_exp, pres_stdev = perform_binning(pressure)

end = time.time()       # End time


#----------------------------Autocorrelation time function and plots----------------------------|

potparms = fitting(autocorr, mlist, pot_stdev)      # Fit the autocorrelation function for each magnitude
plot_stdev(mlist, pot_stdev, potparms, 'PotE')

kinparms = fitting(autocorr, mlist, kin_stdev)      
plot_stdev(mlist, kin_stdev, kinparms, 'KinE')

totparms = fitting(autocorr, mlist, tot_stdev)
plot_stdev(mlist, tot_stdev, totparms, 'TotE')

momparms = fitting(autocorr, mlist, mom_stdev)
plot_stdev(mlist, mom_stdev, momparms, 'Mom')

tempparms = fitting(autocorr, mlist, temp_stdev)
plot_stdev(mlist, temp_stdev, tempparms, 'Temp')

presparms = fitting(autocorr, mlist, pres_stdev)
plot_stdev(mlist, pres_stdev, presparms, 'Pres')


#---------------------------------Write binning parameters to file-------------------------------|

if os.path.exists('averages.dat'):
    os.remove('averages.dat')

with open('averages.dat', 'w') as f:
    f.write('\# Statistical analysis of the Verlet simulation\n')
    f.write('\n')
    f.write(f'Data points:{len(energy)}\n')
    f.write(f'bin size used: {mlist[-2]}\n')
    f.write('\n')
    f.write('Averages:\n')
    f.write(f'Potential energy: {pot_exp[-1]:.5f} ± {pot_stdev[-2]:.5f} kJ/mol\n')
    f.write(f'Kinetic energy: {kin_exp[-1]:.5f} ± {kin_stdev[-2]:.5f} kJ/mol\n')
    f.write(f'Total energy: {tot_exp[-1]:.5f} ± {tot_stdev[-2]:.5f} kJ/mol \n')
    f.write('\n')
    f.write(f'Momentum: {mom_exp[-1]:.5f} ± {mom_stdev[-2]:.5f} kJ/mol\n')
    f.write(f'Temperature: {temp_exp[-1]:.5f} ± {temp_stdev[-2]:.5f} K\n')
    f.write(f'Pressure: {pres_exp[-1]:.5f} ± {pres_stdev[-2]:.5f} MPa\n')
    
    f.write('\n')
    f.write(f'Time elapsed: {end-start:.5f} s')



print('|--------Statistical results of the Verlet simulation-------|\n')
print()
print(f'Data points:{len(energy)}')
print(f'bin size used: {mlist[-2]}\n')
print()
print('Averages:')
print(f'Potential energy: {pot_exp[-1]:.5f} ± {pot_stdev[-2]:.5f} kJ/mol')
print(f'Kinetic energy: {kin_exp[-1]:.5f} ± {kin_stdev[-2]:.5f} kJ/mol')
print(f'Total energy: {tot_exp[-1]:.5f} ± {tot_stdev[-2]:.5f} kJ/mol')
print()
print(f'Momentum: {mom_exp[-1]:.5f} ± {mom_stdev[-2]:.5f} kJ/mol')
print(f'Temperature: {temp_exp[-1]:.5f} ± {temp_stdev[-2]:.5f} K')
print(f'Pressure: {pres_exp[-1]:.5f} ± {pres_stdev[-2]:.5f} MPa')
print(f'Time elapsed: {end-start:.5f} s')
print('\n')
print('Data written to file: averages.dat')
