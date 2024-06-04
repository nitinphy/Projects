"""Created on Tue Nov 26 10:28:22 2019  @ author: dadhikar"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
from scipy import signal
from scipy.stats.kde import gaussian_kde
from scipy.stats import norm


# Specifying the path for the file
file_path = " "

# reading data from the file (filename) in file_path
def read_file(filename, a, b):
    x0, x1 = np.loadtxt(file_path + os.sep+ filename, skiprows=1,
                        usecols=[a,b], unpack=True)
    return x0, x1

file = input("Enter the name of the file with format:  ") # enter the file name

# t is the time 
# dvx in phase voltage fluctuation  data
# dvy out of phase voltage fluctuation data
t, dvx = read_file(file, 0, 1)
#plt.plot(t, dvx, label = 'Before detrending')
#plt.legend()
#plt.show()
# detrending the voltage fluctuation time series
dvx = signal.detrend(dvx, type='linear')
plt.plot(t, dvx/0.015, label = 'After detrending')
plt.legend()
plt.show()

#sys.exit()

# Oscillating voltage 
vosc = float(input("Enter the oscillating voltage value:  "))

# Normalized voltage fluctuation
dvx_n = dvx/vosc


dvx_n *= 10**(5)

# creating figure and axis instance
fig, ax = plt.subplots(nrows=1, ncols=1, sharex=False, sharey=False, figsize = (5,5), dpi = 200)

# non-parametric pdf estimation
kde_density = gaussian_kde(dvx_n, bw_method=None, weights=None)
bin_range  = np.linspace(dvx_n.min()-2, dvx_n.max()+2, 500)
#bin_range  = np.linspace(-1.5, 2, 500)
#bin_range  = np.linspace(-5, 5, 100)
nparam_density = kde_density(bin_range)
ax.plot(bin_range, nparam_density, 'r-', label='KDE')

# parametric fit: assume normal distribution
parameter = norm.fit(dvx_n)    # parameters : mean and standard deviation
param_density = norm.pdf(bin_range, loc=parameter[0], scale=parameter[1])
ax.plot(bin_range, param_density, 'k--', label='Gaussian fit')

# histogram for the time series data
#bins = np.linspace(dvx_n.min(), dvx_n.max(), 100)
#bin_centers = 0.5*(bins[1:] + bins[:-1])
ax.hist(dvx_n, bins=bin_range, density=True, weights=None, cumulative=False, bottom=None,
         histtype='stepfilled', align='mid', orientation='vertical', rwidth=None, log=False, 
        stacked=False, normed=None, data=None,  color='green', label = 'Histogram', alpha = 0.5 )

ax.set_ylabel('Probability density function')
ax.set_xlabel(r'$\frac{\delta R}{R}[* 10^{-5}]$')
ax.legend()
plt.show()


# -------------------------------------------------------------------------------------------------- 
# writing estimated pdf data into a  file

File_Name = input('Enter the filename:  ')
FileToWrite = open ("/home/dadhikar/Desktop/__CuIr2S4_V2__/Thermal driven study/pdf/" + os.sep + 
                    File_Name, 'w')
FileToWrite.write('bins (*10**(-5))' +'\t'+ 'npara_density' +'\t'+ 'para_density' + '\n')

for i in range(len(bin_range)):
    FileToWrite.write(str(bin_range[i]) +'\t'+ str(nparam_density[i]) +'\t'+ str(param_density[i]) 
                                                                                              +'\n')

FileToWrite.close()

# --------------------------------------------------------------------------------------------------















