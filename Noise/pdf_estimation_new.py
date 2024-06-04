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
class file_dir:
    t_series_dir = os.getcwd() +os.sep+ 'time_series_data/heating_cycle'
       


# reading data from the file (filename) in file_path
def read_file(filenumber, a, b):
    x0, x1 = np.loadtxt(file_dir.t_series_dir + os.sep+ 't_series_up'+str(filenumber)+'.txt',
                        skiprows=1, usecols=[a,b], unpack=True)
    return x0, x1


# enter the file name
filenumber = input("Enter file number in 001 format:  ") 

# t is the simulation time
# dR the detrended normalized resistance fluctuation
t, dR = read_file(filenumber, 0, 1)


plot = False

if plot:
    plt.plot(t, dR, label = 'Normalized resistance fluctuation')
    plt.legend()
    plt.show()



# kernel-density estimate using Gaussian kernels.
kde_density = gaussian_kde(dR, bw_method=None, weights=None)
# defining the bin size in the range of data
bin_range  = np.linspace(dR.min(), dR.max(), 32)
# kde estimation in the bin_range
nparam_density = kde_density(bin_range)

# parametric fit: assume normal distribution
parameter = norm.fit(dR)    # parameters : mean and standard deviation
param_density = norm.pdf(bin_range, loc=parameter[0], scale=parameter[1])

# creating figure and axis instance
fig, ax = plt.subplots(1,1, figsize = (5,5), dpi = 200)
ax.plot(bin_range, nparam_density, 'r-', label='Density estimation')
ax.plot(bin_range, param_density, 'k--', label='Gaussian fit')
# histogram for the time series data
ax.hist(dR, bins=bin_range, density=True, cumulative=False, label = 'Histogram', color='green', alpha = 1.0)
ax.set_ylabel('Probability density function (pdf)')
ax.set_xlabel(r'$\frac{\delta R}{R}$')
ax.legend(loc='best', frameon=True, framealpha=0.5, prop={'size': 6})
#plt.legend()
plt.show()


# writing estimated pdf data into a  file
pdf_filename = 'pdf'+str(filenumber)+'_h.txt'
FileToWrite = open (file_dir.t_series_dir +os.sep+ pdf_filename, 'w')

bins_ = 'bins_'+str(filenumber)+'_h'
nonpar_den_ = 'nonpar_den_'+str(filenumber)+'_h'
par_den_ = 'par_den_'+str(filenumber)+'_h'
FileToWrite.write(bins_ +'\t'+ nonpar_den_ +'\t'+ par_den_ +'\n')

for i in range(len(bin_range)):
    FileToWrite.write(str(bin_range[i]) +'\t'+ str(nparam_density[i])
                      +'\t'+ str(param_density[i]) +'\n')

FileToWrite.close()














