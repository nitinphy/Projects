"""Created on Tue Nov 13 14:30:27 2018 @author: Dasharath Adhikari"""

import os
import sys
import pandas as pd 
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

# file location
file_path = " "


def read_file(filename, a, b):
    x0, x1 = np.loadtxt(file_path + os.sep + filename, skiprows=1,
                        usecols=[a, b], unpack=True)
    return x0, x1


t_165K, v_165K = read_file('time_series_165K_2.txt', 0, 1)
t_200K, v_200K = read_file('time_series_200K_1.txt', 0, 1)
# Plotting time series
fig, ax = plt.subplots(1, 1, figsize=(10, 10), dpi=250)
ax.plot(t_165K, v_165K, "C0", label=r"T = 165 K")
# ax.plot(t_200K, v_200K-2, "C1", label = r"T = 200 K")
ax.set_xlabel(r't (s)',  fontsize=18)
ax.set_ylabel(r'$\frac{\delta R}{R}[10^{-3}]$ ', fontsize=18)
# ax.set_xlim(0.001,1)
ax.set_ylim(-3, 3)
ax.set_title("Time Series", fontsize=18)
ax.legend()
ax.set_facecolor('white')
plt.show()
