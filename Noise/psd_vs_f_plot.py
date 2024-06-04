"""Author @ Dasharath Adhikari"""
import os
import sys
from scipy import signal
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from psd_slope_alpha import alpha_calculate
fontP = FontProperties()
fontP.set_size('xx-small')
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True


# specify the location for the estimated PSD file
path_psd = r"C:\Users\Nicholas\Desktop\Sample_39\psd_file"

def read_psd_file(filenumber):
    file_format = 'psd_'+str(filenumber).zfill(3)+'.txt'
    f, psd = np.loadtxt(path_psd + os.sep + file_format, skiprows=1,
                        usecols=[0, 1], unpack=True)
    return f, psd

f1, s1 = read_psd_file(filenumber= input('filenumber: '))
f2, s2 = read_psd_file(filenumber= input('filenumber: '))
f3, s3 = read_psd_file(filenumber= input('filenumber: '))
f4, s4 = read_psd_file(filenumber= input('filenumber: '))
f5, s5 = read_psd_file(filenumber= input('filenumber: '))
f6, s6 = read_psd_file(filenumber= input('filenumber: '))
# Plotting total spectral density and alpha
fig, ax1 = plt.subplots(1, 1, figsize=(8,8), facecolor="white", dpi=150)

ax1.loglog(f1, s1, label='100')
ax1.loglog(f2, s2, label='101')
ax1.loglog(f3, s3, label='102')
ax1.loglog(f4, s4, label='103')
ax1.loglog(f5, s5, label='104')
ax1.loglog(f6, s6, label='105')
ax1.set_xlabel(r'f (Hz)',  fontsize=12)
ax1.set_ylabel(r'S$_{I}$[A$^{2}$/Hz]', fontsize=12)
plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', prop=fontP)
plt.show()

# alpha_calculate gives frequency exponent, slope 
f_range = [0.02, 3]
# f_range = [0.01, 8.0]
slope, intercept, std_err = alpha_calculate(f=f1, S=s1, f_range=f_range)