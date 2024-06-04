"""
Created on Tue Nov 13 14:30:27 2018 
@author: Dasharath Adhikari
"""

import os
import sys 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm
from matplotlib import cm

v_dc = np.array([0.0, 0.0, 0.2, 0.2, 0.4, 0.4, 0.6, 0.6, 0.8, 0.8, 1.0, 1.0,
                 1.2, 1.2, 1.4, 1.4, 1.6, 1.6, 1.8, 1.8, 2.0, 2.0, 2.4, 2.8,
                 3.0, 3.0, 3.2, 3.2, 3.4, 3.4, 3.6, 3.6, 3.8, 3.8, 4.0, 4.0,
                 4.5, 4.5, 5.0, 5.0])

f_new = np.linspace(0.001, 1.0, 256)

ff, vv = np.meshgrid(f_new, v_dc)

file_path = "/Users/dadhikar/Box Sync/CuIr2S4/voltage_noise/D1/time_series/230K/psd/"
            
file_list = os.listdir(file_path)

psd_matrix = np.zeros((len(v_dc), len(f_new)))            
# print(file_list)
create_psd_matrix = True
if create_psd_matrix:
    for filenumber in range(len(file_list)):
        filename = "PSD_VT"+str(filenumber+1) + ".txt"
        print(filename)
        f, psd = np.loadtxt(file_path + os.sep + filename, usecols = [0,1],
                           skiprows=1, unpack=True)
        psd_matrix[filenumber] = psd*f**1.8*10**8
                

cmap=cm.PiYG
       
fig, ax = plt.subplots(figsize=(8,8)) 
p = ax.pcolor(ff, vv, psd_matrix, 
            norm=LogNorm(vmin=psd_matrix.min(), vmax=psd_matrix.max()),
              cmap=cmap) 
ax.axis('tight')
ax.set_xscale('log')
# ax.set_zscale('log')
ax.set_xlabel(r"f [Hz]", fontsize=14)
ax.set_ylabel(r'V$_{dc}$[V]', fontsize=14)
cb = fig.colorbar(p, ax=ax)
cb.set_label(r'$\frac{f*S_R}{R^{2}} [no. unit]$', fontsize=14)


            