"""Created on Fri Apr 26 17:16:40 2019  @ author: dasharath"""

import os
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
import numpy as np
import seaborn as sb

# Specifying the path for the file
file_path = r"C:\Users\Nicholas\Desktop\3-25-2022\dts_file"


def read_vx(filename):
    vx = np.loadtxt(file_path + os.sep + filename, skiprows=1, usecols=1)
    # calculating average and standard deviation for time series
    avg = np.mean(vx)
    sd = np.std(vx)
    vx -= avg
    vx /= sd
    return vx

filenumber = input("filenumber >>  ")
file = 'DTS_VT'+str(filenumber).zfill(3)+'.txt'

vx = read_vx(file)
kde_kws={"color": "red", "lw": 3}
hist_kws={"histtype": "stepfilled", "linewidth": 1,"alpha": 1, "color": "g"}
       
fig, ax = plt.subplots(1,1, figsize=(7,7), sharex=False, facecolor="white", dpi= 170)
ax = sb.distplot(vx, bins=None, hist=True, kde=True, color='C6',fit=None, hist_kws= None,
                         norm_hist=True, kde_kws= kde_kws, label = r'42 $\Omega$-cm')
ax.set_xlabel(r"$\frac{\delta R - \overline{\delta R}}{\sigma_{\delta R}}$")
ax.set_ylabel(r"Estimated PDF")
ax.legend(loc='lower center', bbox_to_anchor=(0.80, 0.85), shadow=False, ncol=1,  
           fontsize='small')
#ax.set_xticks([-2, -1, 0, 1, 2])
#ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
#ax.plot(vx, stats.norm.pdf(vx1, avg, sd))
plt.show()
