import os
import sys
import numpy as np
import matplotlib.pyplot as plt

file_dir = r"C:\Users\Nicholas\Desktop\Sample_39\psd_file"

f, s, s_fit  = np.loadtxt(file_dir+ os.sep+ 'PSD_325K.txt', usecols=[0,1, 2],skiprows=1, unpack=True)

fig, ax = plt.subplots(figsize=(5, 5), dpi=150)
ax.loglog(f, s, 'o', label='PSD')
ax.loglog(f, f*s_fit, 'o', label='PSD')
#ax.set_xlabel('ime(min)',fontsize=25)
#ax.set_ylabel('$V_{x}(V)$', fontsize=25)
#ax.legend()
#ax.xaxis.set_ticks_position('both')
#ax.yaxis.set_ticks_position('both')
#ax[1].xaxis.set_ticks_position('both')
#ax[1].yaxis.set_ticks_position('both')
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax[1].plot(time_m, V_y, color='red', lw=2, ls='-', label='$325$')
#ax[1].set_xlabel('Time(min)',fontsize=14)
#ax[1].set_ylabel('$V_{y}(V)$', fontsize=14)
#ax[1].legend()
plt.show()

