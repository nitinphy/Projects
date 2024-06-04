"""
Created on Wed Feb 14 16:18:47 2018
@author: Dasharath Adhikari
"""
import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


file_path = "/Users/dadhikar/Box Sync/CuIr2S4/voltage_noise/D1/time_series/230K/psd/"
            
file_list = os.listdir(file_path)            
# print(file_list)           
for file in len(file_list):
    if file.endswith('.txt'):
        print(file)
           
"""         
#print(current_matrix.shape) 
#print(current_matrix[3]) 
Tem = np.linspace(200,243, 44, endpoint=True)
Temperature = np.repeat(Tem,2)
Voltage = np.linspace(0,5,2501, endpoint=True)
X, Y = np.meshgrid(Voltage, Temperature)
#print(Y.shape)
#levels = MaxNLocator(nbins=15).tick_values(current_matrix.min(), current_matrix.max())
cmap = plt.get_cmap('inferno')
#norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
fig, ax = plt.subplots(figsize=(8,8))  
p = ax.pcolor(X, Y, current_matrix, cmap=cmap)
ax.axis('tight')
ax.set_ylabel(r"$Temperature(K)$", fontsize=15)
ax.set_xlabel(r"$Voltage(V)$", fontsize=15)
cb = fig.colorbar(p, ax=ax)
cb.set_label(r'$Current(A)$', fontsize=15)
"""
           
            
            
            
        
        
 
           
                    
                    
                    
                    
                    
                    
                        
            
        
        