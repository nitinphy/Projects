"""Created on Wed Mar 27 19:58:41 2019
 @author: dasharath"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sc
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True


def alpha_calculate(f, S, f_range):
    """
    Parameters
    ----------
    f : 1D- array
        freqeuncy (in Hz)
    S : 1D- array
        PSD extimated
    f_range : [f_min, f_max]
        listfor lower and upper 
        frequencies
    Returns
    -------
    slope : number
        slope of the plot
    std_err : number
        error in slope estimation
    
    """
    f_list = []
    S_list = []
    for i in range(len(f)):
        if f[i] >= f_range[0] and f[i] <= f_range[1]:
            f_list.append(f[i])
            S_list.append(S[i])
        else:
            continue
    f = np.asarray(f_list)
    S = np.asarray(S_list)  
    f = np.log10(f)
    S = np.log10(S)
    slope, intercept, r_value, p_value, std_err = sc.stats.linregress(f, S)
    S_fit = slope* f + intercept
    # f_new = np.linspace(0.001, 1.0, 256)
    # f_new = np.log10(f_new)
    # S_new = slope*f_new + intercept
    slope = round(slope, 3); std_err = round(std_err, 3)
    print("."*30)
    print(r'Fit with Slope: {}+-{}'.format(slope, std_err))
    print('Intercept: ', intercept)
    print('R-Value: ',r_value)
    print('R**2-Value: ',r_value**2)
    print("."*30)
    fig, ax = plt.subplots(1,1, figsize=(7, 7), dpi=150)
    ax.loglog(10**f, 10**S, "bD", label = "PSD")
    ax.loglog(10**f, 10**S_fit, color='k', lw=2.0, ls='--',
            label = r'Fit with Slope: {}+-{}'.format(slope, std_err))
    # ax.loglog(10**f_new, 10**S_new, "g", label = "PSD_new")
    ax.set_xlabel(r"f [Hz]",  fontsize=15)
    ax.set_ylabel(r"PSD [A$^{2}$/Hz]", fontsize=15)
    # ax.set_xlim(f.min(), f.max())
    # ax.set_ylim(10e-14, 10e-1)
    # ax.set_title("")
    ax.legend(loc='best')
    plt.show()
    return (slope, intercept, std_err)



    
   
