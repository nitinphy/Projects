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
# location of csv data from the Keithley
base_directory = r"X:\Personal Data\Nbo2\09202023\206\noise final set\16min\time series"

# Construct full paths to different directories
path_csv = base_directory
path_text = os.path.join(base_directory, "text_file")
path_detrend = os.path.join(base_directory, "detrend_file")
path_DTS = os.path.join(base_directory, "dts_file")
path_DTS_padding = os.path.join(base_directory, "dts_padding_file")
path_psd = os.path.join(base_directory, "psd_file")
plots_folder = os.path.join(base_directory, "plots")
#path_psdn = r"C:\Users\Nicholas\Desktop\Sample_102\psdn_file"
# Create the "plots" folder in the same directory as the CSV files
plots_folder = os.path.join(path_csv, "plots")
if not os.path.exists(plots_folder):
    os.makedirs(plots_folder)
plt.rcParams['axes.linewidth'] = 2
#..........................................................................
filenumber = input('Enter the psd file : ')

# Update the font settings
font = {'family' : 'Times New Roman',
        'weight' : 'normal',
        'size'   : 18}

plt.rc('font', **font)

# Update the line thickness
plt.rc('lines', linewidth=2)

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
    global rsq
    rsq = r_value**2
    # f_new = np.linspace(0.001, 1.0, 256)
    # f_new = np.log10(f_new)
    # S_new = slope*f_new + intercept
    slope = round(slope, 3); std_err = round(std_err, 3)
    #global foooi,fooi,foi,fi
    foooi=(0.0039**(slope))*(10**(intercept))
    fooi=(0.01**(slope))*(10**(intercept))
    foi=(0.1**(slope))*(10**(intercept))
    fi=(1**(slope))*(10**(intercept))
    print("."*30)
    print(r'Fit with Slope: {}+-{}'.format(slope, std_err))
    print('Intercept: ', intercept)
    print('R-Value: ',r_value)
    print('R**2-Value: ',rsq)
    print("."*30)
    print('Noise at 3.9 mHz', foooi)
    print('Noise at 10 mHz', fooi)
    print('Noise at 100 mHz', foi)
    print('Noise at 1 Hz', fi)
    print("."*30)
    fig, ax = plt.subplots(1,1, figsize=(10,8), dpi=150)
    ax.loglog(10**f, 10**S, "bD", label = "PSD")
    ax.loglog(10**f, 10**S_fit, color='k', lw=2.0, ls='--',
            label = r'Fit with Slope: {}+-{}'.format(slope, std_err))
    # ax.loglog(10**f_new, 10**S_new, "g", label = "PSD_new")
    ax.set_xlabel(r"f [Hz]",  fontsize=18)
    ax.set_ylabel(r"PSD [A$^{2}$/Hz]", fontsize=18)
    # ax.set_xlim(f.min(), f.max())
    # ax.set_ylim(10e-14, 10e-1)
    # ax.set_title("")
    ax.legend(loc='best')
    plt.show()
    file_name_4 = 'psd_plot_slope_' + filenumber + '.png'
    plt.savefig(os.path.join(plots_folder, file_name_4), dpi=300, bbox_inches='tight')
    return (slope, intercept, std_err)



    
   
