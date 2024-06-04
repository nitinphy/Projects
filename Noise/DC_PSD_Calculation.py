# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 12:56:48 2021

@author: dasharath, edited by Nicholas
"""
import os
from scipy import signal
from psd_slope_alpha import alpha_calculate
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
import numpy as np

# Specifying the path for the file (Must be padded decimated time series file location)
file_path = r"C:\Users\Nicholas\Desktop\1-29-2023\dts_padding_file"
path_psdnfit = r"C:\Users\Nicholas\Desktop\1-29-2023\psdnfit_file"
# specify the location for the estimated PSD file
#..........................................................................
filenumber = input('Enter the txt file number: ')
#..........................................................................
file_to_read = 'dts_padding_it_'+str(filenumber).zfill(3)+'.txt'

# reading time, vx and vy from the file
def read_dts(filepath, filename, a, b):
    x0, x1 = np.loadtxt(filepath + os.sep + filename, skiprows=1,
                            usecols=[a, b], unpack=True)
    return x0, x1

time, fluctuation = read_dts(file_path, file_to_read, 0, 1)

#avgfluctuation = fluctuation.mean()

fs = 512.0  # data sampling frequency
D1 = 4.0  # first stage decimation factor
D2 = 4.0  # second stage decimation factor
D3 = 2.0  # third stage decimation factor
decimation_factor = D1*D2*D3
f_pass_final = (0.5/decimation_factor) * fs

nperseg = 1024*4

def power_spectrum_density(x, spectrum=False):
    """
    Caculate spectral density using Welch periodogram method.
    x : array_like, time series of measurement values
    fs :sampling frequency of the x time series in units of Hz
    window : desired window, defaults to ‘hanning’
    nperseg : length of each segment. Defaults to 256
    noverlap: number of points to overlap between segments
    If None, noverlap = nperseg / 2
    nfft : Length of the FFT used, if a zero padded FFT is desired
    If None, the FFT length is nperseg. Defaults to None.
    detrend: specifies how to detrend each segment. defaults to ‘constant’.
    return_onesided : bool,if True, return a one-sided spectrum for real data
    if False return a two-sided spectrum
    scaling : { ‘density’, ‘spectrum’ }, optional
    """
    if spectrum:
        f, S = signal.welch(x, fs=fs/(D1*D2*D3), window='hann',
                            nperseg=nperseg, noverlap=0.5*nperseg, nfft=None,
                            detrend='constant', return_onesided=True,
                            scaling='spectrum', axis=-1)
    else:
        f, S = signal.welch(x, fs=fs/(D1*D2*D3), window='hann',
                            nperseg=nperseg, noverlap=0.5*nperseg, nfft=None,
                            detrend='constant', return_onesided=True,
                            scaling='density', axis=-1)
    return f, S


# spectral density estimation from vx data
fx, Sx = power_spectrum_density(fluctuation)

# not selecting zero frequency component
fx = fx[1:]
Sx = Sx[1:]
#Sxn = Sx[1:]/(1000**2)
Sxn = (Sx[1:])/(iavg)

def psd_plot(plot=False):
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 7), dpi=150)
        ax.loglog(fx, Sx, c="r", lw=2.0, label="PSD (unnormalized)")
        #ax.loglog(fx, Sxn, c="g", lw=2.0, label="PSD (normalized)")
        ax.set_xlabel(r'f [Hz]', fontsize=15)
        ax.set_ylabel(r'S$_{I}$ [$\frac{I^{2}}{Hz}$]', fontsize=15)
        # ax.set_xlim(f.min(), f.max())
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_title("")
        ax.set_title('Power Spectral Density')
        ax.legend(loc='best')
        plt.show()

psd_plot(plot=False)

# alpha_calculate gives frequency exponent, slope 
f_range = [0.001, 1.0]
# f_range = [0.01, 8.0]
slope, intercept, std_err = alpha_calculate(f=fx, S=Sx, f_range=f_range)

f_new = np.linspace(0.001, 1.0, 256)
f_new = np.log10(f_new)
S_new = slope*f_new + intercept
S_new = 10**S_new
f_new = 10**f_new

S_new = S_new[1:]
f_new = f_new[1:]

foooi=(0.0039**(slope))*(10**(intercept))
fooi=(0.01**(slope))*(10**(intercept))
foi=(0.1**(slope))*(10**(intercept))
fi=(1**(slope))*(10**(intercept))

Sxn_Fit = [foooi,fooi,foi,fi]
f_Fit = [0.0039,0.01,0.1,1.0]

def psdnfit_file(_write=False):
    if _write:
        # creates the text file of calculated spectral density
        psdnfit_file_format = 'psdnfit_'+str(filenumber.zfill(3))+'.txt'
        FileToWrite = open(path_psdnfit + os.sep + psdnfit_file_format, 'w')
        FileToWrite.write("f_Fit" + "\t" + "Sxn_Fit" + "\n" )
        for idx, freq in enumerate(f_Fit):
            FileToWrite.write(str(round(freq, 5)) +"\t"+ str(Sxn_Fit[idx]) + "\n")       
        FileToWrite.close()

psdnfit_file(_write=True) 