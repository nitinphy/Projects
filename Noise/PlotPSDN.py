# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 14:17:54 2022

@author: Nicholas
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

# Specifying the path for the file
file_path = r"C:\Users\Nicholas\Desktop\3-25-2022\dts_padding_file"
#..........................................................................
filenumber = input('Enter the txt file number: ')
#..........................................................................
file_to_read = 'dts_padding_VT'+str(filenumber).zfill(3)+'.txt'

# reading time, vx and vy from the file
def read_dts(filepath, filename, a, b, c):
    x0, x1, x2 = np.loadtxt(filepath + os.sep + filename, skiprows=1,
                            usecols=[a, b, c], unpack=True)
    return x0, x1, x2

time, fluctuation_x, fluctuation_y = read_dts(file_path, file_to_read, 0, 1, 2)

#avgfluctuation = fluctuation.mean()

fs = 256.0  # data sampling frequency
D1 = 4.0  # first stage decimation factor
D2 = 2.0  # second stage decimation factor
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
fx, Sx = power_spectrum_density(fluctuation_x)
# spectral density estimation form vy data
fy, Sy = power_spectrum_density(fluctuation_y)
# not selecting zero frequency component
fx = fx[1:]
fy = fy[1:]
Sx = Sx[1:]
Sy = Sy[1:]
#Sxn = Sx[1:]/(1.3078**2)

psd_n_plot(plot=True)

# alpha_calculate gives frequency exponent, slope 
f_range = [0.001, 4.0]
# f_range = [0.01, 8.0]
#slope, intercept, std_err = alpha_calculate(f=fx, S=Sx, f_range=f_range)
slope, intercept, std_err = alpha_calculate(f=fx, S=Sx_n, f_range=f_range)