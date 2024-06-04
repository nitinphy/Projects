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
#from psd_slope_alpha import alpha_calculate
#import psd_slope_alpha

# specify the location of the raw data
path_raw = r"C:\Users\Nicholas\Desktop\3-25-2022poly"
# Specifying the path for the file
file_path = r"C:\Users\Nicholas\Desktop\3-25-2022poly\dts_padding_file"
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

# Calculating normalized spectrum

def psd_plot(plot=False):
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 7), dpi=150)
        ax.loglog(fx, Sx, c="r", lw=2.0, label="PSD (Signal Unnormalized)")
        ax.loglog(fy, Sy, c="g", lw=2.0, label="PSD (Background Unnormalized)") #Can comment out backgrund PSD
        ax.set_xlabel(r'f [Hz]', fontsize=15)
        ax.set_ylabel(r'S$_{I}$ [$\frac{I^{2}}{Hz}$]', fontsize=15)
        # ax.set_xlim(f.min(), f.max())
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_title("")
        ax.set_title('Power Spectral Density')
        ax.legend(loc='best')
        plt.show()

psd_plot(plot=False)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import xlrd
from xlrd import open_workbook
from os.path import join

wb= open_workbook(join(path_raw,'NoiseBalanceValues.xls'))
sheet = wb.sheet_by_index(0)

Vosc = sheet.cell_value(int(filenumber), 2)
RB = sheet.cell_value(int(filenumber), 3)
R1 = sheet.cell_value(int(filenumber), 4)


print("")
print("The Value of Vosc is: ", Vosc)
print("The Value of R Balance is: ", RB)
print("The Value of R1/R2 is: ", R1)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gives the oscillating voltage drop in the sample
Vs = Vosc / (1 + (R1 / RB))
# Calculates the normalized spectrum density   
Sx_n = Sx/(Vs**2) 
Sy_n = Sy/(Vs**2)

def psd_n_plot(plot=False):
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 7), dpi=150)
        ax.loglog(fx, Sx_n, c="r", lw=1.5, label="PSD (Signal Normalized)")
        ax.loglog(fy, Sy_n, c="k", lw=1.5, label="PSD (Background Normalized)") #Can comment out backgrund PSD
        ax.set_xlabel(r'f [Hz]', fontsize=14)
        ax.set_ylabel(r'$\frac{S_{R}}{R^{2}}$ [Hz$^{-1}$]', fontsize=14)
        # ax.set_xlim(f.min(), f.max())
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_title("")
        ax.set_title('Normalized Power Spectral Density')
        ax.legend(loc='best')
        plt.show()

psd_n_plot(plot=False)

# alpha_calculate gives frequency exponent, slope 
f_range = [0.0001, 3]
# f_range = [0.01, 8.0]
#slope, intercept, std_err = alpha_calculate(f=fx, S=Sx, f_range=f_range)
slope, intercept, std_err = alpha_calculate(f=fx, S=Sx_n, f_range=f_range)


#--------------------------------------------------------------------------
# Workbook is created
# directory for Noise Analysis Parameters
''' To export data into excel file'''
#import xlwt
#import xlutils
#from xlwt import Workbook
from xlutils.copy import copy 
from os.path import join
from psd_slope_alpha import rsq
path_NoiseData = r"C:\Users\Nicholas\Desktop\3-25-2022poly"

rb=xlrd.open_workbook(join(path_NoiseData,'File1.xls'))
wb = copy(rb)

foooi=(0.0039**(slope))*(10**(intercept))
fooi=(0.01**(slope))*(10**(intercept))
foi=(0.1**(slope))*(10**(intercept))
fi=(1**(slope))*(10**(intercept))


# add_sheet is used to create sheet.
sheet1 = wb.get_sheet('sheet1')

sheet1.write(0, 0, 'File number')
sheet1.write(0, 1, 'Vosc(V)')
sheet1.write(0, 2, 'Balancing Resistance (Ohm)')
sheet1.write(0, 3, 'R1 or R2 (Ohm)')
sheet1.write(0, 4, 'Slope')
sheet1.write(0, 5, 'Error in slope')
sheet1.write(0, 6, 'R Square ')
sheet1.write(0, 7, 'Noise at 3.9 mHz ')
sheet1.write(0, 8, 'Noise at 10 mHz ')
sheet1.write(0, 9, 'Noise at 100 mHz ')
sheet1.write(0, 10, 'Noise at 1 Hz ')

sheet1.write(int(filenumber), 0, filenumber)
sheet1.write(int(filenumber), 1, Vosc)
sheet1.write(int(filenumber), 2, RB)
sheet1.write(int(filenumber), 3, R1)
sheet1.write(int(filenumber), 4, slope)
sheet1.write(int(filenumber), 5, std_err)
sheet1.write(int(filenumber), 6, rsq)
sheet1.write(int(filenumber), 7, foooi)
sheet1.write(int(filenumber), 8, fooi)
sheet1.write(int(filenumber), 9, foi)
sheet1.write(int(filenumber), 10, fi)
#File name of excel file
wb.save(join(path_NoiseData,'File1.xls'))