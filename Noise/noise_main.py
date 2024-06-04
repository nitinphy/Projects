"""Initial creation date Tue Nov 13 14:30:27 2018
   @author: Dasharath Adhikari
   
   Important:
       you need to specify
       1) location of experimental data files = 'path_raw'
       2) location where you want to store the converted files = 'path_text' 
       3) location where you want to stor PSD files = 'path_PSD'
   
"""

import os
import sys
from scipy import signal
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
from file_format_converter import itx_to_txt_converter
from low_pass_kaiser_window_design import  three_stage_decimation
from psd_slope_alpha import alpha_calculate

# specify the location of the raw data
path_raw = r"C:\Users\Nicholas\Desktop\12-11-2020"
# specify the file location for the converted text file
path_text = r'C:\Users\Nicholas\Desktop\12-11-2020\Text_file'


#print(os.listdir(path_text))                  
file_convert = False

if file_convert:
  print("Converting file format")
  itx_to_txt_converter(raw_data_path=path_raw, text_data_path=path_text)


file_list = sorted(os.listdir(path_text))
print('-'*25)
print('The file list: ', file_list)
print('-'*25)

#sys.exit() #Set file_convert True/False to just convert files

file_to_read = input('Enter the name of the file to read data from >>  ')


# reading time, vx and vy from the file
def read_file(filepath, filename, a, b, c):
    x0, x1, x2 = np.loadtxt(filepath + os.sep + filename, skiprows=1,
                            usecols=[a, b, c], unpack=True)
    return x0, x1, x2


# t = time
# vx = in-phase fluctuation signal
# vy = out-of-phase fluctuation signal
t, vx, vy = read_file(path_text, file_to_read, 0, 1, 2)


# plot the time series if plot is set True
def plot_signal_time_series(time, vx, vy, plot=False):
    if plot:
        # Plotting time series
        fig, ax = plt.subplots(1, 2, figsize=(10, 10), dpi=150)
        # Remove any trend present on  vx and vy fluctuations.
        ax[0].plot(time, vx, c="r", label="Signal + Background")
        ax[0].plot(time, vy, c="g", label="Background only")
        ax[0].set_xlabel(r't (s)',  fontsize=15)
        ax[0].set_ylabel(r'$\delta$V (V)', fontsize=15)
        # ax.set_xlim(0.001,1)
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_facecolor('xkcd:salmon')
        ax[0].set_title('Signal before detrending')
        ax[0].legend(loc='best')
        vx = signal.detrend(vx, axis=-1, type='constant', bp=0) # Can choose linear/constant
        vy = signal.detrend(vy, axis=-1, type='constant', bp=0)
        ax[1].plot(time, vx, c="r", label="Signal + Background")
        ax[1].plot(time, vy, c="g", label="Background only")
        ax[1].set_xlabel(r't (s)',  fontsize=15)
        ax[1].set_ylabel(r'$\delta$V (V)', fontsize=15)
        # ax.set_xlim(0.001,1)
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_facecolor('xkcd:salmon')
        ax[1].set_title('Signal after detrending')
        ax[1].legend(loc='best')
        plt.show()


plot_signal_time_series(t, vx, vy, plot=True) #Can set True/False

fs = 256.0  # data sampling frequency
D1 = 4.0  # first stage decimation factor
D2 = 2.0  # second stage decimation factor
D3 = 2.0  # third stage decimation factor
decimation_factor = D1*D2*D3
f_pass_final = (0.5/decimation_factor) * fs
print(r'f_sampling = {} Hz, D1= {}, D2= {}, and D3= {}'.format(fs, D1, D2, D3))


# first detrend and do threee stage (antialising + decimation)
# of vx and vy signal
vx = three_stage_decimation(vx, fs, f_pass_final, D1, D2, D3)
vy = three_stage_decimation(vy, fs, f_pass_final, D1, D2, D3)
t_new = np.linspace(start=0, stop=int(t.max()),
                    num=len(vx), endpoint=True)


def plot_decimated_signal_time_series(time, vx, vy, plot=False):
    if plot:
        # Plotting time series
        fig, ax = plt.subplots(1, 1, figsize=(10, 10), dpi=150)
        # Remove any trend present on  vx and vy fluctuations.
        ax.plot(time, vx, c="r", label="Signal + Background")
        ax.plot(time, vy, c="g", label="Background only")
        ax.set_xlabel(r't (s)',  fontsize=15)
        ax.set_ylabel(r'$\delta$V (V)', fontsize=15)
        # ax.set_xlim(0.001,1)
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_facecolor('xkcd:salmon')
        ax.set_title('Signal after decimation')
        ax.legend(loc='best')
        plt.show()


plot_decimated_signal_time_series(t_new, vx, vy, plot=True)

n = len(vx)
# zero padding vx and vy in order to make data, power of 2
n1 = int(np.log2(n))
n1 = np.power(2, n1 + 1)
padding = np.zeros((n1 - n))
vx = np.concatenate((vx, padding))
vy = np.concatenate((vy, padding))

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
fx, Sx = power_spectrum_density(vx)
# spectral density estimation form vy data
fy, Sy = power_spectrum_density(vy)

# not selecting zero frequency component
fx = fx[1:]
fy = fy[1:]
Sx = Sx[1:]
Sy = Sy[1:]


def psd_plot(plot=True):
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10), dpi=150)
        ax.loglog(fx, Sx, c="r", lw=2.0, label="PSD (signal + background)")
        ax.loglog(fy, Sy, c="g", lw=2.0, label="Background only")
        ax.set_xlabel(r'f [Hz]', fontsize=15)
        ax.set_ylabel(r'S$_{v}$ [$\frac{V^{2}}{Hz}$]', fontsize=15)
        # ax.set_xlim(f.min(), f.max())
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_title("")
        ax.set_title('Power Spectral Density')
        ax.legend(loc='best')
        plt.show()

psd_plot(plot=True)



# Power spectral density of the sample signal
# S = abs(Sx - Sy)  
# Calculating normalized spectrum
Vosc = float(input("Enter the value of oscillating voltage >>  "))
RB = float(input("Enter the value of balancing resistance >>  "))
R1 = float(input("Enter the value of R1/R2 resistance >>  "))
# Gives the oscillating voltage drop in the sample
Vs = Vosc / (1 + (R1 / RB))
# Calculates the normalized spectrum density   
S_n = Sx/(Vs**2) 
Sy_n = Sy/(Vs**2)

def psd_n_plot(plot=True):
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 7), dpi=150)
        ax.loglog(fx, S_n, c="r", lw=1.5, label="PSD (X-channel)")
        ax.loglog(fy, Sy_n, c="k", lw=1.5, label="PSD (Y-channel)") #Can comment out backgrund PSD
        ax.set_xlabel(r'f [Hz]', fontsize=14)
        ax.set_ylabel(r'$\frac{S_{R}}{R^{2}}$ [Hz$^{-1}$]', fontsize=14)
        # ax.set_xlim(f.min(), f.max())
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_title("")
        ax.set_title('Normalized Power Spectral Density')
        ax.legend(loc='best')
        plt.show()


psd_n_plot(plot=True)
              

# alpha_calculate gives frequency exponent, slope 
f_range = [fx.min(), 0.95]
# f_range = [0.01, 8.0]
slope, intercept, std_err = alpha_calculate(f=fx, S=S_n, f_range=f_range)

print('*'*20)
print('PSD at f = 0.01 is {}'.format(S_n[np.where(fx==0.01)]))
print('PSD at f = 0.1 is {}'.format(S_n[np.where(fx==0.1)]))
print('PSD at f = 1.0 is {}'.format(S_n[np.where(fx==1.0)]))
print('*'*20)
#sys.exit()

f_new = np.linspace(0.001, 1.0, 256)
f_new = np.log10(f_new)
S_new = slope*f_new + intercept
S_new = 10**S_new
f_new = 10**f_new

# directory for estimted power spectral density files
path_PSD = r"C:\Users\Nicholas\Desktop\12-11-2020\psd_file"

def psd_file(_write=True):
    if _write:
        # creates the text file of calculated spectral density
        psd_file_name = 'PSD_'+file_to_read
        FileToWrite = open(path_PSD + os.sep + psd_file_name, 'w')
        FileToWrite.write("f" + "\t" + "Sx" "\t" + "S_n" + "\n" )
        for idx, freq in enumerate(fx):
            FileToWrite.write(str(round(freq, 5)) +"\t"+ str(Sx[idx]) +"\t"+
                              str(S_n[idx]) + "\n")       
        FileToWrite.close()
psd_file(_write=True)
    