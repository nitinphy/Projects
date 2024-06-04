"""Initial creation date Tue Nov 13 14:30:27 2018
This code will work for DC method using 2450 Keithely 


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
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
# from file_format_converter import itx_to_txt_converter
from low_pass_kaiser_window_design import  three_stage_decimation
from psd_slope_alpha import alpha_calculate

# location of csv data from the Keithley
path_csv = r"C:\Users\Nicholas\Desktop\1-26-2023\csv_file"
# specify the location of the converted text file
path_text = r"C:\Users\Nicholas\Desktop\1-26-2023\text_file"
# specify the location of the normalized time series
path_norm_ts = r"C:\Users\Nicholas\Desktop\1-26-2023\norm_ts_file"
# directory for detrended time series
path_detrend = r"C:\Users\Nicholas\Desktop\1-26-2023\detrend_file"
# directory for decimated time series
path_DTS = r"C:\Users\Nicholas\Desktop\1-26-2023\dts_file"
# directory for decimated time series AFTER ZERO PADDING
path_DTS_padding = r"C:\Users\Nicholas\Desktop\1-26-2023\dts_padding_file"
# specify the location for the estimated PSD file
path_psd = r"C:\Users\Nicholas\Desktop\1-26-2023\psd_file"
# specify the location for the estimated PSD file
path_psdn = r"C:\Users\Nicholas\Desktop\1-26-2023\psdn_file"
# specify the location for the estimated NORMALIZED PSD file
path_psdnfit = r"C:\Users\Nicholas\Desktop\1-26-2023\psdnfit_file"

for file in os.listdir(path_csv):
    print(file)

#..........................................................................
filenumber = input('Enter the csv file number: ')
#..........................................................................
file_format = 'it_'+str(filenumber).zfill(3)+'.csv'

temperature = input('Enter the temperature of current run: ')
temp = temperature
bias_input = input('Enter the Voltage/Current Supplied (V/A): ')
bias = bias_input

# read csv file as a  pandas dataframe
dataframe = pd.read_csv(path_csv + os.sep+ file_format, delimiter=',',
                            header=None, error_bad_lines=False,
                            skiprows=9, skipfooter=0, nrows=None, low_memory=False)  

# get current time series from the csv file
i = dataframe.iloc[:, 0].to_numpy()

print('Average current: ', i.mean())

iavg = i.mean()

inorm = ((i-iavg)/iavg)

fs = 512.0 # sampling frequency
t_list = []
t0 = 0.0
for n in range(len(inorm)):    
    t_list.append(t0)
    dt = 1/fs
    t0 += dt 
t_list = np.asarray(t_list) 

print('Creating a text file of current time series...')
text_file_format = 'it_'+str(filenumber.zfill(3))+'.txt'
file = open(path_text + os.sep + text_file_format, 'w')
file.write('current(A)' + '\t' + 'time(s)' + '\n')
for m in range(len(t_list)):
    file.write(str(i[m]) +'\t'+ str(t_list[m]) +'\n')  
file.close()
print('Creating text file completed')

print('Now lets calculate PSD from the current time series....')

# t = time
# i = current 
t = np.loadtxt(path_text + os.sep + text_file_format, skiprows=1,
                            usecols=[1], unpack=True)

ipolyfit = np.polyfit(t,inorm,2)

polyfitcurvei = np.polyval(ipolyfit,t)

idetrend = inorm #- polyfitcurvei

#idetrend = signal.detrend(i, axis=-1, type='linear', bp=[0,len(i)])
#idetrend = signal.detrend(i, axis=-1, type='constant', bp=0) # Can choose linear/constant

def detrend_file(_write=False):
    if _write:
        # creates the text file of the detrended time series
        detrend_file_name = 'detrend_'+text_file_format
        FileToWrite = open(path_detrend + os.sep + detrend_file_name, 'w')
        FileToWrite.write("t" + "\t" + "i_detrend" + "\n" )
        for idx, times in enumerate(t):
            FileToWrite.write(str(round(times, 5)) +"\t"+ str(idetrend[idx]) + "\n")       
        FileToWrite.close()
detrend_file(_write=True)

# plot the time series if plot is set True
def plot_signal_time_series(time, inorm, idetrend, plot=False):
    if plot:
        # Plotting time series
        fig, ax = plt.subplots(1, 2, figsize=(7, 7), dpi=100)
        # Remove any trend present on  vx and vy fluctuations.
        ax[0].plot(time, inorm, c="r", label="Before detrending")
        #ax[0].plot(time, vy, c="g", label="Background only")
        ax[0].set_xlabel(r't (s)',  fontsize=10)
        ax[0].set_ylabel(r'$\delta$I (A)', fontsize=10)
        # ax.set_xlim(0.001,1)
        # ax.set_ylim(10e-14, 10e-1)
        ax[0].legend(loc='best')
        #i = signal.detrend(i, axis=-1, type='linear', bp=[41458,458542]) # Can choose linear/constant
        #vy = signal.detrend(vy, axis=-1, type='linear', bp=0)
        ax[1].plot(time, idetrend, c="r", label="After detrending")
        #ax[1].plot(time, vy, c="g", label="Background only")
        ax[1].set_xlabel(r't (s)',  fontsize=10)
        #ax[1].set_ylabel(r'$\delta$I (A)', fontsize=15)
        # ax.set_xlim(0.001,1)
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_facecolor('xkcd:salmon')
        ax[1].legend(loc='best')
        plt.show()

plot_signal_time_series(t, inorm, idetrend, plot=False) #Can set True/False

fs = 512.0  # data sampling frequency
D1 = 4.0  # first stage decimation factor
D2 = 4.0  # second stage decimation factor
D3 = 2.0  # third stage decimation factor
decimation_factor = D1*D2*D3
f_pass_final = (0.5/decimation_factor) * fs
print(r'f_sampling = {} Hz, D1= {}, D2= {}, and D3= {}'.format(fs, D1, D2, D3))

inorm = three_stage_decimation(idetrend, fs, f_pass_final, D1, D2, D3)
t_new = np.linspace(start=0, stop=int(t.max()),
                    num=len(inorm), endpoint=True)

def plot_decimated_signal_time_series(time, inorm, plot=False):
    if plot:
        # Plotting time series
        fig, ax = plt.subplots(1, 1, figsize=(7, 7), dpi=150)
        # Remove any trend present on  vx and vy fluctuations.
        ax.plot(time, inorm, c="r", label="Signal + Background")
        #ax.plot(time, vy, c="g", label="Background only")
        ax.set_xlabel(r't (s)',  fontsize=15)
        ax.set_ylabel(r'$\delta$I (A)', fontsize=15)
        # ax.set_xlim(0.001,1)
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_facecolor('xkcd:salmon')
        ax.set_title('Signal after decimation')
        ax.legend(loc='best')
        plt.show()

plot_decimated_signal_time_series(t_new, inorm, plot=False)

def dts_file(_write=False):
    if _write:
        # creates the text file of calculated spectral density
        #dts_file_name = 'DTS_'+file_to_read
        dts_file_name = 'DTS_'+text_file_format
        FileToWrite = open(path_DTS + os.sep + dts_file_name, 'w')
        FileToWrite.write("t_new" + "\t" + "i" "\n" )
        for idx, times in enumerate(t_new):
            FileToWrite.write(str(round(times, 5)) +"\t"+ str(inorm[idx]) + "\n")       
        FileToWrite.close()
dts_file(_write=True)

n = len(inorm)
m = len(t_new)
# zero padding vx and time in order to make data, power of 2
n1 = int(np.log2(n))
n1 = np.power(2, n1 + 1)
m1 = int(np.log2(m))
m1 = np.power(2, m1 + 1)
padding = np.zeros((n1 - n))
ipad = np.concatenate((inorm, padding))
tpad = np.concatenate((t_new, padding))
#vy = np.concatenate((vy, padding))

def dts_padding_file(_write=False):
    if _write:
        # creates the text file of padded decimates time series
        dts_padding_file_name = 'dts_padding_'+text_file_format
        FileToWrite = open(path_DTS_padding + os.sep + dts_padding_file_name, 'w')
        FileToWrite.write("t_new" + "\t" + "i" "\n" )
        for idx, times in enumerate(tpad):
            FileToWrite.write(str(round(times, 5)) +"\t"+ str(ipad[idx]) + "\n")       
        FileToWrite.close()
dts_padding_file(_write=True)

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
fx, Sx = power_spectrum_density(ipad)
# not selecting zero frequency component
fx = fx[1:]
Sx = Sx[1:]
Sxn = (Sx[1:])/(iavg)

def psd_plot(plot=False):
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 7), dpi=150)
        ax.loglog(fx, Sx, c="r", lw=2.0, label="PSD (unnormalized)")
        #ax.loglog(fy, Sxn, c="g", lw=2.0, label="PSD (normalized)")
        ax.set_xlabel(r'f [Hz]', fontsize=15)
        ax.set_ylabel(r'S$_{I}$ [$\frac{I^{2}}{Hz}$]', fontsize=15)
        # ax.set_xlim(f.min(), f.max())
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_title("")
        ax.set_title('Power Spectral Density')
        ax.legend(loc='best')
        plt.show()

psd_plot(plot=False)

def psd_file(_write=False):
    if _write:
        # creates the text file of calculated spectral density
        psd_file_format = 'psd_'+str(filenumber.zfill(3))+'.txt'
        FileToWrite = open(path_psd + os.sep + psd_file_format, 'w')
        FileToWrite.write("f" + "\t" + "Sx" + "\n" )
        for idx, freq in enumerate(fx):
            FileToWrite.write(str(round(freq, 5)) +"\t"+ str(Sx[idx]) + "\n")       
        FileToWrite.close()

psd_file(_write=True) 

def psdn_plot(plot=False):
    if plot:
        fig, ax = plt.subplots(1, 1, figsize=(7, 7), dpi=150)
        ax.loglog(fx, Sx/(iavg), c="g", lw=2.0, label="PSD (normalized)")
        #ax.loglog(fy, Sxn, c="g", lw=2.0, label="PSD (normalized)")
        ax.set_xlabel(r'f [Hz]', fontsize=15)
        ax.set_ylabel(r'S$_{I}$ [$\frac{1}{Hz}$]', fontsize=15)
        # ax.set_xlim(f.min(), f.max())
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_title("")
        ax.set_title('Normalized Power Spectral Density')
        ax.legend(loc='best')
        plt.show()
  
psdn_plot(plot=False)

def psdn_file(_write=False):
    if _write:
        # creates the text file of calculated normalized spectral density
        psdn_file_format = 'psdn_'+str(filenumber.zfill(3))+'.txt'
        FileToWrite = open(path_psdn + os.sep + psdn_file_format, 'w')
        FileToWrite.write("f" + "\t" + "Sxn" + "\n" )
        for idx, freq in enumerate(fx):
            FileToWrite.write(str(round(freq, 5)) +"\t"+ str((Sx[idx])/(iavg)) + "\n")       
        FileToWrite.close()

psdn_file(_write=True) 

# alpha_calculate gives frequency exponent, slope 
f_range = [0.001, 1]
# f_range = [0.01, 8.0]
slope, intercept, std_err = alpha_calculate(f=fx, S=Sxn, f_range=f_range)
#slope, intercept, std_err = alpha_calculate(f=fx, S=Sx/(iavg**2), f_range=f_range)

f_new = np.linspace(0.001, 1.0, 256)
f_new = np.log10(f_new)
S_new = slope*f_new + intercept
S_new = 10**S_new
f_new = 10**f_new

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

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
''' To export data into excel file'''
#import xlwt
#import xlrd
#import xlutils
#from xlwt import Workbook
import xlrd
#from xlrd import open_workbook
from xlutils.copy import copy 
from os.path import join
from psd_slope_alpha import rsq

# Workbook is created

# directory for Noise Analysis Parameters
path_NoiseData = r"C:\Users\Nicholas\Desktop\1-26-2023"

rb=xlrd.open_workbook(join(path_NoiseData,'Noise_IV_Analysis.xls'))
wb = copy(rb)

foooi=(0.0039**(slope))*(10**(intercept))
fooi=(0.01**(slope))*(10**(intercept))
foi=(0.1**(slope))*(10**(intercept))
fi=(1**(slope))*(10**(intercept))


# add_sheet is used to create sheet.
sheet1 = wb.get_sheet('NoiseData')

sheet1.write(0, 0, 'File number')
sheet1.write(0, 1, 'Avg of Fluc)')
sheet1.write(0, 2, 'Slope')
sheet1.write(0, 3, 'Error in slope')
sheet1.write(0, 4, 'R Square ')
sheet1.write(0, 5, 'Noise at 3.9 mHz ')
sheet1.write(0, 6, 'Noise at 10 mHz ')
sheet1.write(0, 7, 'Noise at 100 mHz ')
sheet1.write(0, 8, 'Noise at 1 Hz ')
sheet1.write(0, 9, 'Temp (K)')
sheet1.write(0, 10, 'Volt/Curr Supplied (V/A)')

sheet1.write(int(filenumber), 0, filenumber)
sheet1.write(int(filenumber), 1, iavg)
sheet1.write(int(filenumber), 2, slope)
sheet1.write(int(filenumber), 3, std_err)
sheet1.write(int(filenumber), 4, rsq)
sheet1.write(int(filenumber), 5, foooi)
sheet1.write(int(filenumber), 6, fooi)
sheet1.write(int(filenumber), 7, foi)
sheet1.write(int(filenumber), 8, fi)
sheet1.write(int(filenumber), 9, temp)
sheet1.write(int(filenumber), 10, bias)
#File name of excel file
wb.save(join(path_NoiseData,'Noise_IV_Analysis.xls'))

print('Data Saved')