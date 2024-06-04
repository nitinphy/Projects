# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 14:06:34 2022

@author: Nicholas
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
path_raw = r"C:\Users\Nicholas\Desktop\3-25-2022"
# specify the file location for the converted text file
path_text = r'C:\Users\Nicholas\Desktop\3-25-2022\Text_file'
# directory for detrended time series
path_detrend = r"C:\Users\Nicholas\Desktop\3-25-2022\detrend_file"
# directory for decimated time series
path_DTS = r"C:\Users\Nicholas\Desktop\3-25-2022\dts_file"
# directory for decimated time series AFTER ZERO PADDING
path_DTS_padding = r"C:\Users\Nicholas\Desktop\3-25-2022\dts_padding_file"

#..........................................................................
filenumber = input('Enter the text file number: ')
#..........................................................................
file_to_read = 'VT'+str(filenumber).zfill(3)+'.txt'

# reading time, vx and vy from the file
def read_file(filepath, filename, a, b, c):
    x0, x1, x2 = np.loadtxt(filepath + os.sep + filename, skiprows=1,
                            usecols=[a, b, c], unpack=True)
    return x0, x1, x2

# t = time
# vx = in-phase fluctuation signal
# vy = out-of-phase fluctuation signal
t, vx, vy = read_file(path_text, file_to_read, 0, 1, 2)

# Remove any trend present on  vx and vy fluctuations.
# Can choose type=='linear' or type=='constant'
vxdetrend = signal.detrend(vx, axis=-1, type='linear', bp=[0,len(vx)])
vydetrend = signal.detrend(vy, axis=-1, type='linear', bp=[0,len(vy)])

def detrend_file(_write=False):
    if _write:
        # creates the text file of the detrended time series
        detrend_file_name = 'detrend_'+file_to_read
        FileToWrite = open(path_detrend + os.sep + detrend_file_name, 'w')
        FileToWrite.write("t" + "\t" + "vx_detrend" "\t" + "vy_detrend" + "\n" )
        for idx, times in enumerate(t):
            FileToWrite.write(str(round(times, 5)) +"\t"+ str(vxdetrend[idx]) +"\t"+
                              str(vydetrend[idx]) + "\n")       
        FileToWrite.close()
detrend_file(_write=False)

# plot the time series if plot is set True
def plot_signal_time_series(time, vx, vy, vxdetrend, vydetrend, plot=False):
    if plot:
        # Plotting time series
        fig, ax = plt.subplots(1, 2, figsize=(10, 10), dpi=150)
        ax[0].plot(time, vx, c="r", label="Signal + Background")
        ax[0].plot(time, vy, c="g", label="Background only")
        ax[0].set_xlabel(r't (s)',  fontsize=15)
        ax[0].set_ylabel(r'$\delta$V (V)', fontsize=15)
        # ax.set_xlim(0.001,1)
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_facecolor('xkcd:salmon')
        ax[0].set_title('Signal before detrending')
        ax[0].legend(loc='best')
        ax[1].plot(time, vxdetrend, c="r", label="Signal + Background")
        #ax[1].plot(time, vydetrend, c="g", label="Background only")
        ax[1].set_xlabel(r't (s)',  fontsize=15)
        ax[1].set_ylabel(r'$\delta$V (V)', fontsize=15)
        # ax.set_xlim(0.001,1)
        # ax.set_ylim(10e-14, 10e-1)
        # ax.set_facecolor('xkcd:salmon')
        ax[1].set_title('Signal after detrending')
        ax[1].legend(loc='best')
        plt.show()

plot_signal_time_series(t, vx, vy, vxdetrend, vydetrend, plot=True) #Can set True/False

fs = 256.0  # data sampling frequency
D1 = 4.0  # first stage decimation factor
D2 = 2.0  # second stage decimation factor
D3 = 2.0  # third stage decimation factor
decimation_factor = D1*D2*D3
f_pass_final = (0.5/decimation_factor) * fs
print(r'f_sampling = {} Hz, D1= {}, D2= {}, and D3= {}'.format(fs, D1, D2, D3))

# first detrend and do threee stage (antialising + decimation)
# of vx and vy signal
vx = three_stage_decimation(vxdetrend, fs, f_pass_final, D1, D2, D3)
vy = three_stage_decimation(vydetrend, fs, f_pass_final, D1, D2, D3)
t_new = np.linspace(start=0, stop=int(t.max()),
                    num=len(vx), endpoint=True)

def dts_file(_write=False):
    if _write:
        # creates the text file of the decimated time series after detrending
        dts_file_name = 'DTS_'+file_to_read
        FileToWrite = open(path_DTS + os.sep + dts_file_name, 'w')
        FileToWrite.write("t_new" + "\t" + "vx" "\t" + "vy" + "\n" )
        for idx, times in enumerate(t_new):
            FileToWrite.write(str(round(times, 5)) +"\t"+ str(vx[idx]) +"\t"+
                              str(vy[idx]) + "\n")       
        FileToWrite.close()
dts_file(_write=False)

def plot_decimated_signal_time_series(time, vx, vy, plot=False):
    if plot:
        # Plotting decimated time series
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

plot_decimated_signal_time_series(t_new, vx, vy, plot=False)

# zero padding vx and vy in order to make data, power of 2
n = len(vx)
n1 = int(np.log2(n))
n1 = np.power(2, n1 + 1)
padding = np.zeros((n1 - n))
vx = np.concatenate((vx, padding))
vy = np.concatenate((vy, padding))
tpad = np.concatenate((t_new, padding))

def dts_padding_file(_write=False):
    if _write:
        # creates the text file of the zero padded decimates time series
        dts_padding_file_name = 'dts_padding_'+file_to_read
        FileToWrite = open(path_DTS_padding + os.sep + dts_padding_file_name, 'w')
        FileToWrite.write("t_new" + "\t" + "vx" "\t" + "vy" + "\n")
        for idx, times in enumerate(tpad):
            FileToWrite.write(str(round(times, 5)) +"\t"+ str(vx[idx]) + "\t" + str(vy[idx]) + "\n")       
        FileToWrite.close()
dts_padding_file(_write=False)