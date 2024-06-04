"""Initial creation date Tue Nov 13 14:30:27 2018
   @author: Dasharath Adhikari
   
   Important:
       you need to specify
       1) location of experimental data files = 'path_raw'
       2) location where you want to store the converted files = 'path_text' 
       3) location where you want to stor PSD files = 'path_PSD'
    Also...
       4) can set XXX_file(_write=True) from True to False to not overwrite data
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

from low_pass_kaiser_window_design import  three_stage_decimation
from psd_slope_alpha import alpha_calculate
from file_format_converter import itx_to_txt_converter
#import psd_slope_alpha
base_directory = r"F:\New folder"

path_raw = base_directory
path_csv = base_directory
path_text = os.path.join(base_directory, "text_file")
path_detrend = os.path.join(base_directory, "detrend_file")
path_DTS = os.path.join(base_directory, "dts_file")
path_DTS_padding = os.path.join(base_directory, "dts_padding_file")
path_PSD = os.path.join(base_directory, "psd_file")
path_psdnfit = os.path.join(base_directory, "psdn_file")
plots_folder = os.path.join(base_directory, "plots")

plots_folder = os.path.join(path_csv, "plots")
if not os.path.exists(plots_folder):
    os.makedirs(plots_folder)
    
# Create directories if they don't exist
os.makedirs(path_text, exist_ok=True)
os.makedirs(path_detrend, exist_ok=True)
os.makedirs(path_DTS, exist_ok=True)
os.makedirs(path_DTS_padding, exist_ok=True)
os.makedirs(path_PSD, exist_ok=True)
os.makedirs(path_psdnfit, exist_ok=True)
os.makedirs(plots_folder, exist_ok=True)
    
    

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
#..........................................................................
start_file_number = 1  # Change this to the starting file number
end_file_number = 78   # Change this to the ending file number

for filenumber in range(start_file_number, end_file_number + 1):
    temperature = str(filenumber).zfill(3)
    #bias_input = str(filenumber).zfill(3)

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
    vxavg = vx.mean()
    vyavg = vy.mean()
    
    vxdetrend = vx #- vxavg
    vydetrend = vy#- vyavg
    
    #vxpolyfit = np.polyfit(t,vxdetrend,10)
    #vypolyfit = np.polyfit(t,vy,10)
    
    #polyfitcurvex = np.polyval(vxpolyfit,t)
    #polyfitcurvey = np.polyval(vypolyfit,t)
    
    #vxdetrend = vxdetrend - polyfitcurvex
    #vydetrend = vy - polyfitcurvey
    #vydetrend = signal.detrend(vy, axis=-1, type='linear', bp=[0,len(vy)])
    
    def detrend_file(_write=False):
        if _write:
            # creates the text file of the detrended time series
            detrend_file_name = 'detrend_'+file_to_read
            FileToWrite = open(path_detrend + os.sep + detrend_file_name, 'w')
            FileToWrite.write("t"+str(filenumber).zfill(3) + "\t" + "vx_detrend"+str(filenumber).zfill(3) + 
                              "\t" + "vy_detrend"+str(filenumber).zfill(3) + "\n" )
            for idx, times in enumerate(t):
                FileToWrite.write(str(round(times, 5)) +"\t"+ str(vxdetrend[idx]) +"\t"+
                                  str(vydetrend[idx]) + "\n")       
            FileToWrite.close()
    detrend_file(_write=True)
    
    # plot the time series if plot is set True
    def plot_signal_time_series(time, vx, vy, vxdetrend, vydetrend, plot=False):
        if plot:
            # Plotting time series
            fig, ax = plt.subplots(2, figsize=(10, 10), dpi=150)
            ax[0].plot(time, vx, c="r", label="Signal + Background")
            ax[0].plot(time, vy, c="g", label="Background")
            #ax[0].set_xlabel(r't (s)',  fontsize=15)
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
            file_name_1 = 'time_series_plot_' + str(filenumber) + '.png'
            plt.savefig(os.path.join(plots_folder, file_name_1), dpi=300, bbox_inches='tight')
            plt.close()

    plot_signal_time_series(t, vx, vy, vxdetrend, vydetrend, plot=True) #Can set True/False
    
    fs = 512.0  # data sampling frequency
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
            FileToWrite.write("t_new"+str(filenumber).zfill(3) + "\t" + "vx"+str(filenumber).zfill(3) + "\t" 
                              + "vy"+str(filenumber).zfill(3) + "\n" )
            for idx, times in enumerate(t_new):
                FileToWrite.write(str(round(times, 5)) +"\t"+ str(vx[idx]) +"\t"+
                                  str(vy[idx]) + "\n")       
            FileToWrite.close()
    dts_file(_write=True)
    
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
            file_name_2 = 'plot_decimated_signal_time_series' + str(filenumber) + '.png'
            plt.savefig(os.path.join(plots_folder, file_name_2), dpi=300, bbox_inches='tight')
            plt.close()
    plot_decimated_signal_time_series(t_new, vx, vy, plot=True)
    
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
            FileToWrite.write("t_new"+str(filenumber).zfill(3) + "\t" + "vx"+str(filenumber).zfill(3) + "\t" 
                              + "vy"+str(filenumber).zfill(3) + "\n")
            for idx, times in enumerate(tpad):
                FileToWrite.write(str(round(times, 5)) +"\t"+ str(vx[idx]) + "\t" + str(vy[idx]) + "\n")       
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
    fx, Sx = power_spectrum_density(vx)
    # spectral density estimation form vy data
    fy, Sy = power_spectrum_density(vy)
    # not selecting zero frequency component
    fx = fx[1:]
    fy = fy[1:]
    Sx = Sx[1:]
    Sy = Sy[1:]
    
    def psd_plot(plot=False):
        if plot:
            fig, ax = plt.subplots(1, 1, figsize=(10, 10), dpi=150)
            ax.loglog(fx, Sx, c="r", lw=2.0, marker = ".",label="PSD (signal + background)")
            ax.loglog(fy, Sy, c="g", lw=2.0, label="Background")
            ax.set_xlabel(r'f [Hz]', fontsize=15)
            ax.set_ylabel(r'S$_{v}$ [$\frac{V^{2}}{Hz}$]', fontsize=15)
            # ax.set_xlim(f.min(), f.max())
            # ax.set_ylim(10e-14, 10e-1)
            # ax.set_title("")
            ax.set_title('Power Spectral Density')
            ax.legend(loc='best')
            plt.show()
            file_name_3 = 'psd_plot' + str(filenumber) + '.png'
            plt.savefig(os.path.join(plots_folder, file_name_3), dpi=300, bbox_inches='tight')
            plt.close()
    psd_plot(plot=True)
    
    # Power spectral density of the sample signal
    #S = abs(Sx - Sy)
    # Calculating normalized spectrum
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    import xlrd
    from xlrd import open_workbook
    from os.path import join
    
    wb= open_workbook(join(path_raw,'Noise_IV_Analysis.xls'))
    sheet = wb.sheet_by_index(0)
    
    Temperature = sheet.cell_value(int(filenumber), 1)
    Vosc = sheet.cell_value(int(filenumber), 2)
    RB = sheet.cell_value(int(filenumber), 3)
    R1 = sheet.cell_value(int(filenumber), 4)
    
    
    print("")
    print("The Value of Vosc is: ", Vosc)
    print("The Value of R Balance is: ", RB)
    print("The Value of R1/R2 is: ", R1)
    
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # Gives the oscillating voltage drop in the sample
    Vs = Vosc / ((1 + ((R1 / RB))))
    # Calculates the normalized spectrum density   
    S_n = Sx/(Vs**2) 
    Sy_n = Sy/(Vs**2)
    
    def psd_n_plot(plot=False):
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
            file_name_4 = 'psd_n_plot' + str(filenumber) + '.png'
            plt.savefig(os.path.join(plots_folder, file_name_4), dpi=300, bbox_inches='tight')
            plt.close()
    psd_n_plot(plot=True)
            
    def psd_file(_write=False):
        if _write:
            # creates the text file of calculated spectral density
            psd_file_name = 'PSD_'+file_to_read
            FileToWrite = open(path_PSD + os.sep + psd_file_name, 'w')
            FileToWrite.write("f"+str(filenumber).zfill(3) + "\t" + "Sx"+str(filenumber).zfill(3) + "\t" 
                              + "Sy"+str(filenumber).zfill(3) + "\t" +"S_n"+str(filenumber).zfill(3) + "\t" 
                              + "Sy_n"+str(filenumber).zfill(3) + "\n")
            for idx, freq in enumerate(fx):
                FileToWrite.write(str(round(freq, 5)) +"\t"+ str(Sx[idx]) +"\t"+
                                  str(Sy[idx]) + "\t"+ str(S_n[idx]) + "\t" + str(Sy_n[idx]) + "\n")       
            FileToWrite.close()
    psd_file(_write=True)
          
    # alpha_calculate gives frequency exponent, slope 
    # f_range = [0.01, 8.0]
    f_range = [0.0001, 1.0]
    slope, intercept, std_err = alpha_calculate(f=fx, S=S_n, f_range=f_range)
    
    f_new = np.linspace(0.001, 1.0, 256)
    f_new = np.log10(f_new)
    S_new = slope*f_new + intercept
    S_new = 10**S_new
    f_new = 10**f_new
    
    foooi=(0.0039**(slope))*(10**(intercept))
    fooi=(0.01**(slope))*(10**(intercept))
    foi=(0.1**(slope))*(10**(intercept))
    fi=(1**(slope))*(10**(intercept))
    
    S_n_Fit = [foooi,fooi,foi,fi]
    f_Fit = [0.0039,0.01,0.1,1.0]
    
    #print('*'*20)
    #print('PSD at f = 0.01 is {}'.format(S_n[np.where(fx==0.01)]))
    #print('PSD at f = 0.1 is {}'.format(S_n[np.where(fx==0.1)]))
    #print('PSD at f = 1.0 is {}'.format(S_n[np.where(fx==1.0)]))
    #print('*'*20)
    #sys.exit()
    
    def psdnfit_file(_write=False):
        if _write:
            # creates the text file of calculated spectral density
            psdnfit_file_format = 'psdnfit_'+str(filenumber).zfill(3)+'.txt'
            FileToWrite = open(path_psdnfit + os.sep + psdnfit_file_format, 'w')
            FileToWrite.write("f_Fit"+str(filenumber).zfill(3) + "\t" + "Sxn_Fit"+str(filenumber).zfill(3) + "\n" )
            for idx, freq in enumerate(f_Fit):
                FileToWrite.write(str(round(freq, 5)) +"\t"+ str(S_n_Fit[idx]) + "\n")       
            FileToWrite.close()
    
    psdnfit_file(_write=True) 
    
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    ''' To export data into excel file'''
    #import xlwt
    #import xlrd
    #import xlutils
    #from xlwt import Workbook
    from xlutils.copy import copy 
    from os.path import join
    from psd_slope_alpha import rsq
    
    # Workbook is created
    
    # directory for Noise Analysis Parameters
    path_NoiseData = r"F:\New folder"
    
    rb=xlrd.open_workbook(join(path_NoiseData,'Noise_IV_Analysis.xls'))
    wb = copy(rb)
    
    foooi=(0.0039**(slope))*(10**(intercept))
    fooi=(0.01**(slope))*(10**(intercept))
    foi=(0.1**(slope))*(10**(intercept))
    fi=(1**(slope))*(10**(intercept))
    
    
    # add_sheet is used to create sheet.
    sheet1 = wb.get_sheet('NoiseData')
    
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
    sheet1.write(0, 11, 'Temperature ')
    
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
    sheet1.write(int(filenumber), 11, Temperature)
    #File name of excel file
    wb.save(join(path_NoiseData,'Noise_IV_Analysis.xls'))