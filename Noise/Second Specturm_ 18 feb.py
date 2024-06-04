# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 11:03:55 2024

@author: Nitin Kumar
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal import butter, lfilter
import pandas as pd
import glob
import os

def calculate_power_spectrum(signal, fs):
    # Number of data points
    N = len(signal)
    
    # Perform FFT
    freqs = fftfreq(N, 1/fs)
    fft_vals = fft(signal)
    
    # Compute power spectrum
    power_spectrum = np.abs(fft_vals)**2 / N
    
    return freqs, power_spectrum

def integrate_power_spectrum(freqs, power_spectrum, bandpass):
    # Apply bandpass filter
    b, a = butter(4, bandpass, btype='bandpass', fs=fs)
    filtered_power_spectrum = lfilter(b, a, power_spectrum)
    
    # Integrate over chosen bandpass
    integrated_power = np.trapz(filtered_power_spectrum, freqs)
    
    return integrated_power

# Folder path
folder_path = r'F:\NbO2 manuscript\Data\325 current driven noise\Unormalized Noise\16 Min\TEST FOR SS'

# Get a list of all CSV files in the folder
csv_files = glob.glob(os.path.join(folder_path, '*001.csv'))

# Loop over all CSV files
for file_path in csv_files:
    # Read the CSV file
    df = pd.read_csv(file_path)

    # Assuming the signal is in the first column, if not, change it to the correct column name
    V = df.iloc[:, 1].values

    # Example parameters
    fs = 500  # Sampling frequency (Hz)
    N1 = 4800   # Number of points per segment
    N2 = 100    # Number of segments
    T1 = N1 / fs  # Time duration of each segment
    T = N2 / fs   # Total data collection time

    # Break signal into segments
    segments = np.split(V, N2)

    # Calculate power spectrum for each segment
    power_spectrum_segments = []
    for segment in segments:
        freqs, power_spectrum = calculate_power_spectrum(segment, fs)
        power_spectrum_segments.append(power_spectrum)

    # Integrate power spectrum over a chosen bandpass for each segment
    bandpass = (10, 100)  # Example bandpass (Hz)
    integrated_powers = []
    for power_spectrum in power_spectrum_segments:
        integrated_power = integrate_power_spectrum(freqs, power_spectrum, bandpass)
        integrated_powers.append(integrated_power)

    # Calculate second spectrum
    freqs_second_spectrum, power_spectrum_second_spectrum = calculate_power_spectrum(integrated_powers, 1/T1)

    # Line and points plot for power spectrum of each segment (positive frequencies only) on log-log scale
    plt.figure(figsize=(10, 6))
    for i, power_spectrum in enumerate(power_spectrum_segments):
        plt.plot(freqs[:N1//2], power_spectrum[:N1//2], marker='o', linestyle='-', markersize=2, label=f'Segment {i+1}')  # Only positive frequencies
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power')
    plt.title('Power Spectrum of Segments (Positive Frequencies, Log-Log Scale)')
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(folder_path, f'Power_Spectrum_{os.path.basename(file_path)}.png'))
    plt.show()
    print(power_spectrum_second_spectrum)
    # Line and points plot for second spectrum (positive frequencies only) on log-log scale
    plt.figure(figsize=(10, 6))
    plt.plot(freqs_second_spectrum[:N2//2], power_spectrum_second_spectrum[:N2//2], marker='o', linestyle='-', markersize=2)  # Only positive frequencies
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power')
    plt.title('Second Spectrum (Positive Frequencies, Log-Log Scale)')
    plt.grid(True)
    plt.savefig(os.path.join(folder_path, f'Second_Spectrum_{os.path.basename(file_path)}.png'))
    plt.show()
