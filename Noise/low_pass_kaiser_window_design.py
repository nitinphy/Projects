"""
Created on Sat Nov 24 11:43:57 2018 
@author: Dasharath Adhikari
"""

from __future__ import division
import numpy as np
import math
from scipy import signal


def kaiser_beta(ripple):
    
    """
    delta_factor : float
        The desired attenuation in the stopband and maximum ripple in
        the passband.For example delta_factor = 0.001. 
        This should be a *positive* number.
    """ 
    ripple = float(ripple)
    a = -20*math.log10(ripple)       # Attenuation in dB
    
    if a > 50:
        beta = 0.1102 * (a - 8.7)
    elif a > 21:
        beta = 0.5842 * (a-21) ** 0.4 + 0.07886 * (a-21)
    else:
        beta = 0.0
        
    return a, beta


def transition_band(fp, fs, f_S):
    
    """
    fp = pass band frequency
    fs = stop band frequency
    f_S = sampling frequency
    fn = Nyquist frequency  = 0.5 * f_S
    
    """
    TB = (fs - fp)/ (0.5 * f_S)
    
    return TB


def kaiser_window_length(ripple, fp, fs, f_s):
    
    """
    For given ripple and the normalized transition width (fs-fp)/fN,
    fp = pass band frequency
    fs = stop band frequency
    f_S = sampling frequency
    fn = Nyquist frequency  = 0.5 * f_s
    returns the length for the Kaiser window's impulse response function.
    
    """
    a = -20 * np.log10(ripple)
    TB = (fs - fp)/ (0.5 * f_s) # Normalized transition band
    M = int(np.ceil((a - 8)/(2.285*np.pi*TB))) # Order of the Kaiser window
    N = M + 1                                  # Length of the Kaiser window
    
    return N 

def kaiser_FIR(ripple, f_s, fp, fs):
    """
    ripple : float
    (The desired attenuation in the stopband and maximum ripple in
    the passband.For example delta_factor = 0.001.
    This should be a *positive* number.)
    For given ripple and the normalized transition width (fs-fp)/fN,
    fp = pass band frequency
    fs = stop band frequency
    f_s = sampling frequency
    fn = Nyquist frequency  = 0.5 * f_s
    returns the length for the Kaiser window's impulse response function.
    """
    ripple = float(ripple)
    a = -20 * math.log10(ripple)  # Attenuation in dB
    if a > 50:
        beta = 0.1102 * (a - 8.7)
    elif a > 21:
        beta = 0.5842 * (a-21) ** 0.4 + 0.07886 * (a-21)
    else:
        beta = 0.0
    # Normalized transition band
    TB = (fs - fp) / (f_s)
    # cut-off frequency and is always defined
    # as the fraction of sampling frequency
    fc = (fp + fs) / (2 * f_s)
    # order of the Kaiser window
    M = int((a - 8) / (2.285 * 2 * np.pi * TB))        
    # Length of the Kaiser window
    N = M + 1                                         
    alpha = M/2
    # create 1-D array of lenght N
    n = np.arange(N)                                
    # Sinc filter
    h = np.sinc(2 * fc * (n - (alpha / 2)))
    # Kaiser window function
    w = np.i0( beta * np.sqrt (1 - ((n/alpha) - 1) ** 2)) / np.i0 (beta)
    h = h*w # Windowed sinc-filter 
    h = h/np.sum(h)
    return h  # Returns an windnow-sinc filter of length N. 

def three_stage_decimation(x, f_s, f_pass_final, D1, D2, D3):
    
    """
    The initial sampling rate is f_s Hz.
    Downsample the signal by D1*D2*D3
    I will be doing three stages of signal decimation
    (antialiasing filtering + decimation) with  
    downsampling factors of D1, D2, and D3 respectively. 
    """ 
    # detrend the signal first
    x = signal.detrend(x, axis=-1, type='linear', bp=0)
    # ---First Stage of Signal Decimation------------------------------
    # ripple = 0.01
    fs = f_s/2   # stop band
    window = kaiser_FIR(0.01, f_s, f_pass_final, fs)
    x = np.convolve(x, window)
    x = x[::int(D1)]
    f_s = f_s / D1     # new sampling frequency
    # ---Second Stage of Signal Decimation -----------------------------
    # ripple = 0.01
    fs = f_s/2
    window = kaiser_FIR(0.01, f_s, f_pass_final, fs)
    x = np.convolve(x, window)
    x = x[::int(D2)]
    f_s = f_s/D2  # new sampling frequency   
    # ----Third Stage of Signal Decimation-------------------------------
    # ripple = 0.01
    fs = f_s/2
    window = kaiser_FIR(0.01, f_s, f_pass_final, fs)
    x = np.convolve(x, window)
    x = x[::int(D3)]
    # final sampling frequency will be f_s /D3
    return x
