# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 13:48:22 2018

@author: das_d
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats.kde import gaussian_kde
fig, ax = plt.subplots(figsize=(6, 6), dpi=50)
# read first column from file and store as an time_sec array
time_sec = np.loadtxt('VT001.txt', usecols=0)
# delete first 120000 elements from array and store as time_sec1 array
# time_sec1 = np.delete(time_sec, np.s_[0:120000])
# change time from sec to min
time_m = (1/60.)*time_sec
V_x = np.loadtxt('VT001.txt', usecols=1)
# V_xnew = np.delete(V_x, np.s_[0:120000])
# calculate means of V_xnew
mean = np.mean(V_x)
variance = np.var(V_x)
std = np.sqrt(variance)

# Generate new V_x1 array of values with zero mean value
V_x1 = V_x - mean
# print(len(V_x1))
# Estimating the pdf and plotting it
KDEpdf = gaussian_kde(V_x1)
x = np.linspace(-(mean-3*std), (mean+3*std), 2000)
ax.plot(x, KDEpdf(x), 'r', label="KDE estimation", color="blue")
ax.legend()
plt.show()
