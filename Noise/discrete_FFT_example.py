import numpy as np
from math import pi, sin
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import fftpack
from scipy import signal


def signal_sample(t):
    # Generate 1Hz and 22Hz sine signal with noise background
    f_1 = 2*np.sin(2*pi*t)
    f_2 = 3*np.sin(22*2*pi*t)
    f_3 = 2*np.random.randn(*np.shape(t))
    return f_1 + f_2 + f_3


B = 30.0   # Band frequnecy
f_s = 80.0      # Sampling frequency
delta_f = 0.01  # Resolution in spectrum
N = int(f_s/delta_f)  # Number of sample point taken
T = N/f_s            # Sampling time period

# Generates function in time space
t = np.linspace(0, T, N)
f_t = signal_sample(t)


# Plotting function
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(t, f_t)
ax[0].set_xlabel('Time(s)')
ax[0].set_ylabel('Sample_Signal')
ax[1].plot(t, f_t)
ax[1].set_xlim(0, 10)
ax[1].set_xlabel('Time(s)')

# fourier transformation
F = fftpack.fft(f_t)
f = fftpack.fftfreq(N, 1.0/f_s)
mask = np.where(f >= 0)
fig, ax = plt.subplots(3, 1, figsize=(10, 10))
ax[0].plot(f[mask], np.log((F[mask])**2), label="Real")
# ax[0].plot(B, 0, 'r*', markersize=10)
ax[0].set_ylabel(r'$\log(|F|)$', fontsize=15)
ax[1].plot(f[mask], abs(F[mask])/N, label="Real")
ax[1].set_xlim(0.5, 1.5)
ax[1].set_ylabel('$|F|$', fontsize=15)
ax[2].plot(f[mask], abs(F[mask])/N, label="Real")
ax[2].set_xlim(21, 23)
ax[2].set_ylabel('$|F|$', fontsize=15)
ax[2].set_xlabel('Frequency(Hz)', fontsize=15)
