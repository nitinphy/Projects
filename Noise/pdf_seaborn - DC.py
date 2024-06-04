"""Created on Fri Apr 26 17:16:40 2019  @ author: dasharath"""

import os
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
import numpy as np
import seaborn as sb

# Specifying the path for the file
file_path = r"X:\Personal Data\Nbo2\09 Sept\09262023\224\Noise pos\16min\dts_file"
# specify the location for the estimated PDF file
path_pdf = r"X:\Personal Data\Nbo2\09 Sept\09262023\224\Noise pos\16min\pdf_file"

for file in os.listdir(file_path):
    print(file)
    
start_file_number = 1  # Change this to the starting file number
end_file_number = 129   # Change this to the ending file number


def read_vx(filename):
    vx = np.loadtxt(file_path + os.sep + filename, skiprows=1, usecols=1)
    # calculating average and standard deviation for time series
    avg = np.mean(vx)
    sd = np.std(vx)
    vx -= avg
    vx /= sd
    return vx
for filenumber in range(start_file_number, end_file_number + 1):
   
    file = 'DTS_it_'+str(filenumber).zfill(3)+'.txt'
    
    vx = read_vx(file)
    kde_kws={"color": "red", "lw": 3}
    hist_kws={"histtype": "stepfilled", "linewidth": 1,"alpha": 1, "color": "g"}
    hist_kws = np.array(hist_kws)    
    fig, ax = plt.subplots(1,1, figsize=(7,7), sharex=False, facecolor="white", dpi= 170)
    ax = sb.distplot(vx, bins=None, hist=True, kde=True, color='C6',fit=None, hist_kws= None,
                             norm_hist=True, kde_kws= kde_kws, label = r'PDF')
    ax.set_xlabel(r"$\frac{\delta R - \overline{\delta R}}{\sigma_{\delta R}}$")
    ax.set_ylabel(r"Estimated PDF")
    ax.legend(loc='lower center', bbox_to_anchor=(0.80, 0.85), shadow=False, ncol=1,  
               fontsize='small')
    #ax.set_xticks([-2, -1, 0, 1, 2])
    #ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
    #ax.plot(vx, stats.norm.pdf(vx1, avg, sd))
    # Save the plot with the file number as the filename
    filenumbers = str(filenumber)
    plot_filename = os.path.join(path_pdf, f'plot_{filenumbers.zfill(3)}.png')
    plt.savefig(plot_filename)  # Save the plot as an image file
    
    plt.close()  # Close the
    
    
    #hist_kws = hist_kws[1:]
    #vx = vx[1:]
    
   # def pdf_file(_write=False):
       # if _write:
            # creates the text file of calculated spectral density
            
           # pdf_file_format = 'pdf_'+str(filenumbers.zfill(3))+'.txt'
           #FileToWrite.write("Bin" + "\t" + "PDF_Estimation" + "\n" )
            #for idx, freq in enumerate(vx):
                #FileToWrite.write(str(round(freq, 5)) +"\t"+ str(hist_kws[idx]) + "\n")       
           # FileToWrite.close()
    
    #pdf_file(_write=True) 
