"""Created on Sat Oct  5 17:36:19 2019 @author: dadhikar"""

import os
import numpy as np
import pandas as pd


path = '/home/dadhikar/Desktop/__CuIr2S4_V2__/noise_electric field driven/time series/230K/'


def csv_file_read(filename, a, b, c):
    """
    Reading column 'a' and 'b' in filename using read_csv() 
    and returns the columns data as x and y respectively
    """
    dataframe = pd.read_csv(path + os.sep+ filename, delimiter=None, header=None, names=None, 
                index_col=None, usecols=[a, b, c], skiprows=2, skipfooter=0, nrows=None)  
    x1 = dataframe.iloc[:, 0]
    x2 = dataframe.iloc[:, 1] 
    x3 = dataframe.iloc[:, 2]          
    return x1, x2, x3

print('Enter the name of the file:')
print('~'*20)
file_name = input('>>')
print('~'*20)

psd1, psd2, psd3 = csv_file_read(file_name, 9, 10,11) 

file = open(path + os.sep + 'psd.txt', 'w')
file.write('avg' +'\t'+ 'sd' + '\n')

x = psd2
for i in range(len(x)):
    if i%2 == 0:
        avg = np.mean([x[i], x[i+1]])
        sd =  np.std([x[i], x[i+1]])
        file.write(str(avg) +'\t'+ str(sd) + '\n')
    else:
        continue
     
file.close()    