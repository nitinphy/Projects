#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Author @ Dasharath Adhikari
"""
import os
import sys
import pandas as pd 
import numpy as np

# file location
pathin = r"C:\Users\Nicholas\Desktop\Kappa-CNC\csv_file"
pathout = r"C:\Users\Nicholas\Desktop\Kappa-CNC\text_file"

for file in os.listdir(pathin):
    print(file)
    
    
# sys.exit()    

def csv_file_read(filename, a):
    """
    Reading column of filename using read_csv() 
    and return them
    """
    dataframe = pd.read_csv(pathin + os.sep+ filename, delimiter=None, header=None, names=None, 
                            index_col=None, usecols=[a], skiprows=0, skipfooter=0, nrows=None)  
    x1 = dataframe.iloc[:, 0]           
    return x1 
#..............................................................................................
input_file = input('Please enter the csv filename you want to convert: ')
#..............................................................................................
i = csv_file_read(input_file, 0)

print('Lenth of current time series: ', len(i))
# print(len(i))

t_list = []
t0 = 0.0
for n in range(len(i)):    
    t_list.append(t0)
    dt = 0.001953125
    t0 += dt 
t_list = np.asarray(t_list) 


print('Creating a text file...')
# creating new file and write i  and t 
#..............................................................................................
write_into_file = input('Please enter the text filename  you want to create: ')
#..............................................................................................
file = open(pathout + os.sep + write_into_file, 'w')
file.write('current(A)' + '\t' + 'time(s)' + '\n')

for m in range(len(t_list)):
    file.write(str(i[m]) +'\t'+ str(t_list[m]) +'\n')  

file.close()
print('Creating text file completed')


"""
i_list = []
for n in range(len(i)):
    if n <= 66439:
        i_list.append(i[n])
    else:
        i_list.append((i[n]-50.0*10**(-8)))

i = np.asarray(i_list)
"""



