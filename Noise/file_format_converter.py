"""Created on Thu May 17 14:15:43 2018 @author: Dasharath Adhikari"""

import os
import sys


def itx_to_txt_converter(raw_data_path, text_data_path):
    """
    raw_data_path is the location fo raw data from the experiment
    text_data_path is the file location for the output text data
    This function converts igor file (.itx file) form noise spectroscopy
    measurement data file (Igor files) into text file.
    """
    # creates list from the sorted files in the raw_data_path directory
    filelist = sorted(os.listdir(raw_data_path))
    linecount = int(input('please enter the first file number >> '))
    for file in filelist:
        if file.endswith('.itx'):
            # sys.exit(-1)
            print('-'*20)
            print('File number: ', linecount)
            print('-'*20)
            with open(raw_data_path + os.sep + file, 'r',
                      encoding='ISO-8859-1') as myfile:
                readlines = myfile.read()
                split_lines = readlines.split('\n')
                # create text file in the text_data_path folder
                # file name will be VT followed by line count value
                FileToWriteInto = open(text_data_path + os.sep +
                                       'VT'+ str(linecount).zfill(3)+'.txt',
                                       'w')
                FileToWriteInto.write('time_s' + '\t' + 'Vx_V' + '\t' +
                                      'v_y (V)' + '\n')
                for line in range(len(split_lines)):
                    if((line < 3) or (line > len(split_lines) - 9)):
                        continue
                        # print("Comment lines are encountered")
                    else:
                        each_data = split_lines[line].split('\t')
                        t_s = each_data[0]
                        Vx_V = each_data[1]
                        Vy_V = each_data[2]
                        final_string = t_s + '\t' + Vx_V + '\t' + Vy_V
                        FileToWriteInto.write(final_string + '\n')
                FileToWriteInto.close()
            linecount += 1
