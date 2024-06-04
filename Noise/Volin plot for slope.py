# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 14:29:14 2024

@author: sndkp
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#%matplotlib inline

# Specify the Excel file path
excel_file = r'F:\Personal Data\NbO2\Time sereis for voltage driven noise till april 2024/All slopes in one column.xlsx'

# Read the data from the single sheet
df = pd.read_excel(excel_file)

# Fill NaN values with the mean of each column
df_cleaned = df.fillna(df.mean())

# Create a bean plot for S1, S2, and S3 with custom labels
plt.figure(figsize=(8, 6))

# Customize the colors
sns.violinplot(data=df_cleaned[['S206A', 'S206B', 'S206C']], palette=['skyblue', 'salmon', 'lightgreen'])

# Set custom labels for the x-axis
plt.xticks(ticks=[0, 1, 2], labels=['Before transitions', 'Around Transitions', 'After transitions'])

plt.title(r"Value of $\alpha$ from 5 samples")
#plt.xlabel("Columns")
plt.ylabel(r"$\alpha$")


plt.show()
