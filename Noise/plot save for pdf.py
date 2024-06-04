# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 19:24:50 2024

@author: sndkp
"""
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Folder containing the CSV files
folder_path = r"F:\APS march meeting 2024\Data\325 current driven noise\Noise\16 Min"

# Loop through all files in the folder
for filename in os.listdir(folder_path):
    if filename.startswith("it_") and filename.endswith(".csv"):
        # Load the CSV file into a DataFrame
        file_path = os.path.join(folder_path, filename)
        df = pd.read_csv(file_path, skiprows=0, encoding_errors='ignore', on_bad_lines='skip')

        # Calculate the mean of the "Reading" column
        mean_reading = df["Reading"].iloc[1:].mean()

        # Subtract the mean from the "Reading" column
        adjusted_reading = df["Reading"] - mean_reading

        # Plot the KDE for the adjusted "Reading" column
        sns.kdeplot(adjusted_reading, label=f"Value: {df['Value'][0]}")
                    
                # Add legend outside the loop to avoid multiple legends
        #plt.legend()
                        
                # Set labels and title
        plt.xlabel("Adjusted Reading")
        plt.ylabel("Density")
        plt.title("Kernel Density Estimate of Adjusted Reading")
                        
                # Save the plot in the same folder with a unique identifier in the filename
        plot_filename = f"kde_plot_with_legend_{filename[:-4]}.png"
        plt.savefig(os.path.join(folder_path, plot_filename))
                        
                # Show the plot
        plt.show()
        plt.close()
