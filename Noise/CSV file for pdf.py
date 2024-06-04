# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 20:11:05 2024

@author: sndkp
"""

import os
import pandas as pd

# Folder containing the CSV files
folder_path = r"F:\APS march meeting 2024\Data\325 current driven noise\Noise\Raw data"




# Create an empty list to store adjusted DataFrames
adjusted_dfs = []

# Loop through all files in the folder
for filename in os.listdir(folder_path):
    if filename.startswith("it_") and filename.endswith(".csv"):
        # Load the CSV file into a DataFrame
        file_path = os.path.join(folder_path, filename)
        try:
            df = pd.read_csv(file_path, skiprows=8, encoding_errors='ignore', on_bad_lines='skip')

            # Calculate the mean of the "Reading" column
            mean_reading = df["Reading"].mean()

            # Subtract the mean from the "Reading" column
            df["Adjusted_Reading"] = df["Reading"] - mean_reading

            # Add the adjusted DataFrame to the list
            adjusted_dfs.append(df)
        except pd.errors.ParserError as e:
            print(f"Error occurred while processing file: {filename}")
            print(f"Error details: {e}")
            continue  # Skip to the next file if an error occurs

# Concatenate all DataFrames into a single DataFrame
adjusted_readings_df = pd.concat(adjusted_dfs, ignore_index=True)

# Save the adjusted readings to a CSV file
adjusted_readings_df.to_csv(os.path.join(folder_path, "adjusted_readings.csv"), index=False)