import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.ticker import AutoMinorLocator
import matplotlib.ticker as ticker

# Matplotlib configuration
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.weight'] = 'bold'
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['axes.titlesize'] = 16

# Folder path
folder_path = r"F:\Personal Data\NbO2\Size dependent study\Line 2\D4"

# List all CSV files
csv_files = [file for file in os.listdir(folder_path) if file.startswith('') and file.endswith('.csv')]

# Function to plot a batch of files
def plot_batch(file_batch, batch_num):
    plt.figure()
    
    for file_num, csv_file in enumerate(file_batch, start=1):
        file_path = os.path.join(folder_path, csv_file)
        try:
            df = pd.read_csv(file_path, skiprows=8)
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            continue
        
        # Select only the first half of the data
        half_index = len(df) // 2
        df = df.iloc[:half_index]

        current = df['Reading'].values
        voltage = df['Value'].values
        current_np = np.array(current)
        voltage_np = np.array(voltage)

        plt.plot(voltage_np, current_np * 1000, label=f'Run {csv_file[:-4]}', linestyle='-', marker='.', linewidth=1)

    plt.ylabel('Current (mA)')
    plt.xlabel('Voltage (V)')
    plt.legend()
    plt.tick_params(axis="y", direction="in", length=6, width=2)
    plt.tick_params(axis="x", direction="in", length=6, width=2)

    minor_locator_x = ticker.AutoMinorLocator(4)
    minor_locator_y = ticker.AutoMinorLocator(4)
    plt.gca().xaxis.set_minor_locator(minor_locator_x)
    plt.gca().yaxis.set_minor_locator(minor_locator_y)
    plt.gca().tick_params(axis='both', which='minor', direction='in', length=6, width=2)

    plot_file_name = f'IV_Plot_Batch_{batch_num}.png'
    plot_file_path = os.path.join(folder_path, plot_file_name)
    plt.savefig(plot_file_path, bbox_inches='tight', dpi=300)
    plt.close()

# Process files in batches of three
batch_size = 3
for i in range(0, len(csv_files), batch_size):
    file_batch = csv_files[i:i+batch_size]
    plot_batch(file_batch, i // batch_size + 1)

print("Plots saved successfully.")
