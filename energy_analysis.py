import pandas as pd
import glob
import matplotlib.pyplot as plt

# Define the pattern for the .dat files
file_pattern = "ene*.dat"

# Get a list of all files matching the pattern
file_list = glob.glob(file_pattern)

# Initialize accumulators
sum_df = None
count_df = None

# Loop over the list of files and read/process each one
for file in file_list:
    df = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=1,
                     names=['t', 'Mean-Field Energy', 'PE Difference', 'Total Energy', 'norm','Entropy'])
    
    if sum_df is None:
        sum_df = df.copy()
        count_df = pd.DataFrame({'t': df['t'], 'count': 1})
    else:
        sum_df['Entropy'] += df['Entropy']
        count_df['count'] += 1

# Compute the averages
average_df = sum_df.copy()
average_df['Entropy'] /= count_df['count']

# Select the required columns
average_df_selected = average_df[['t', 'Entropy']]

# Write the average data to tot_energy.dat
average_df_selected.to_csv('avg_entropy.dat', sep='\t', index=False, header=['t', 'Entropy'])

# Plot state 1 pop with t
plt.figure(figsize=(10, 6))
plt.plot(average_df_selected['t'], average_df_selected['Entropy'], label='Average Entropy')
plt.xlabel('Time (t)')
plt.ylabel('Average Energy')
plt.title('Average Entropy vs Time (Adiabatic)')
plt.legend()
plt.grid(True)
plt.show()
