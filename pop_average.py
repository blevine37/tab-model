import pandas as pd
import glob
import matplotlib.pyplot as plt

# Define the pattern for the .dat files
file_pattern = "pop*.dat"

# Get a list of all files matching the pattern
file_list = glob.glob(file_pattern)

# Initialize accumulators
sum_df = None
count_df = None

# Loop over the list of files and read/process each one
for file in file_list:
    df = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=1,
                     names=['t', 'state_0_pop', 'state_1_pop', 'poptot'])
    
    if sum_df is None:
        sum_df = df.copy()
        count_df = pd.DataFrame({'t': df['t'], 'count': 1})
    else:
        sum_df['state_0_pop'] += df['state_0_pop']
        sum_df['state_1_pop'] += df['state_1_pop']
        sum_df['poptot'] += df['poptot']
        count_df['count'] += 1

# Compute the averages
average_df = sum_df.copy()
average_df['state_0_pop'] /= count_df['count']
average_df['state_1_pop'] /= count_df['count']
average_df['poptot'] /= count_df['count']

# Write the average data to avg_pop.dat
average_df.to_csv('avg_apop.dat', sep='\t', index=False, header=['t', 'state_0_pop', 'state_1_pop', 'poptot'])

# Plot state 1 pop with t
plt.figure(figsize=(10, 6))
plt.plot(average_df['t'], average_df['state_1_pop'], label='dState 1 Pop')
plt.xlabel('Time (t)')
plt.ylabel('aState 1 Population')
plt.title('aState 1 Population vs Time')
plt.legend()
plt.grid(True)
plt.show()
