import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
# Specify the paths to the directories containing the ene*.dat files
data_directory1 = '/home/adurden/Ari/TAB/tab-model/data6'
data_directory2 = '/home/adurden/Ari/TAB/tab-model/data4/old'
# Define the column names

def process_data2(directory):
    columns2 = ["t", "S0 pop", "S1 pop", "S2 pop", "poptot"]
    file_paths = glob.glob(os.path.join(directory, 'dpop*.dat'))
    dfs = []        
    count=0
    for file_path in file_paths:
        df = pd.read_csv(file_path, delim_whitespace=True, names=columns2, skiprows=1)
        dfs.append(df)
        count+=1
    print(count)
    combined_df = pd.concat(dfs)
    average_df = combined_df.groupby('t').mean().reset_index()
    return average_df, len(file_paths)
avg_df1, file_count1 = process_data2(data_directory1)
avg_df2, file_count2 = process_data2(data_directory2)
#avg_df3, file_count3 = process_data2(data_directory3)
# Plot the data

columns_pop = ["t", "S0 pop", "S1 pop", "S2 pop", "poptot"]
pop_df = pd.read_csv('/home/adurden/Ari/TAB/tab-model/data4/pop.dat', delim_whitespace=True, names=columns_pop, index_col=False)
pop_df['t'] = pop_df['t'] / 100


plt.figure(figsize=(12, 6))
plt.plot(avg_df1['t'], avg_df1['S0 pop'], label='new')
plt.plot(avg_df2['t'], avg_df2['S0 pop'], label='old')
#plt.plot(avg_df3['t'], avg_df3['S0 pop'], label='rescales all directions')
plt.scatter(pop_df['t'], pop_df['S0 pop'], label='exact', color='black', marker='o', s=7)  # Plot as dots with smaller size
plt.xlabel('Time (fs)', fontsize=24)
plt.ylabel('S0 Population', fontsize=24)

#plot x axis up to 300 fs
plt.xlim(0, 300)

plt.title('S0 pop vs Time', fontsize=24)
plt.legend(fontsize=18)
plt.grid(True)
plt.savefig('momentum_rescale.pdf', format='pdf')
#plt.show()