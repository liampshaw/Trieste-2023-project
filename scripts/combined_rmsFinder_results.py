import pandas as pd
import glob
import os
import re

# Directory path
directories = glob.glob('genomes/n100/*')

dataframes = []
for directory in directories:
    # List to store DataFrames
    genus = directory.split('/')[-1]
    print(genus)
    # Loop through each file in the directory that ends with 'RMS.csv'
    for filename in glob.glob(os.path.join(directory, '*RMS.csv')):
        # Extract the 'GCF' part from the filename
        genome_id = '_'.join(filename.split('/')[-1].split('_')[0:2])
        # Read the CSV file into a DataFrame
        if os.path.getsize(filename)>1:
            df = pd.read_csv(filename)
            # Check if the DataFrame is not empty
            if not df.empty:
                # Add the 'genome' column
                df['genome'] = genome_id
                df['genus'] = genus
                # Append the DataFrame to the list
                dataframes.append(df)

# Concatenate all DataFrames into a single DataFrame
final_df = pd.concat(dataframes)

# Display the final DataFrame (optional)
final_df.to_csv('combined.csv', index=False)


# Get number of genomes to normalise by
total_genomes = {re.sub('.*/', '', d):len(glob.glob(os.path.join(d, '*.gz'))) for d in directories}
print(total_genomes)


# Summarise by genus/sequence
summary_df = final_df.groupby(['genus', 'sequence']).agg({'genome': pd.Series.nunique})
# Reset the index 
summary_df.reset_index(inplace=True)

print(summary_df)
# Add genus counts
summary_df['total_genomes'] = summary_df['genus'].map(total_genomes)
summary_df['normalized_count'] = summary_df['genome'] / summary_df['total_genomes'] 
summary_df.to_csv('target_db.csv', index=False)