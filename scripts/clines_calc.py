import pandas as pd
import pickle
import numpy as np
import statsmodels.api as sm


#inputs
parameter_space = snakemake.input['parameter_space'] 
allele_counts_df_file = snakemake.input['allele_counts_df_file']
allele_freq_df_file = snakemake.input['allele_freq_df_file']  
population_counts_file = snakemake.input['population_counts_file'] 
hash = snakemake.params['hash']


#outputs
cline_file = snakemake.output["cline_file"]

with open(population_counts_file, 'rb') as f:
    population_counts = pickle.load(f)

total_pop_size = sum(population_counts.values())

all_alt_allele_count = pd.read_csv(allele_counts_df_file).set_index('Unnamed: 0')
all_alt_allele_freq = pd.read_csv(allele_freq_df_file).set_index('Unnamed: 0')

total_genomes = total_pop_size * 2
min_freq = 0.05
min_count = total_genomes * min_freq

all_alt_allele_freq = all_alt_allele_freq[all_alt_allele_count.sum(axis=1) > min_count]

#sample_counts = all_alt_allele_freq.sample(100)
sample_counts = all_alt_allele_freq.copy()

sample_counts = sample_counts.reset_index() #.reset_index()
# Melt the DataFrame to long format
long_df = sample_counts.melt(id_vars=['Unnamed: 0'], var_name='Population', value_name='Allele_Count')

# Rename 'index' to 'SNP'
long_df.rename(columns={'Unnamed: 0': 'SNP'}, inplace=True)

results = []

# Run linear models for each SNP
for snp in long_df['SNP'].unique():
    snp_data = long_df[long_df['SNP'] == snp]
    X = sm.add_constant(snp_data['Population'].str.extract('(\d+)$').astype(int))  # Extract numeric part of population names
    y = snp_data['Allele_Count']
    model = sm.OLS(y, X).fit()
    slope = model.params[0]
    p_value = model.pvalues[0]
    results.append((snp, slope, p_value))

# Convert results to DataFrame
results_df = pd.DataFrame(results, columns=['SNP', 'Slope', 'P_Value'])

# Determine significance and assign colors
alpha = 0.05
results_df['Color'] = np.where((results_df['P_Value'] < alpha) & (results_df['Slope'] > 0), 'green', 
                               np.where((results_df['P_Value'] < alpha) & (results_df['Slope'] < 0), 'red', 'blue'))


cline = results_df[results_df['Color'].isin(['red', 'green'])]

cline.to_csv(cline_file)