import pandas as pd
import json 

#inputs
parameter_space = snakemake.input['parameter_space'] 
allele_freq_gen0 = snakemake.input['allele_freq_gen0'] 
allele_freq_df_file = snakemake.input['allele_freq_df_file'] 
population_counts_file = snakemake.input['population_counts_file']
allele_counts_df_file = snakemake.input['allele_counts_df_file'] 

#outputs
delta_allele_freq_lfmm = snakemake.output["delta_allele_freq_lfmm"]
delta_p_env_var_lfmm = snakemake.output["delta_p_env_var_lfmm"]

# Open and read the JSON file
with open(parameter_space, 'r') as f:
    parameters = json.load(f)

common_garden_number = parameters['common_garden_number']

print(allele_freq_df_file)

gen0 = pd.read_csv(allele_freq_gen0).set_index('Unnamed: 0')

## gen 1
gen1 = pd.read_csv(allele_freq_df_file).set_index('Unnamed: 0')

## only use from gen 0 the populations that survived to gen 1
gen0 = gen0[gen1.columns]

#gen0.columns = gen1.columns
gen1_allele_counts = pd.read_csv(allele_counts_df_file).set_index('Unnamed: 0')

with open(population_counts_file, 'rb') as f:
    population_countsgen1 = pickle.load(f)
total_genomesgen1 = sum(population_countsgen1.values()) * 2

min_freq = 0.05
min_countgen1 = total_genomesgen1 * min_freq

delta_pgen1 = gen1 - gen0
delta_pgen1 = delta_pgen1[gen1_allele_counts.sum(axis=1) > min_countgen1]
delta_pgen1.to_csv(delta_allele_freq_lfmm, index=None)

env_var = pd.read_csv(f'env_var_acg_{common_garden_number}.txt', header=None)[0]
range_start = len(env_var) - 1
range_ends = range_start + len(env_var)
pop_initial = [f'pop{i}' for i in range(range_start, range_ends)]
pop_env = {key: value for key, value in zip(pop_initial, env_var)}
env_var = pd.Series(list(delta_pgen1.columns.map(pop_env)))

print(delta_pgen1.shape)
print(env_var.shape)
env_var.to_csv(delta_p_env_var_lfmm,index=None, header=None)
