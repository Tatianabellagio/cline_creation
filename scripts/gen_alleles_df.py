import pandas as pd
import pickle
import allel 
import numpy as np

#inputs
parameter_space = snakemake.input['parameter_space'] 
vcf_file = snakemake.input['vcf_file'] 
population_counts_file = snakemake.input['population_counts_file'] 
hash = snakemake.params['hash']



#outputs
allele_counts_df_file = snakemake.output["allele_counts_df"]
allele_freq_df_file = snakemake.output["allele_freq_df"]
allele_counts_lfmm = snakemake.output["allele_counts_lfmm"]
env_var_lfmm = snakemake.output["env_var_lfmm"]
left_pos_lfmm = snakemake.output["left_pos_lfmm"]

vcf = allel.read_vcf(vcf_file)
samples = vcf['samples']
geno_array = vcf['calldata/GT']
pos = vcf['variants/POS']
chrom = vcf['variants/CHROM']

mask = np.any(geno_array == 2, axis=(1, 2))
# Invert the mask to select the entries that do NOT have a '2'
inverted_mask = ~mask
# Use the mask to filter out the entries with at least one '2'
geno_array = geno_array[inverted_mask]
pos = pos[inverted_mask]


with open(population_counts_file, 'rb') as f:
    population_counts = pickle.load(f)

total_pop_size = sum(population_counts.values())

# Check if the geno_array has enough individuals
assert geno_array.shape[1] >= total_pop_size, "geno_array does not have enough individuals."

# Initialize variables to keep track of slicing
subpops = {}
current_position = 0

# Slice the geno_array into subpopulations based on population_counts
for pop_id, pop_size in population_counts.items():
    if pop_size > 0:  # Only consider populations with non-zero size
        start = current_position
        end = start + pop_size
        subpops[f'pop{pop_id}'] = geno_array[:, start:end, :]
        current_position = end


all_alt_allele_count = {}
all_alt_allele_freq = {}
col_names = []
col_todel = []
#for order, pop in enumerate([pop1, pop2, pop3, pop3]):
for pop_name, pop in subpops.items():
    if pop.shape[1] > 0 :
        total_alleles = pop.shape[1] * 2 
        alt_count = pop.sum(axis=2).sum(axis=1)
        alt_freq = alt_count / total_alleles
        #alt_freq = alt_freq.round(4)
        #chrom_pos = pd.Series(chrom.astype(str)) + '_' +  pd.Series(pos.astype(str))
        col_names.append('chrom_pos' + pop_name)
        col_todel.append('chrom_pos' + pop_name)
        col_names.append(pop_name)
        alt_allele_count = pd.DataFrame(data = {'chrom_pos': pos, pop_name: alt_count})
        all_alt_allele_count[pop_name] = alt_allele_count
        alt_allele_freq = pd.DataFrame(data = {'chrom_pos': pos, pop_name: alt_freq})    
        all_alt_allele_freq[pop_name] = alt_allele_freq


all_alt_allele_freq = pd.concat(all_alt_allele_freq,axis=1)
all_alt_allele_count = pd.concat(all_alt_allele_count,axis=1)

all_alt_allele_freq.columns = col_names 
all_alt_allele_count.columns = col_names 

og_positions = all_alt_allele_count[col_todel[0]]
og_positions_dict = og_positions.to_dict()

all_alt_allele_freq = all_alt_allele_freq.drop(col_todel,axis=1)
all_alt_allele_count = all_alt_allele_count.drop(col_todel,axis=1)

all_alt_allele_count.index = all_alt_allele_count.index.map(og_positions_dict)
all_alt_allele_freq.index = all_alt_allele_freq.index.map(og_positions_dict)

all_alt_allele_count.to_csv(allele_counts_df_file)
all_alt_allele_count.to_csv(allele_freq_df_file)


## for lfmmm
total_genomes = geno_array.shape[1]*2
min_freq = 0.05
min_count = total_genomes * min_freq

print(all_alt_allele_count.shape)
all_alt_allele_count = all_alt_allele_count[all_alt_allele_count.sum(axis=1) > min_count]
all_alt_allele_count.to_csv(allele_counts_lfmm, index=None)

# read teh sequence of envrionemnts 
env_var = pd.read_csv('env_var.txt', header=None)[0]
# Create the dictionary
pop_env = {f'pop{index+1}': value for index, value in enumerate(env_var)}

env_var = pd.Series(list(all_alt_allele_count.columns.map(pop_env)))

env_var.to_csv(env_var_lfmm,index=None, header=None)

left_pos = all_alt_allele_count.index
left_pos = pd.Series(left_pos).reset_index()
left_pos.columns = ['ignore', 'chrom_pos_left']

left_pos.to_csv(left_pos_lfmm,index=None)