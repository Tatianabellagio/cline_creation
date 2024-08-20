
######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/Users/tatiana/mambaforge/envs/simulations/lib/python3.11/site-packages', '/Users/tatiana/Library/Caches/snakemake/snakemake/source-cache/runtime-cache/tmp0j85vu5q/file/Users/tatiana/Documents/grenenet/simulations/test_slim/cline_creation/scripts', '/Users/tatiana/Documents/grenenet/simulations/test_slim/cline_creation/scripts']); import pickle; snakemake = pickle.loads(b'\x80\x04\x95e\x08\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8c\x18parameters/a389d997.json\x94\x8c\x1dresults/vcfs/vcf_a389d997.vcf\x94\x8c8results/population_counts/population_counts_a389d997.pkl\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\x0fparameter_space\x94K\x00N\x86\x94\x8c\x08vcf_file\x94K\x01N\x86\x94\x8c\x16population_counts_file\x94K\x02N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x18\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x1e)}\x94\x8c\x05_name\x94h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bh\x10h\nh\x12h\x0bh\x14h\x0cub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c-results/alleles_df/allele_counts_a389d997.csv\x94\x8c+results/alleles_df/allele_freq_a389d997.csv\x94\x8c,results/lfmm/allele_counts_lfmm_a389d997.csv\x94\x8c*results/lfmm/env_vars/env_var_a389d997.csv\x94\x8c+results/lfmm/lest_pos/left_pos_a389d997.csv\x94e}\x94(h\x0e}\x94(\x8c\x10allele_counts_df\x94K\x00N\x86\x94\x8c\x0eallele_freq_df\x94K\x01N\x86\x94\x8c\x12allele_counts_lfmm\x94K\x02N\x86\x94\x8c\x0cenv_var_lfmm\x94K\x03N\x86\x94\x8c\rleft_pos_lfmm\x94K\x04N\x86\x94uh\x16]\x94(h\x18h\x19eh\x18h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bh3h,h5h-h7h.h9h/h;h0ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94\x8c\x08a389d997\x94a}\x94(h\x0e}\x94\x8c\x04hash\x94K\x00N\x86\x94sh\x16]\x94(h\x18h\x19eh\x18h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bhMhJub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94hJa}\x94(h\x0e}\x94\x8c\x04hash\x94K\x00N\x86\x94sh\x16]\x94(h\x18h\x19eh\x18h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bhMhJub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c0/var/folders/89/m0n7cpqn6153r2j98t2n5yq40000gr/T\x94M\x00xMqre}\x94(h\x0e}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94\x8c\x06mem_mb\x94K\x03N\x86\x94\x8c\x07mem_mib\x94K\x04N\x86\x94uh\x16]\x94(h\x18h\x19eh\x18h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bhqK\x01hsK\x01huhnhwM\x00xhyMqrub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0e}\x94h\x16]\x94(h\x18h\x19eh\x18h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x18sNt\x94bh\x19h\x1ch\x1e\x85\x94R\x94(h\x1e)}\x94h"h\x19sNt\x94bub\x8c\x06config\x94}\x94\x8c\nparameters\x94]\x94(}\x94(\x8c\x0brepr_scheme\x94K\x01\x8c\x07mutrate\x94\x8c\x041e-8\x94\x8c\x0fcapacity_charge\x94G?\xb9\x99\x99\x99\x99\x99\x9a\x8c\x10initial_pop_size\x94M\xe8\x03\x8c\rsubpop_number\x94K\x05\x8c\x0emigration_rate\x94G?\x84z\xe1G\xae\x14{\x8c\x08BURNIN_1\x94M\xe8\x03\x8c\x08BURNIN_2\x94M\xd0\x07\x8c\x08last_gen\x94M\x10\'\x8c\x12final_sel_strength\x94G?\xb9\x99\x99\x99\x99\x99\x9a\x8c\x0bactive_loci\x94K2\x8c\x02h2\x94G?\xe0\x00\x00\x00\x00\x00\x00u}\x94(\x8c\x0brepr_scheme\x94K\x02\x8c\x07mutrate\x94\x8c\x041e-9\x94\x8c\x0fcapacity_charge\x94G?\xc9\x99\x99\x99\x99\x99\x9a\x8c\x10initial_pop_size\x94M\xd0\x07\x8c\rsubpop_number\x94K\n\x8c\x0emigration_rate\x94G?\xa9\x99\x99\x99\x99\x99\x9a\x8c\x08BURNIN_1\x94M\xd0\x07\x8c\x08BURNIN_2\x94M\xa0\x0f\x8c\x08last_gen\x94M N\x8c\x12final_sel_strength\x94G?\xc9\x99\x99\x99\x99\x99\x9a\x8c\x0bactive_loci\x94Kd\x8c\x02h2\x94G?\xe3333333ues\x8c\x04rule\x94\x8c\x0egen_alleles_df\x94\x8c\x0fbench_iteration\x94K\x00\x8c\tscriptdir\x94\x8cN/Users/tatiana/Documents/grenenet/simulations/test_slim/cline_creation/scripts\x94ub.'); from snakemake.logging import logger; logger.printshellcmds = False; __real_file__ = __file__; __file__ = '/Users/tatiana/Documents/grenenet/simulations/test_slim/cline_creation/scripts/gen_alleles_df.py';
######## snakemake preamble end #########
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
all_alt_allele_count.to_csv(allele_counts_lfmm)


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