import pyslim 
import pandas as pd
import msprime
import tskit

import allel
import numpy as np

import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns

import sys

import pickle

# Access the runid argument
runid = sys.argv[1]
migration_rate = sys.argv[2]

mutrate = 1e-6
causal_pos = 15005
initial_ne = 1000
# Load the .trees file
output_slim_tree = tskit.load('results/'+ runid + "_tree_final.trees")

pop_names = []
for pops in output_slim_tree.tables.populations:
    if pops.metadata is not None:
        # Decode the metadata (assuming it's stored as JSON by SLiM)
        pop_names.append(pops.metadata['name'])

# Create the demography from the tree sequence
demography = msprime.Demography.from_tree_sequence(output_slim_tree)

# Set the initial size for pop_0
for pop in demography.populations:
    if pop.name == "pop_0":
        pop.initial_size = initial_ne

# Add migration rates using population indices
subpop_number = len(pop_names)
for i in range(1, subpop_number):
    # Set migration rate from subpop i to subpop i+1
    demography.set_migration_rate(source=i, dest=i + 1, rate=migration_rate)
    
    # Set migration rate from subpop i+1 to subpop i (bi-directional)
    demography.set_migration_rate(source=i + 1, dest=i, rate=migration_rate)

rts = pyslim.recapitate(output_slim_tree,
            recombination_rate=3e-6,
            demography=demography,
            #ancestral_Ne=200, 
            random_seed=5)

next_id = pyslim.next_slim_mutation_id(rts)
ts = msprime.sim_mutations(
           rts,
           rate=mutrate,
           model=msprime.SLiMMutationModel(type=0, next_id=next_id),
           keep=True,
)

def generate_nucleotides_until_diverse(ts, causal_pos):
    all_same = True
    while all_same:
        print(all_same)
        # Generate nucleotides
        nts = pyslim.generate_nucleotides(ts)
        nts = pyslim.convert_alleles(nts)
        
        # Get the sites table
        sites_table = nts.tables.sites

        # Find the index of the site at the causal position
        index_causal_pos = None
        for index, site in enumerate(sites_table):
            if site.position == causal_pos:
                index_causal_pos = index
                break

        if index_causal_pos is None:
            raise ValueError(f"Site at position {causal_pos} not found.")

        # Get the derived states at the causal position
        derived_states = []
        mutations_table = nts.tables.mutations
        for mut in mutations_table:
            if mut.site == index_causal_pos:
                derived_states.append(mut.derived_state)

        # Check if all derived states are the same
        all_same = (pd.Series(derived_states) == derived_states[0]).all()

    # Return the nts once we have diverse derived states
    return nts, derived_states

# Call the function with your time series and causal position
nts, derived_states = generate_nucleotides_until_diverse(ts, causal_pos)


sites_table = nts.tables.sites

# Print the ancestral state of all mutations
for index, site in enumerate(sites_table):
    if site.position == causal_pos:
        index_causal_pos = index
        print(f"Position: {site.position}, Ancestral State: {site.ancestral_state}")

## now there is an ancestral state 
mutations_table = nts.tables.mutations

# Print the ancestral state of all mutations
for index, mut in enumerate(mutations_table):
    if mut.site == index_causal_pos:
        #print(mut)
        #print(mut)
        #print(index)
        pick_as_ancestral = mut.derived_state
        ##print(f"Position: {site.position}, Ancestral State: {site.ancestral_state}")

## assign the ancestral sites to one of the 2 mutations 0 or 1 
sites_table = nts.tables.sites
# Create a new site table to hold the modified entries
new_sites_table = tskit.SiteTable()


# Iterate through the existing sites, modify the ancestral state, and add them to the new table
for site in sites_table:
    # Modify the ancestral state of the causal mutation '
    if site.position == causal_pos:
        print(site)
        new_sites_table.add_row(
            position=site.position,
            ancestral_state= pick_as_ancestral, 
            metadata=site.metadata
        )  
    else:  
        new_sites_table.add_row(
            position=site.position,
            ancestral_state= site.ancestral_state, 
            metadata=site.metadata
        )

tables = nts.dump_tables()
tables.sites.replace_with(new_sites_table)
# Now you can save the modified tree sequence
nts = tables.tree_sequence()


# Identify multiallelic sites
multiallelic_sites = set()
invariant_sites = set()

for site in nts.sites():
    alleles = set()
    for mutation in site.mutations:
        alleles.add(mutation.derived_state)

    if len(alleles) > 2:  # Multiallelic site
        #print(site)
        multiallelic_sites.add(site.id)

    elif len(alleles) == 0:
        print('invariant site')
        invariant_sites.add(site.id)

# Filter out the multiallelic sites
def site_is_biallelic(site):
    return site.id not in multiallelic_sites

filtered_ts = nts.delete_sites(
    [site.id for site in nts.sites() if site.id in multiallelic_sites]
)




# Write the filtered tree sequence to a VCF file
name_vcf = f"results/vcfs/vcf_w_neutral_mutid_{runid}.vcf"
with open(name_vcf, 'w') as file:
    # Pass the file object as the output parameter
    filtered_ts.write_vcf(output=file) #allow_position_zero = True)


population_counts = {pop.id: 0 for pop in nts.populations()}

for ind in nts.individuals():
    # Each individual has a population
    pop_id = ind.population
    population_counts[pop_id] += 1

# Save the dictionary as a pickle file
with open(f'results/pop_counts/population_counts_{runid}.pkl', 'wb') as file:
    pickle.dump(population_counts, file)

total_pop_size = sum(population_counts.values())

vcf_file = f'results/vcfs/vcf_w_neutral_mutid_{runid}.vcf'
vcf = allel.read_vcf(vcf_file)
samples = vcf['samples']
geno_array = vcf['calldata/GT']
pos = vcf['variants/POS']
chrom = vcf['variants/CHROM']

## delete all the sites with thwo alternative loci 

mask = np.any(geno_array == 2, axis=(1, 2))
# Invert the mask to select the entries that do NOT have a '2'
inverted_mask = ~mask
# Use the mask to filter out the entries with at least one '2'
geno_array = geno_array[inverted_mask]
pos = pos[inverted_mask]

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

all_alt_allele_count.to_csv(f'results/alleles_df/all_alt_allele_count_{runid}.csv')
all_alt_allele_count.to_csv(f'results/alleles_df/all_alt_allele_freq_{runid}.csv')

total_genomes = geno_array.shape[1]*2

total_genomes
min_freq = 0.01
min_count = total_genomes * min_freq

print(all_alt_allele_count.shape)
all_alt_allele_freq = all_alt_allele_freq[all_alt_allele_count.sum(axis=1) > min_count]

#sample_counts = all_alt_allele_freq.sample(100)
sample_counts = all_alt_allele_freq.copy()

sample_counts = sample_counts.reset_index() #.reset_index()
# Melt the DataFrame to long format
long_df = sample_counts.melt(id_vars=['index'], var_name='Population', value_name='Allele_Count')

# Rename 'index' to 'SNP'
long_df.rename(columns={'index': 'SNP'}, inplace=True)

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


results_df.to_csv(f'results/snps_w_clines/snps_w_clines_{runid}.csv')

# Plot the dynamics of alleles for each SNP across populations
plt.figure(figsize=(10, 6))
sns.lineplot(data=long_df, x='Population', y='Allele_Count', hue='SNP', palette='gray', alpha = 0.2)

# Overlay significant SNPs
for _, row in results_df.iterrows():
    if row['Color'] in ['green', 'red']:
        if row['SNP'] != 15005:
            snp_data = long_df[long_df['SNP'] == row['SNP']]
            line = sns.lineplot(data=snp_data, x='Population', y='Allele_Count', color=row['Color'], alpha=0.2)
            sns.scatterplot(data=snp_data, x='Population', y='Allele_Count',  s = 50, color=row['Color'], edgecolor=row['Color'], linewidth = 0, alpha=0.2)

        if row['SNP'] == 15005:
            snp_data = long_df[long_df['SNP'] == row['SNP']]
            line = sns.lineplot(data=snp_data, x='Population', y='Allele_Count', color='Black')
            sns.scatterplot(data=snp_data, x='Population', y='Allele_Count',  s = 50, color='Black', edgecolor='Black', linewidth = 0)

plt.title('Dynamics of Alleles in Each SNP Across Populations')
plt.xlabel('Population')
plt.ylabel('Allele Count')
plt.ylim(-0.2, 1.1)
plt.legend(title='SNP', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(rotation=90)  # Rotate x-axis labels if needed
plt.tight_layout()
plt.savefig(f'results/plots_af/allele_freq_runid{runid}.png')
plt.show()