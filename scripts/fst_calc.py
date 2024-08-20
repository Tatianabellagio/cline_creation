import pandas as pd
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import allel

#inputs
parameter_space = snakemake.input['parameter_space'] 
vcf_file = snakemake.input['vcf_file'] 
population_counts_file = snakemake.input['population_counts_file'] 
causal_loci_file = snakemake.input['causal_loci_file'] 
hash = snakemake.params['hash']

#outputs
mean_fst_file = snakemake.output["mean_fst_file"]
fst_per_pos_file = snakemake.output["fst_per_pos"]
fst_manhattan_file = snakemake.output["fst_manhattan"]
fst_across_pops_file = snakemake.output["fst_across_pops"]
fst_heatmap_file = snakemake.output["fst_heatmap"]

causal_loci = pd.read_csv(causal_loci_file,header=None)[0].values


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

# Create the subpops list
subpops = []
current_index = 0

for pop, count in population_counts.items():
    # Create a list of indices for this population
    if count > 0:
        subpop_indices = list(range(current_index, current_index + count))
        subpops.append(subpop_indices)
    current_index += count

a, b, c = allel.weir_cockerham_fst(geno_array, subpops)

fst = np.sum(a) / (np.sum(a) + np.sum(b) + np.sum(c))


with open(mean_fst_file, 'w') as file:
    file.write(str(fst))


fst_per_pos = (np.sum(a, axis=1) /(np.sum(a, axis=1) + np.sum(b, axis=1) + np.sum(c, axis=1)))

fst_per_pos = pd.DataFrame({'pos':pos, 'fst':fst_per_pos})

fst_per_pos.to_csv(fst_per_pos_file, index=False)

df = fst_per_pos.copy()

# Set the color palette
colors = sns.color_palette("crest", n_colors=5)

# Initialize the plot
plt.figure(figsize=(20, 6))

# Plot the data
sns.scatterplot(
    data=df,
    x='pos',
    y='fst',
    edgecolor=None,
    s=40,
    alpha=0.5
)

# Aesthetics
plt.xlabel('Position')
plt.ylabel('Fst')
plt.title('Manhattan Plot of Fst Values')

# Add vertical lines at specific intervals
intervals = []
for i in causal_loci:
    intervals.append((i-1000, i+1000))

for interval in intervals:
    middle = (interval[0] + interval[1]) / 2
    plt.axvline(x=middle, color='grey', linestyle='dashed')

# Remove the top and right spines for a cleaner look
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Show the plot
plt.tight_layout()
plt.savefig(fst_manhattan_file)


g = allel.GenotypeArray(geno_array)

# Initialize a matrix to store the pairwise FST values
n_pops = len(subpops)
fst_matrix = np.zeros((n_pops, n_pops))

# Loop over all pairs of populations
for i in range(n_pops):
    for j in range(i + 1, n_pops):
        # Get the allele counts for the two populations
        ac1 = g.count_alleles(subpop=subpops[i])
        ac2 = g.count_alleles(subpop=subpops[j])

        # Calculate Weir & Cockerham's FST between the two populations
        num, den = allel.hudson_fst(ac1, ac2)
        fst_value = np.sum(num) / np.sum(den)
        
        # Store the FST value in the matrix
        fst_matrix[i, j] = fst_value
        fst_matrix[j, i] = fst_value  # Symmetric matrix

# Convert the matrix to a DataFrame for better readability
fst_df = pd.DataFrame(fst_matrix, columns=[f'pop{i}' for i in range(n_pops)], index=[f'pop{i}' for i in range(n_pops)])
fst_df.to_csv(fst_across_pops_file)

sns.heatmap(fst_df.iloc[1:, 1:], square=True)
plt.savefig(fst_heatmap_file)
