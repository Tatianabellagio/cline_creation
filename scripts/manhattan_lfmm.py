import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

#inputs
lfmm_pvalues_file = snakemake.input['lfmm_pvalues'] 
left_pos_lfmm_file = snakemake.input['left_pos_lfmm'] 
causal_loci_file = snakemake.input['causal_loci_file'] 

#outputs
lfmm_manhattan = snakemake.output["lfmm_manhattan"]
sign_pos_df = snakemake.output["sign_pos_df"]

causal_loci = pd.read_csv(causal_loci_file,header=None)[0].values

left_pos = pd.read_csv(left_pos_lfmm_file)
p_values = pd.read_csv(lfmm_pvalues_file).reset_index().drop( 'Unnamed: 0', axis=1)

p_values = p_values.merge(left_pos, left_on = 'index', right_on = 'ignore', how = 'inner').drop( ['index', 'ignore'], axis=1)

p_values.columns = ['pvalue', 'pos']

threshold_value = 0.05 / len(p_values)

p_values['chrom'] =1 

sign_pos = p_values[p_values['pvalue'] < threshold_value]
sign_pos.to_csv(sign_pos_df,index=None)

## here someitmes the p values might be 0 so i will replace it by the lowest pvalue or he treshold for visual purposes 
if (p_values['pvalue'] == 0).any():

    lowerst_pvalue = p_values['pvalue'].min()   

    if lowerst_pvalue < threshold_value:
        lowest_value = lowerst_pvalue
    elif lowerst_pvalue > threshold_value:
        lowest_value = threshold_value

    p_values[p_values['pvalue'] == 0, 'pvalue'] = lowest_value

df = p_values.copy()

colors = sns.color_palette("crest", n_colors = 5)

# Parsing chromosome number and position
df['chromosome'] = df['chrom']
df['position'] = df['pos']
df['-log10(pvalue)'] = -np.log10(df['pvalue'])

# Calculate the offset for each chromosome to prevent overlap
chromosome_offsets = {}
offset = 0
for chrom in sorted(df['chromosome'].unique()):
    chromosome_offsets[chrom] = offset
    max_position = df[df['chromosome'] == chrom]['position'].max()
    offset += max_position + 1000000  # Adding 1 million as a buffer between chromosomes

# Apply offsets to positions
df['adjusted_position'] = df.apply(lambda row: row['position'] + chromosome_offsets[row['chromosome']], axis=1)

# Creating the Manhattan plot
plt.figure(figsize=(20, 6))

for chrom in sorted(df['chromosome'].unique()):
    subset = df[df['chromosome'] == chrom]
    plt.scatter(subset['adjusted_position'], subset['-log10(pvalue)'], c=colors[chrom % len(colors)], label=f'Chr {chrom}', s=10)


# Highlight clumped SNPs

# Aesthetics
plt.xlabel('Adjusted Position')
plt.ylabel('-log10(pvalue)')
#plt.title('Manhattan Plot')
#plt.grid(axis='y')
#plt.legend(title="Chromosome", bbox_to_anchor=(1.05, 1), loc='upper left')
ax = plt.gca()  # Get current axes
ax.spines['top'].set_visible(False)  # Remove the top spine
ax.spines['right'].set_visible(False)
# Threshold line (optional)
threshold = -np.log10(threshold_value)
plt.axhline(y=threshold, color='grey', linestyle='dashed')
#plt.title(f'{biovar}')  # Set the title

significant_snps = df[df['-log10(pvalue)'] > threshold]

# Highlight significant SNPs with red circles
plt.scatter(significant_snps['adjusted_position'], significant_snps['-log10(pvalue)'],
            facecolors='none', edgecolors='red', s=100, linewidths=2, label='Significant SNPs')

intervals = []
for i in causal_loci:
    intervals.append((i-1000, i+1000))
    
for interval in intervals:
    middle = (interval[0] + interval[1]) / 2
    plt.axvline(x=middle, color='grey', linestyle='dashed')

# Show the plot
plt.tight_layout()
plt.savefig(lfmm_manhattan)