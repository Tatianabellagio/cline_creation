import pyslim 
import pandas as pd
import msprime
import tskit
import random
import json
import pickle

pos_nucl = [0, 1, 2, 3]

#inputs
parameter_space = snakemake.input['parameter_space'] 
tree_seq_file = snakemake.input['tree_seq_file'] 
causal_loci_file = snakemake.input['causal_loci_file'] 

hash = snakemake.params['hash']

#outputs
vcf_file = snakemake.output["vcf_files"]
population_counts_file = snakemake.output["population_counts_files"]
genetic_diversity_file = snakemake.output["genetic_diversity_files"]

print(vcf_file)

# read the causal loci 
causal_loci = pd.read_csv(causal_loci_file,header=None)[0].values

# Open and read the JSON file
with open(parameter_space, 'r') as f:
    parameters = json.load(f)

initial_ne = parameters['initial_ne']

## if the tree is after the common garden then there is no migration
if 'acg' in tree_seq_file:
    migration_rate = 0
elif 'bcg' in tree_seq_file:
    migration_rate = parameters['migration_rate']

print('migration rate is:', migration_rate)

mutrate = parameters['mutrate']
output_slim_tree = tskit.load(tree_seq_file)

## create demography 
pop_names = []
for pops in output_slim_tree.tables.populations:
    if pops.metadata is not None:
        # Decode the metadata (assuming it's stored as JSON by SLiM)
        pop_names.append(pops.metadata['name'])

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


## convert the causal loci into nucleotide type so they dont randomly get assing thte same nucleotide
print(tree_seq_file)
mutation_tables_mod = ts.tables.mutations
sites_table = ts.tables.sites
print(causal_loci)
for causal_locus in causal_loci:
    print(causal_locus)
    # Print the ancestral state of all mutations
    for index, site in enumerate(sites_table):
        if site.position == causal_locus:
            index_causal_locus = index

    replacement_dic = {}
    derived_states = [] 
    for mut in mutation_tables_mod:
        if mut.site == index_causal_locus:
            derived_states.append(mut.derived_state)
            replacement_dic[mut.derived_state] = mut.metadata


    print(pd.Series(derived_states).unique())
    random_nucleotide = random.sample(pos_nucl, 4)

    for order, key in enumerate(replacement_dic.keys()):
        replacement_dic[key]['mutation_list'][0]['nucleotide'] = random_nucleotide[order]

    
    for i, mut in enumerate(mutation_tables_mod):
        if mut.site == index_causal_locus:
            # Replace the mutation metadata and create a new mutation object
            new_mut = mut.replace(metadata=replacement_dic[mut.derived_state])
            # Update the mutation table with the new mutation object
            mutation_tables_mod[i] = new_mut

tables = ts.dump_tables()
tables.mutations.replace_with(mutation_tables_mod)
# Now you can save the modified tree sequence
ts = tables.tree_sequence()


## now conver tthe rest into nucleotide type
nts = pyslim.generate_nucleotides(ts)
nts = pyslim.convert_alleles(nts)


## create ancestral state to match derived state 
index_causal_loci = {}

sites_table = ts.tables.sites
# Print the ancestral state of all mutations
for index, site in enumerate(sites_table):
    if site.position in causal_loci:
        index_causal_loci[index] = site.position

pos_derived= {}

## now there is an ancestral state 
mutations_table = nts.tables.mutations

# Print the ancestral state of all mutations
for index, mut in enumerate(mutations_table):
    if mut.site in index_causal_loci.keys():
        #print(index)
        pos = index_causal_loci[mut.site]
        pos_derived[pos] = mut.derived_state
        ##print(f"Position: {site.position}, Ancestral State: {site.ancestral_state}")

## assign the ancestral sites to one of the 2 mutations 0 or 1 
sites_table = nts.tables.sites
# Create a new site table to hold the modified entries
new_sites_table = tskit.SiteTable()


# Iterate through the existing sites, modify the ancestral state, and add them to the new table
for site in sites_table:
    # Modify the ancestral state of the causal mutation '
    if site.position in pos_derived.keys():
        new_sites_table.add_row(
            position=site.position,
            ancestral_state= pos_derived[site.position], 
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

##### check 
for causal_locus in causal_loci:
    print(causal_locus)
   # Get the sites table
    sites_table = nts.tables.sites

    # Find the index of the site at the causal position
    index_causal_pos = None
    for index, site in enumerate(sites_table):
        if site.position == causal_locus:
            index_causal_pos = index
            print(f"Position: {site.position}, Ancestral State: {site.ancestral_state}")
            break

    if index_causal_pos is None:
        raise ValueError(f"Site at position {causal_locus} not found.")

    # Get the derived states at the causal position
    derived_states = []
    mutations_table = nts.tables.mutations
    for mut in mutations_table:
        if mut.site == index_causal_pos:
            derived_states.append(mut.derived_state)
    print(pd.Series(derived_states).unique())

    # Check if all derived states are the same
    all_same = (pd.Series(derived_states) == derived_states[0]).all()
    if all_same:
        print("All derived states are the same.")
    else:
        print("Derived states are not the same.")


    if (site.ancestral_state in derived_states) == False:
        print('derived state not in ancestral state')
    elif (site.ancestral_state in derived_states) == True:
        print('derived state in ancestral state')

## calculate diversity and save it 
genetic_diversity = nts.diversity()

with open(genetic_diversity_file, 'w') as file:
    file.write(str(genetic_diversity))

# Write the filtered tree sequence to a VCF file
with open(vcf_file, 'w') as file:
    # Pass the file object as the output parameter
    nts.write_vcf(output=file, allow_position_zero = True)


population_counts = {pop.id: 0 for pop in nts.populations()}

for ind in nts.individuals():
    # Each individual has a population
    pop_id = ind.population
    population_counts[pop_id] += 1

# Save the dictionary as a pickle file
with open(population_counts_file, 'wb') as file:
    pickle.dump(population_counts, file)
