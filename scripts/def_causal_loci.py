import random
import json


parameter_space = snakemake.input['parameter_space'] 
causal_loci_output = snakemake.output['causal_loci'] 

# Open and read the JSON file
with open(parameter_space, 'r') as f:
    parameters = json.load(f)

number_of_causal_loci = parameters['causal_loci']
genome_length = parameters['genome_length']

min_distance = 2  # Minimum distance between loci if not it will create problems in slim 

# Generate causal loci ensuring they are not consecutive
causal_loci = set()
while len(causal_loci) < number_of_causal_loci:
    new_locus = random.randint(0, genome_length-1)
    # Ensure the new locus doesn't violate the minimum distance rule
    if all(abs(new_locus - locus) >= min_distance for locus in causal_loci):
        causal_loci.add(new_locus)

# Sort the loci
causal_loci = sorted(causal_loci)

with open(causal_loci_output, 'w') as f:
    for loci in causal_loci:
        f.write(f"{loci}\n")
