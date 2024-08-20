import random
import json


parameter_space = snakemake.input['parameter_space'] 
causal_loci_output = snakemake.output['causal_loci'] 

# Open and read the JSON file
with open(parameter_space, 'r') as f:
    parameters = json.load(f)

number_of_causal_loci = parameters['causal_loci']

causal_loci = [random.randint(0, 500000) for _ in range(number_of_causal_loci)]
causal_loci.sort()

with open(causal_loci_output, 'w') as f:
    for loci in causal_loci:
        f.write(f"{loci}\n")
