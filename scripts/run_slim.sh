
## This script requires jq to be installed to read JSON files
## Get the hash value

hash_value="${snakemake_params[hash]}"
echo $hash_value
parameter_space="${snakemake_input[parameter_space]}"
echo $parameter_space

# Extract values from the JSON using jq
repr_scheme=$(jq -r '.repr_scheme' "$parameter_space")
mutrate=$(jq -r '.mutrate' "$parameter_space")
recomb_rate=$(jq -r '.recomb_rate' "$parameter_space")
h2=$(jq -r '.h2' "$parameter_space")
capacity_charge=$(jq -r '.capacity_charge' "$parameter_space")
initial_pop_size=$(jq -r '.initial_pop_size' "$parameter_space")
subpop_number=$(jq -r '.subpop_number' "$parameter_space")
migration_rate=$(jq -r '.migration_rate' "$parameter_space")
BURNIN_1=$(jq -r '.BURNIN_1' "$parameter_space")
BURNIN_2=$(jq -r '.BURNIN_2' "$parameter_space")
cg_duration=$(jq -r '.cg_duration' "$parameter_space")
initial_sel_strength=$(jq -r '.initial_sel_strength' "$parameter_space")
final_sel_strength=$(jq -r '.final_sel_strength' "$parameter_space")

##common garden set up
common_garden_number=$(jq -r '.common_garden_number' "$parameter_space")
common_garden_migrants=$(jq -r '.common_garden_migrants' "$parameter_space")
#sel_strength_CG=$(jq -r '.sel_strength_CG' "$parameter_space")
COMMON_GARDEN_CYCLE=$(jq -r '.COMMON_GARDEN_CYCLE' "$parameter_space")


echo $hash_value
# Run SLiM with the extracted parameters as floats
slim \
    -d "hash_value='$hash_value'" \
    -d "repr_scheme='$repr_scheme'" \
    -d "mutrate=$mutrate" \
    -d "recomb_rate=$recomb_rate" \
    -d "h2=$h2" \
    -d "capacity_charge=$capacity_charge" \
    -d "initial_pop_size=$initial_pop_size" \
    -d "subpop_number=$subpop_number" \
    -d "migration_rate=$migration_rate" \
    -d "BURNIN_1=$BURNIN_1" \
    -d "BURNIN_2=$BURNIN_2" \
    -d "cg_duration=$cg_duration" \
    -d "initial_sel_strength=$initial_sel_strength" \
    -d "final_sel_strength=$final_sel_strength" \
    -d "common_garden_number=$common_garden_number" \
    -d "common_garden_migrants=$common_garden_migrants" \
    -d "COMMON_GARDEN_CYCLE=$COMMON_GARDEN_CYCLE" \
    scripts/cline_creation.slim


    #-d "sel_strength_CG=$sel_strength_CG" \