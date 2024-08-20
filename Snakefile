# =================================================================================================
#     Dependencies Cline Creation 
# =================================================================================================

#
# Simulation parameters 
#
# repr_scheme
# mutrate
# capacity_charge
# initial_pop_size
# subpop_number
# migration_rate
# BURNIN_1
# BURNIN_2
# last_gen
# final_sel_strength
# active_loci
# h2


configfile: "config.yaml"

rule all:
    input:
        expand(
            "results/alleles_df/allele_counts_{hash}.csv",
            hash=open('all_hashes.txt').read().strip().split('\n')),
         expand(
            "results/fst/fst_heatmap_{hash}.png",
            hash=open('all_hashes.txt').read().strip().split('\n')),
        expand(
            "results/clines/clines_df_{hash}.csv",
            hash=open('all_hashes.txt').read().strip().split('\n')),
        expand(
            "results/lfmm/sign_pos/sign_pos_{hash}.csv",
            hash=open('all_hashes.txt').read().strip().split('\n')),



rule def_causal_loci:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
    output:
        causal_loci="results/causal_loci/causal_loci_{hash}.txt",
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/def_causal_loci/causal_loci_{hash}.txt"
    script:
        "scripts/def_causal_loci.py"

rule run_simulations:
    input:
        causal_loci="results/causal_loci/causal_loci_{hash}.txt",
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
    output:
        tree_seq="results/tree_seq/tree_seq_{hash}.trees",
        phenotypes="results/phenotypes/phenotypes_{hash}.txt",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/run_simulations/sim_{hash}.txt"
    script:
        "scripts/run_slim.sh"

rule process_trees:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        tree_seq_file="results/tree_seq/tree_seq_{hash}.trees",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        vcf_file="results/vcfs/vcf_{hash}.vcf",
        population_counts_file="results/population_counts/population_counts_{hash}.pkl",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/process_trees/process_trees_{hash}.txt"
    script:
        "scripts/process_tree.py"


rule gen_alleles_df:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        vcf_file="results/vcfs/vcf_{hash}.vcf",
        population_counts_file="results/population_counts/population_counts_{hash}.pkl",
    output:
        allele_counts_df="results/alleles_df/allele_counts_{hash}.csv",
        allele_freq_df="results/alleles_df/allele_freq_{hash}.csv",
        allele_counts_lfmm="results/lfmm/allele_counts_lfmm_{hash}.csv",
        env_var_lfmm="results/lfmm/env_vars/env_var_{hash}.csv",
        left_pos_lfmm="results/lfmm/left_pos/left_pos_{hash}.csv"

    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/gen_alleles_df/gen_alleles_df_{hash}.txt"
    script:
        "scripts/gen_alleles_df.py"

rule fst_calculation:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        vcf_file="results/vcfs/vcf_{hash}.vcf",
        population_counts_file="results/population_counts/population_counts_{hash}.pkl",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        mean_fst_file="results/fst/mean_fst_{hash}.txt",
        fst_per_pos="results/fst/fst_per_pos_{hash}.csv",
        fst_manhattan="results/fst/fst_manhattan_{hash}.png",
        fst_across_pops="results/fst/fst_across_pops_{hash}.csv",
        fst_heatmap="results/fst/fst_heatmap_{hash}.png",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/fst_calculation/fst_calculation_{hash}.txt"
    script:
        "scripts/fst_calc.py"

rule clines:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        allele_freq_df_file="results/alleles_df/allele_freq_{hash}.csv",
        allele_counts_df_file="results/alleles_df/allele_counts_{hash}.csv",
        population_counts_file="results/population_counts/population_counts_{hash}.pkl",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        cline_file="results/clines/clines_df_{hash}.csv",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/clines/clines_calc_{hash}.txt"
    script:
        "scripts/clines_calc.py"

rule run_lfmm:
    input:
        allele_counts_lfmm="results/lfmm/allele_counts_lfmm_{hash}.csv",
        env_var_lfmm="results/lfmm/env_vars/env_var_{hash}.csv",
    output:
        lfmm_pvalues="results/lfmm/pvalues/pvalues_{hash}.csv",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/run_lfmm/run_lfmm_{hash}.txt"
    script:
        "scripts/run_lfmm.R"

rule manhattan_lfmm:
    input:
        lfmm_pvalues="results/lfmm/pvalues/pvalues_{hash}.csv",
        left_pos_lfmm="results/lfmm/left_pos/left_pos_{hash}.csv",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        lfmm_manhattan="results/lfmm/manhattan/manhattan_{hash}.png",
        sign_pos_df="results/lfmm/sign_pos/sign_pos_{hash}.csv",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/manhattan_lfmm/manhattan_lfmm_{hash}.txt"
    script:
        "scripts/manhattan_lfmm.py"

