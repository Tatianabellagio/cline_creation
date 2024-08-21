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
stages = ["bcg", "acg"]

rule all:
    input:
        expand(
            "results/fst/fst_per_pos_{hash}_{stage}.csv",
            hash=open('all_hashes.txt').read().strip().split('\n'),
            stage=stages
        ),
        expand(
            "results/clines/clines_df_{hash}_{stage}.csv",
            hash=open('all_hashes.txt').read().strip().split('\n'),
            stage=stages
        ),
        expand(
            "results/lfmm/manhattan/manhattan_{hash}_{stage}.png",
            hash=open('all_hashes.txt').read().strip().split('\n'),
            stage=stages
        ),



rule def_causal_loci:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
    output:
        causal_loci="results/causal_loci/causal_loci_{hash}.txt",
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/def_causal_loci/causal_loci_{hash}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/def_causal_loci.py"

rule run_simulations:
    input:
        causal_loci="results/causal_loci/causal_loci_{hash}.txt",
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
    output:
        expand("results/tree_seq/tree_seq_{{hash}}_{stage}.trees", stage=stages),
        expand("results/phenotypes/phenotypes_{{hash}}_{stage}.txt", stage=stages),
        local_adaptation="results/local_adaptation/local_adapt_{hash}.txt",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/run_simulations/sim_{hash}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/run_slim.sh"

rule process_trees:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        tree_seq_file="results/tree_seq/tree_seq_{hash}_{stage}.trees",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        vcf_file="results/vcfs/vcf_{hash}_{stage}.vcf",
        population_counts_file="results/population_counts/population_counts_{hash}_{stage}.pkl",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/process_trees/process_trees_{hash}_{stage}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/process_tree.py"


rule gen_alleles_df:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        vcf_file="results/vcfs/vcf_{hash}_{stage}.vcf",
        population_counts_file="results/population_counts/population_counts_{hash}_{stage}.pkl",
    output:
        allele_counts_df="results/alleles_df/allele_counts_{hash}_{stage}.csv",
        allele_freq_df="results/alleles_df/allele_freq_{hash}_{stage}.csv",
        allele_counts_lfmm="results/lfmm/allele_counts_lfmm_{hash}_{stage}.csv",
        env_var_lfmm="results/lfmm/env_vars/env_var_{hash}_{stage}.csv",
        left_pos_lfmm="results/lfmm/left_pos/left_pos_{hash}_{stage}.csv"

    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/gen_alleles_df/gen_alleles_df_{hash}_{stage}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/gen_alleles_df.py"

rule fst_calculation:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        vcf_file="results/vcfs/vcf_{hash}_{stage}.vcf",
        population_counts_file="results/population_counts/population_counts_{hash}_{stage}.pkl",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        mean_fst_file="results/fst/mean_fst_{hash}_{stage}.txt",
        fst_per_pos="results/fst/fst_per_pos_{hash}_{stage}.csv",
        fst_manhattan="results/fst/fst_manhattan_{hash}_{stage}.png",
        fst_across_pops="results/fst/fst_across_pops_{hash}_{stage}.csv",
        fst_heatmap="results/fst/fst_heatmap_{hash}_{stage}.png",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/fst_calculation/fst_calculation_{hash}_{stage}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/fst_calc.py"

rule clines:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        allele_freq_df_file="results/alleles_df/allele_freq_{hash}_{stage}.csv",
        allele_counts_df_file="results/alleles_df/allele_counts_{hash}_{stage}.csv",
        population_counts_file="results/population_counts/population_counts_{hash}_{stage}.pkl",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        cline_file="results/clines/clines_df_{hash}_{stage}.csv",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/clines/clines_calc_{hash}_{stage}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/clines_calc.py"

rule run_lfmm:
    input:
        allele_counts_lfmm="results/lfmm/allele_counts_lfmm_{hash}_{stage}.csv",
        env_var_lfmm="results/lfmm/env_vars/env_var_{hash}_{stage}.csv",
    output:
        lfmm_pvalues="results/lfmm/pvalues/pvalues_{hash}_{stage}.csv",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/run_lfmm/run_lfmm_{hash}_{stage}.txt"
    conda:
        "envs/lfmm_env.yaml"
    script:
        "scripts/run_lfmm.R"

rule manhattan_lfmm:
    input:
        lfmm_pvalues="results/lfmm/pvalues/pvalues_{hash}_{stage}.csv",
        left_pos_lfmm="results/lfmm/left_pos/left_pos_{hash}_{stage}.csv",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        lfmm_manhattan="results/lfmm/manhattan/manhattan_{hash}_{stage}.png",
        sign_pos_df="results/lfmm/sign_pos/sign_pos_{hash}_{stage}.csv",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/manhattan_lfmm/manhattan_lfmm_{hash}_{stage}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/manhattan_lfmm.py"

