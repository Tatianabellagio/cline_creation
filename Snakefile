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
generations = ["gen0", "gen1", "gen2", "gen3"]
delta_pfiles = ["delta_gen1", "delta_gen3"]

rule all:
    input:
        expand(
            "results/fst/fst_per_pos_{hash}_bcg_gen0.csv",
            hash=open('all_hashes.txt').read().strip().split('\n'),
            stage=["bcg", "acg"],
        ),
        expand(
            "results/clines/clines_df_{hash}_bcg_gen0.csv",
            hash=open('all_hashes.txt').read().strip().split('\n'),
            stage=["bcg", "acg"],
        ),
        expand(
            "results/lfmm/manhattan/manhattan_{hash}_bcg_gen0.png",
            hash=open('all_hashes.txt').read().strip().split('\n'),
        ),
        expand(
            "results/fst/fst_per_pos_{hash}_acg_{gen}.csv",
            hash=open('all_hashes.txt').read().strip().split('\n'),
            gen=generations
        ),
        expand(
            "results/clines/clines_df_{hash}_acg_{gen}.csv",
            hash=open('all_hashes.txt').read().strip().split('\n'),
            gen=generations
        ),
        expand(
            "results/lfmm/manhattan/manhattan_{hash}_acg_{gen}.png",
            hash=open('all_hashes.txt').read().strip().split('\n'),
            gen=generations
        ),
        expand(
            "results/lfmm/manhattan/manhattan_deltap_{hash}_acg_{gen_filter}.png",
            hash=open('all_hashes.txt').read().strip().split('\n'),
            gen_filter=["gen1", "gen3"]
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

rule process_trees_bcg:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        tree_seq_file="results/tree_seq/tree_seq_{hash}_bcg.trees",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        vcf_files=temp("results/vcfs/vcf_{hash}_bcg_gen0_multiallelic.vcf"),
        population_counts_files="results/population_counts/population_counts_{hash}_bcg_gen0.pkl",
        genetic_diversity_files="results/genetic_diversity/genetic_diversity_{hash}_bcg_gen0.txt",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/process_trees/process_trees_{hash}_bcg_gen0.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/process_tree_bcg.py"

rule process_trees_acg:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        tree_seq_file="results/tree_seq/tree_seq_{hash}_acg.trees",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        vcf_files=temp(expand("results/vcfs/vcf_{{hash}}_acg_{gen}_multiallelic.vcf", gen=["gen0", "gen1", "gen2", "gen3"])),
        population_counts_files=expand("results/population_counts/population_counts_{{hash}}_acg_{gen}.pkl", gen=["gen0", "gen1", "gen2", "gen3"]), 
        #genetic_diversity_files=expand("results/genetic_diversity/genetic_diversity_{{hash}}_acg_{gen}.txt", gen=generations), 
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/process_trees/process_trees_{hash}_bcg_gen0.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/process_tree_acg.py"

rule process_vcf_multiallelic:
    input:
        vcf_file_multiallelic="results/vcfs/vcf_{hash}_{stage}_{gen}_multiallelic.vcf",
    output:
        vcf_file="results/vcfs/vcf_{hash}_{stage}_{gen}.vcf",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/process_vcf_multiallelic/process_vcf_multiallelic_{hash}_{stage}_{gen}.txt"
    conda:
        "envs/base_env.yaml"
    shell:
        """
        bcftools norm --multiallelics -any -o {output.vcf_file} {input.vcf_file_multiallelic}
        """

rule gen_alleles_df:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        vcf_file="results/vcfs/vcf_{hash}_{stage}_{gen}.vcf",
        population_counts_file="results/population_counts/population_counts_{hash}_{stage}_{gen}.pkl",
    output:
        input_snp_number_file="results/input_snp_number/input_snp_number_{hash}_{stage}_{gen}.txt",
        allele_counts_df="results/alleles_df/allele_counts_{hash}_{stage}_{gen}.csv",
        allele_freq_df="results/alleles_df/allele_freq_{hash}_{stage}_{gen}.csv",
        allele_freq_lfmm="results/lfmm/allele_freq_lfmm_{hash}_{stage}_{gen}.csv",
        env_var_lfmm="results/lfmm/env_vars/env_var_{hash}_{stage}_{gen}.csv",
        left_pos_lfmm="results/lfmm/left_pos/left_pos_{hash}_{stage}_{gen}.csv"
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/gen_alleles_df/gen_alleles_df_{hash}_{stage}_{gen}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/gen_alleles_df.py"

rule fst_calculation:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        vcf_file="results/vcfs/vcf_{hash}_{stage}_{gen}.vcf",
        population_counts_file="results/population_counts/population_counts_{hash}_{stage}_{gen}.pkl",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        mean_fst_file="results/fst/mean_fst_{hash}_{stage}_{gen}.txt",
        fst_per_pos="results/fst/fst_per_pos_{hash}_{stage}_{gen}.csv",
        fst_manhattan="results/fst/fst_manhattan_{hash}_{stage}_{gen}.png",
        fst_across_pops="results/fst/fst_across_pops_{hash}_{stage}_{gen}.csv",
        fst_heatmap="results/fst/fst_heatmap_{hash}_{stage}_{gen}.png",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/fst_calculation/fst_calculation_{hash}_{stage}_{gen}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/fst_calc.py"

rule clines:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        allele_freq_df_file="results/alleles_df/allele_freq_{hash}_{stage}_{gen}.csv",
        allele_counts_df_file="results/alleles_df/allele_counts_{hash}_{stage}_{gen}.csv",
        population_counts_file="results/population_counts/population_counts_{hash}_{stage}_{gen}.pkl",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        cline_file="results/clines/clines_df_{hash}_{stage}_{gen}.csv",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/clines/clines_calc_{hash}_{stage}_{gen}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/clines_calc.py"

rule gen_delta_p:
    input:
        parameter_space=lambda wildcards: f"parameters/{wildcards.hash}.json",
        allele_freq_gen0="results/alleles_df/allele_freq_{hash}_acg_gen0.csv",
        allele_freq_df_file="results/alleles_df/allele_freq_{hash}_acg_{gen_filter}.csv",
        population_counts_file="results/population_counts/population_counts_{hash}_acg_{gen_filter}.pkl",
        allele_counts_df_file="results/alleles_df/allele_counts_{hash}_acg_{gen_filter}.csv",
    output:
        delta_allele_freq_lfmm="results/lfmm/delta_p_lfmm_{hash}_acg_{gen_filter}.csv",
        delta_p_env_var_lfmm="results/lfmm/env_vars/env_var_{hash}_acg_delta_p_{gen_filter}.csv",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/gen_delta_p/gen_delta_p_{hash}_{gen_filter}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/gen_delta_p.py"

rule run_lfmm:
    input:
        allele_freq_lfmm="results/lfmm/allele_freq_lfmm_{hash}_{stage}_{gens}.csv",
        env_var_lfmm="results/lfmm/env_vars/env_var_{hash}_{stage}_{gens}.csv",
    output:
        lfmm_pvalues="results/lfmm/pvalues/pvalues_{hash}_{stage}_{gens}.csv",
        k_value_file="results/lfmm/k_value/k_value_{hash}_{stage}_{gens}.txt",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/run_lfmm/run_lfmm_{hash}_{stage}_{gens}.txt"
    conda:
        "envs/lfmm_env.yaml"
    script:
        "scripts/run_lfmm.R"

rule manhattan_lfmm:
    input:
        lfmm_pvalues="results/lfmm/pvalues/pvalues_{hash}_{stage}_{gens}.csv",
        left_pos_lfmm="results/lfmm/left_pos/left_pos_{hash}_{stage}_{gens}.csv",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        lfmm_manhattan="results/lfmm/manhattan/manhattan_{hash}_{stage}_{gens}.png",
        sign_pos_df="results/lfmm/sign_pos/sign_pos_{hash}_{stage}_{gens}.csv",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/manhattan_lfmm/manhattan_lfmm_{hash}_{stage}_{gens}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/manhattan_lfmm.py"

rule run_lfmm_deltap:
    input:
        allele_freq_lfmm="results/lfmm/delta_p_lfmm_{hash}_acg_{gen_filter}.csv",
        env_var_lfmm="results/lfmm/env_vars/env_var_{hash}_acg_delta_p_{gen_filter}.csv",
    output:
        lfmm_pvalues="results/lfmm/pvalues/pvalues_deltap_{hash}_acg_{gen_filter}.csv",
        k_value_file="results/lfmm/k_value/k_value_deltap_{hash}_acg_{gen_filter}.txt",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/run_lfmm/run_lfmm_{hash}_acg_{gen_filter}.txt"
    conda:
        "envs/lfmm_env.yaml"
    script:
        "scripts/run_lfmm.R"

rule manhattan_lfmm_deltap:
    input:
        lfmm_pvalues="results/lfmm/pvalues/pvalues_deltap_{hash}_acg_{gen_filter}.csv",
        left_pos_lfmm="results/lfmm/left_pos/left_pos_{hash}_acg_{gen_filter}.csv",
        causal_loci_file="results/causal_loci/causal_loci_{hash}.txt",
    output:
        lfmm_manhattan="results/lfmm/manhattan/manhattan_deltap_{hash}_acg_{gen_filter}.png",
        sign_pos_df="results/lfmm/sign_pos/sign_pos_deltap_{hash}_acg_{gen_filter}.csv",
    params:
        hash=lambda wildcards: wildcards.hash,
    resources:
        mem_mb=30720,
    benchmark:
        "benchmarks/manhattan_lfmm/manhattan_lfmm_{hash}_acg_{gen_filter}.txt"
    conda:
        "envs/base_env.yaml"
    script:
        "scripts/manhattan_lfmm.py"