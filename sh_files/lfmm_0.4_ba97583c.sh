#!/bin/bash
#SBATCH --job-name=lfmm_0.4_ba97583c
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=21  
#SBATCH --output=lfmm_0.4_ba97583c-%j.out       
#SBATCH --error=lfmm_0.4_ba97583c-%j.log         
#SBATCH --mail-user=tbellagio@carnegiescience.edu
#SBATCH --mail-type=FAIL

echo start
source /home/tbellagio/miniforge3/etc/profile.d/conda.sh

export PATH="${PATH}:/home/username/bin"

cd '/carnegie/nobackup/scratch/tbellagio/simulations/cline_creation'

conda activate /home/tbellagio/miniforge3/envs/r-environment
    Rscript lfmm_full_compare_filter.R "results/test_filter_lfmm/allele_freq_filter_0.4_ba97583c.csv" "results/lfmm/env_vars/env_var_ba97583c_acg_gen3.csv" "results/test_filter_lfmm/pvalues_filter_ba97583c_0.4.csv"
    