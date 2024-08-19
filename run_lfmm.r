library("lfmm")
library("dplyr")
library("qvalue")
#library('data.table')

args <- commandArgs(trailingOnly = TRUE)
runid <- args[1]
k <- as.integer(args[2])

delta_p_file <- paste0("allele_counts/all_alt_allele_counts_runid_", runid, ".csv")
env_file = paste0("env_files/env_var_", runid, ".csv")
print(env_file)
# Load training environmental variables dynamically based on the bio_name
env = read.csv(env_file, sep = ',', header = FALSE)
print(env$V1)

# myVector = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)
env_train <- matrix(env$V1, ncol = 1)
print(env_train)
print(dim(env_train))
# Load allele frequencies
allele_freq <- read.csv(delta_p_file, sep = ',', header = TRUE, check.names = FALSE)
colnames(allele_freq) <- NULL
allele_freq <- t(as.matrix(allele_freq))


output_file_wo_calibration_pvalue <- paste0('results/wo_calibration_pvalue_full_genome_' , runid, '_k', k, ".csv")
output_file_effect_sizes_simple <- paste0('results/effect_sizes_simple_full_genome_' , runid, '_k', k,".csv")
print(k)

print(dim(allele_freq))
print(dim(env_train))
# Model training
mod.lfmm <- lfmm_ridge(Y = allele_freq, X = env_train, K = k)

# Association testing
pv <- lfmm_test(Y = allele_freq, X = env_train, lfmm = mod.lfmm, calibrate = "gif")

## save the p values 
p_values <- pv$pvalue
write.csv(p_values, file = output_file_wo_calibration_pvalue)

## save the p values 
#calibrated_pvalues <- pv$calibrated.pvalue
#write.csv(calibrated_pvalues, file = output_file_w_calibration_pvalue)

## save the betas 
beta_i <- pv$B
write.csv(beta_i, file = output_file_effect_sizes_simple)

#print(head(pv$B))
# Beta coefficients for the biovariable
#biovar_beta = matrix(pv$B)

