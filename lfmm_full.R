library("lfmm")
library("dplyr")
library("qvalue")


allele_freq_file = '/Users/tatiana/Documents/grenenet/simulations/test_slim/cline_creation/results/lfmm/delta_p_lfmm_f22a3c94_acg_gen3.csv' #.set_index('Unnamed: 0')
env_var_file = "/Users/tatiana/Documents/grenenet/simulations/test_slim/cline_creation/results/lfmm/env_vars/env_var_f22a3c94_acg_delta_p_gen3.csv"
# Access the inputs and outputs via snakemake

k_value_file <- snakemake@output[["k_value_file"]]
output_file <- snakemake@output[["lfmm_pvalues"]]

# Print the file paths (optional, for debugging)
print(paste("Allele counts file:", allele_freq_file))
print(paste("Environmental variables file:", env_var_file))
print(paste("Output file:", output_file))
print(paste("k_value_file:", k_value_file))

# Load training environmental variables dynamically based on the bio_name
env = read.csv(env_var_file, sep = ',', header = FALSE)
print(env$V1)

# myVector = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)
env_train <- matrix(env$V1, ncol = 1)
print(env_train)
print(dim(env_train))

# Check the number of rows and columns
num_pops <- dim(env_train)[1]

# Determine the value of k based on the number of rows or columns
if (num_pops < 17) {
  k <- 9
} else {
  k <- min(9, max(1, num_pops - 1))
}

print(k)
# Save the value of k to a text file
print(k_value_file)
write(k, file = k_value_file)

# Load allele frequencies
allele_freq <- read.csv(allele_freq_file, sep = ',', header = TRUE, check.names = FALSE)
colnames(allele_freq) <- NULL
allele_freq <- t(as.matrix(allele_freq))

print(dim(allele_freq))
# Model training
mod.lfmm <- lfmm_ridge(Y = allele_freq, X = env_train, K = k)

# Association testing
pv <- lfmm_test(Y = allele_freq, X = env_train, lfmm = mod.lfmm, calibrate = "gif")

## save the p values 
p_values <- pv$pvalue
write.csv(p_values, file = output_file)
