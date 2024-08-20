

# Check if the package is already installed
if (!requireNamespace("lfmm", quietly = TRUE)) {
  # If the package is not installed, install it
  library(devtools)
  devtools::install_github("bcm-uga/lfmm")
}

library("lfmm")
library("dplyr")
library("qvalue")



# Access the inputs and outputs via snakemake
allele_counts_file <- snakemake@input[["allele_counts_lfmm"]]
env_var_file <- snakemake@input[["env_var_lfmm"]]
output_file <- snakemake@output[["lfmm_pvalues"]]

# Print the file paths (optional, for debugging)
print(paste("Allele counts file:", allele_counts_file))
print(paste("Environmental variables file:", env_var_file))
print(paste("Output file:", output_file))

k <- 9

# Load training environmental variables dynamically based on the bio_name
env = read.csv(env_var_file, sep = ',', header = FALSE)
print(env$V1)

# myVector = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)
env_train <- matrix(env$V1, ncol = 1)
print(env_train)
print(dim(env_train))
# Load allele frequencies
allele_freq <- read.csv(allele_counts_file, sep = ',', header = TRUE, check.names = FALSE)
colnames(allele_freq) <- NULL
allele_freq <- t(as.matrix(allele_freq))

print(dim(allele_freq))
print(dim(env_train))
# Model training
mod.lfmm <- lfmm_ridge(Y = allele_freq, X = env_train, K = k)

# Association testing
pv <- lfmm_test(Y = allele_freq, X = env_train, lfmm = mod.lfmm, calibrate = "gif")

## save the p values 
p_values <- pv$pvalue
write.csv(p_values, file = output_file)
