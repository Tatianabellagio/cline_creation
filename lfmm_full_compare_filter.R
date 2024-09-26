library("lfmm")
library("dplyr")
library("qvalue")

# Get the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# args[2] will be the space-separated train samples
allele_freq_file <- args[1]
env_var_file <- args[2]
output_file <- args[3]
# Access the inputs and outputs via snakemake

# Print the file paths (optional, for debugging)
print(paste("Allele counts file:", allele_freq_file))
print(paste("Environmental variables file:", env_var_file))
print(paste("Output file:", output_file))

# Load training environmental variables dynamically based on the bio_name
env = read.csv(env_var_file, sep = ',', header = FALSE)
print(env$V1)

# myVector = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)
env_train <- matrix(env$V1, ncol = 1)
print(env_train)
print(dim(env_train))

# Check the number of rows and columns
num_pops <- dim(env_train)[1]

k = 9 

# Save the value of k to a text file

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
print('done')
