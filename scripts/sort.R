#!/usr/bin/env Rscript
# Load the necessary library
library(data.table)

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]  # Input file
output_file <- args[2]  # Output file

# Read the input file into a data.table
dt <- fread(input_file)

# Sort the data.table by the specified keys (adjust column names as needed)
setorder(dt, V1, V2, V7)

# Write the sorted data.table to the output file
fwrite(dt, output_file, sep = "\t", col.names=F)