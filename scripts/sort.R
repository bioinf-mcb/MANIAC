#!/usr/bin/env Rscript

# Load the necessary library
library(data.table)
# Limit cpus for better parallel handling
setDTthreads(1)

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]  # Input file
input_name <- basename(input_file)
input_path <- dirname(input_file)
output_file <- file.path(input_path, paste0("sorted_", input_name))

# Read the input file into a data.table
dt <- fread(input_file)

# Sort the data.table by the specified keys (adjust column names as needed)
setorder(dt, V1, V2, V7)

# Write the sorted data.table to the output file
fwrite(dt, output_file, sep = "\t", col.names=F)