#!/usr/bin/env Rscript
# Load the necessary library
library(data.table)

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]  # Input file
output_file <- args[2]  # Output file
separator <- args[3]       # Separator for splitting (e.g., "_CDS")

# Read the input file into a data.table
dt <- fread(input_file)

# Create a new column (V9) by splitting V2 on the specified separator
dt[, V9 := sapply(strsplit(as.character(V2), separator), `[`, 1)]

# Sort the data.table by the specified keys (adjust column names as needed)
setorder(dt, V1, V9, V7)

# Write the sorted data.table to the output file
fwrite(dt, output_file, sep = "\t", col.names=F)