import pandas as pd
import os

from Bio import SeqIO

FASTA_PATHS = snakemake.input
OUTPUT_PATH = snakemake.output[0]

names = []
lengths = []

for fasta_path in FASTA_PATHS:
    name = os.path.basename(fasta_path)
    name = ".".join(name.split(".")[:-1])
    names.append(name)
    length = 0

    for record in SeqIO.parse(fasta_path, "fasta"):
        length += len(record.seq)
    
    lengths.append(length)

length_table = pd.DataFrame({"genome": names, "length": lengths})
length_table.to_csv(OUTPUT_PATH, index = False)