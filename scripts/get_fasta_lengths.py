import pandas as pd
import os

from Bio import SeqIO

FASTA_PATHS = snakemake.input
OUTPUT_PATH = snakemake.output[0]

names = []
lengths = []

for seq_record in SeqIO.parse(FASTA_PATHS, "fasta"):
    name = seq_record.id
    name = name.split('.')[0]
    names.append(name)
    length = 0
    length += len(seq_record.seq)
    	lengths.append(length)

length_table = pd.DataFrame({"genome": names, "length": lengths})
length_table_final=length_table.groupby('genome').sum().reset_index()
length_table_final.to_csv(OUTPUT_PATH, index = False)

