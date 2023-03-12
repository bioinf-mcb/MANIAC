import pandas as pd
import os

from Bio import SeqIO

FASTA_PATHS = snakemake.input[0]
OUTPUT_PATH = snakemake.output[0]

names = []
lengths = []

for seq_record in SeqIO.parse(FASTA_PATHS, "fasta"):
    name = seq_record.id
    name = name.split('_')
    #ref_seq sequences usually have NC_.. we dont want to exclude them
    if len(name) <= 2:
        name ='_'.join(name)
        names.append(name)
        length = 0
        length += len(seq_record.seq)
        lengths.append(length)
    else:
        #for our labs convension where we have _ORF_number or _Protein_number
        name ='_'.join(name[:-2])
        names.append(name)
        length = 0
        length += len(seq_record.seq)
        lengths.append(length)

length_table = pd.DataFrame({"genome": names, "length": lengths})
df_prots=length_table['genome'].value_counts()
prots_count=df_prots.to_dict()
length_table_final=length_table.groupby('genome').sum().reset_index()
length_table_final['n_prots']=length_table_final['genome'].map(prots_count)
length_table_final.to_csv(OUTPUT_PATH, index = False)

