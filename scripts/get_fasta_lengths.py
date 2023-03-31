import pandas as pd
import os

from Bio import SeqIO

FASTA_PATHS = snakemake.input[0]
OUTPUT_PATH = snakemake.output[0]
CDS_BASED_BBH = snakemake.params.CDS_BASED_BBH

names = []
lengths = []

for seq_record in SeqIO.parse(FASTA_PATHS, "fasta"):
    
    # convert record ID depending on the pipeline mode
    if CDS_BASED_BBH: 
        name = seq_record.id                        # phage protein/ORF ID
        name = '_'.join(name.split('_')[:-2])       # get phage ID
        length = len(seq_record.seq)                # get protein/ORF ID

    # whole phage is loaded
    else: 
        name = seq_record.id                # unique phageID
        length = len(seq_record.seq)        # whole phage length

    # append
    names.append(name)
    lengths.append(length)
    

# create table
length_table = pd.DataFrame({"genome": names, "length": lengths})
df_prots=length_table['genome'].value_counts()
prots_count=df_prots.to_dict()
length_table_final=length_table.groupby('genome').sum().reset_index()
length_table_final['n_prots']=length_table_final['genome'].map(prots_count)
length_table_final.to_csv(OUTPUT_PATH, index = False)

