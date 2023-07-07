import pandas as pd
import os

from Bio import SeqIO

FASTA_PATHS = snakemake.input[0]
OUTPUT_PATH = snakemake.output[0]
SEPARATOR = snakemake.params.SEPARATOR
CDS_BASED = snakemake.params.CDS_BASED

names = []
lengths = []

for seq_record in SeqIO.parse(FASTA_PATHS, "fasta"):
    
    # convert record ID depending on the pipeline mode
    if CDS_BASED: 
        name = seq_record.id                                    # phage protein/ORF ID
        name = SEPARATOR.join(name.split(SEPARATOR)[:-1])       # get phage ID
        if '|' in name:                                         # WARNING
            print(f"WARNING: '|' character in the fasta header of phage {name}!")
            print(" Please avoid this character in phage identifiers. MMseq2 removes text form identifier preceding '|' character ! ")
            print("MANIAC handles that, but it potentially can lead to errors.")
            name = name.split('|')[-1]                          # remove text preceding '|' character
        length = len(seq_record.seq)                            # get protein/ORF ID

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

# save number ORFs per genome only in cds-based mode
if CDS_BASED: length_table_final.to_csv(OUTPUT_PATH, index = False)
else: length_table_final.iloc[:, :-1].to_csv(OUTPUT_PATH, index = False)


