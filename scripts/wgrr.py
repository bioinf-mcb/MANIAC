""" Calculate whole genome repertoire relatedness (wGRR). """

# import modules
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import numpy as np

# define paths
cds_dir = snakemake.input[0]
ani_file = snakemake.input[1] # get bbh mean ani (from fragments)
best_hits_file = snakemake.input[2] # get number of bbh

wgrr_file = snakemake.output[0]
input_extension = snakemake.params.input_extension


# load data
ani_df = pd.read_csv(ani_file, sep=',')
best_hits = pd.read_csv(best_hits_file)
best_hits.rename(columns = {"query_seq_x": "query_seq", "reference_seq_x": "reference_seq"}, inplace=True)

# get number of bbh used to calculate mean ani
bbh_cols = ['query_seq', 'reference_seq', 'query_fragment_id_x', 'reference_fragment_id_x']
bbh_df = best_hits.groupby(bbh_cols).count().reset_index()[bbh_cols] # group by phages and orfs
bbh_df['bbh_counts'] = 1
bbh_df = bbh_df.groupby(['query_seq', 'reference_seq']).count()['bbh_counts'].reset_index()
bbh_df.rename({'query_seq':'Seq1', 'reference_seq': 'Seq2'}, axis=1, inplace=True)

# get number of prots for each phage
fnames = [path.stem for path in Path(cds_dir).glob(f'*.{input_extension}')]
phageIDs, n_prots = [], []
for fname in fnames:
    cds_fasta = Path(cds_dir, fname + f'.{input_extension}')
    n_prot = len(list(SeqIO.parse(cds_fasta, 'fasta')))
    phageIDs.append(fname)
    n_prots.append(n_prot)

n_prot_df = pd.DataFrame({'Seq1': phageIDs, 'n_prots': n_prots})


### calculate wgrr

# add number of prots
ani_df = pd.merge(ani_df, n_prot_df, how='left', on='Seq1')
ani_df.rename({'n_prots': 'seq1_n_prots'}, axis=1, inplace=True)

n_prot_df.rename({'Seq1': 'Seq2'}, axis=1, inplace=True)
ani_df = pd.merge(ani_df, n_prot_df, how='left', on='Seq2')
ani_df.rename({'n_prots': 'seq2_n_prots'}, axis=1, inplace=True)

ani_df = pd.merge(ani_df, bbh_df, on=['Seq1', 'Seq2'], how='left') # add number of bbh for each pair of genomes
ani_df['min_prots'] = ani_df[['seq1_n_prots', 'seq2_n_prots']].min(axis=1) # min number of prots from pair of genomes
ani_df['ani_sum'] = ani_df['ani_mean'] * ani_df['bbh_counts'] # sum of ani for pair of genomes from bbh
ani_df['wgrr'] = np.round(ani_df['ani_sum'] / ani_df['min_prots'], 3) # wgrr

# round some values
ani_df['ani_mean'] = np.round(ani_df['ani_mean'], 3)
ani_df['ani_sum'] = np.round(ani_df['ani_sum'], 3)

# save wgrr
cols = ['Seq1','Seq2','ani_mean','bbh_counts','ani_sum', \
        'seq1_n_prots','seq2_n_prots','min_prots','wgrr']

ani_df.sort_values(['wgrr','Seq1'], ascending=False, inplace=True)
ani_df[cols].to_csv(wgrr_file, sep=',', index=False)
