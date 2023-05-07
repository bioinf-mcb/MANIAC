# modules
import pandas as pd
import gc
import numpy as np

# paths
INPUT_PATH = snakemake.input[0]
LENGTHS_PATH = snakemake.input[1]
BEST_HITS_PATH = snakemake.output[0]
OUTPUT_PATH = snakemake.output[1]

# params
CDS_BASED = snakemake.params.CDS_BASED
col_dtypes = {'length':'int32','mismatches':'int32',
              'qlen':'int32','gaps':'int32',
              'ani_alnlen':'int32',
              'ani_ids':'int32'}


# load hits (RAM efficiently)
print('Load processed hits... ', end='')
best_hits = pd.read_csv(INPUT_PATH,dtype=col_dtypes) 
print('Done!')

if CDS_BASED:
    print("Finding best bidirectional hits (BBH)...", end='')
    bbh = pd.merge(best_hits, best_hits, how="inner", \
                                         left_on=["query_fragment_id", "reference_fragment_id"], \
                                         right_on=["reference_fragment_id", "query_fragment_id"]) 
    # release RAM
    del best_hits
    gc.collect()

    # deduplicate
    bbh.query("query_seq_x < reference_seq_x",inplace=True)  
    bbh["pident"] = (bbh.pident_x + bbh.pident_y) / 2
    bbh["ani_alnlen"] = (bbh.ani_alnlen_x + bbh.ani_alnlen_y) / 2
    best_hits_final = bbh
    print("Done!")
else:
    best_hits_final = best_hits

    # release RAM
    del best_hits
    gc.collect()


print("Save best hits... ", end='')
best_hits_final.to_csv(BEST_HITS_PATH, index=False)
print("Done!")

# clean best hits table 
print("HERE THERE WAS SOME PROBLEMS WITH IDS AND I DONT REMEMBER WHAT")
best_hits_final = best_hits_final.rename(columns = {"query_seq_x": "query_seq", "reference_seq_x": "reference_seq"})
best_hits_final.query_seq = best_hits_final.query_seq.str.split("_FRAGMENT_", expand = True)[0]
best_hits_final.reference_seq = best_hits_final.reference_seq.str.split("_FRAGMENT_", expand = True)[0]

print("Calculating ANI... ", end='')
# load genome lengths | n_proteins per genome
lengths_df = pd.read_csv(LENGTHS_PATH)  

# genome lengths (ANImm calculation)
genome_length_df = lengths_df[['genome','length']]
genome_length_df = genome_length_df.set_index('genome')

# calculate ANI
ani = best_hits_final.groupby(["query_seq", "reference_seq"]).pident.mean().reset_index()

# calculate different measures
aligned_nucleotides = best_hits_final.groupby(["query_seq", "reference_seq"]).ani_alnlen.sum().reset_index()
aligned_nucleotides["len_1"] = genome_length_df.loc[aligned_nucleotides.query_seq].reset_index(drop=True)
aligned_nucleotides["len_2"] = genome_length_df.loc[aligned_nucleotides.reference_seq].reset_index(drop=True)
aligned_nucleotides["af_1"] = aligned_nucleotides.ani_alnlen / aligned_nucleotides.len_1
aligned_nucleotides["af_2"] = aligned_nucleotides.ani_alnlen / aligned_nucleotides.len_2
aligned_nucleotides["af_mean"] = 2* aligned_nucleotides.ani_alnlen / (aligned_nucleotides.len_1 + aligned_nucleotides.len_2)
aligned_nucleotides["af_min"] = aligned_nucleotides.ani_alnlen / aligned_nucleotides[["len_1", "len_2"]].min(axis=1)
aligned_nucleotides["af_jaccard"] = aligned_nucleotides.ani_alnlen / (aligned_nucleotides.len_1 + aligned_nucleotides.len_2 - aligned_nucleotides.ani_alnlen)

# add measures
merged = pd.merge(ani, aligned_nucleotides, on = ["query_seq", "reference_seq"]) \
           .rename({"query_seq": "Seq1", "reference_seq": "Seq2", "pident": "ANI"}, axis=1)
print("Done!")

# save
if CDS_BASED:
    print("Calculating wGRR... ", end='')

    # proteins number per genome (wGRR calculation)
    n_prot_df = lengths_df[['genome','n_prots']].copy()
    n_prot_df = n_prot_df.rename({'genome': 'Seq1'}, axis=1)

    # get number of bbh used to calculate mean ani
    bbh_cols = ['query_seq', 'reference_seq', 'query_fragment_id_x', 'reference_fragment_id_x']
    bbh_df = best_hits_final.groupby(bbh_cols).count().reset_index()[bbh_cols] # group by phages and orfs
    bbh_df['bbh_counts'] = 1
    bbh_df = bbh_df.groupby(['query_seq', 'reference_seq']).count()['bbh_counts'].reset_index()
    bbh_df.rename({'query_seq':'Seq1', 'reference_seq': 'Seq2'}, axis=1, inplace=True)

    # merge ani df to number of proteins 
    ani_df = pd.merge(merged, n_prot_df, how='left', on='Seq1')
    ani_df.rename({'n_prots': 'seq1_n_prots'}, axis=1, inplace=True)

    n_prot_df.rename({'Seq1': 'Seq2'}, axis=1, inplace=True)
    ani_df = pd.merge(ani_df, n_prot_df, how='left', on='Seq2')
    ani_df.rename({'n_prots': 'seq2_n_prots'}, axis=1, inplace=True)

    # check point if its merged instead of ani_df
    ani_df = pd.merge(ani_df, bbh_df, on=['Seq1', 'Seq2'], how='left') # add number of bbh for each pair of genomes
    ani_df['min_prots'] = ani_df[['seq1_n_prots', 'seq2_n_prots']].min(axis=1) # min number of prots from pair of genomes

    ani_df['ani_sum'] = ani_df['ANI'] * ani_df['bbh_counts'] # sum of ani for pair of genomes from bbh
    ani_df['wgrr'] = np.round(ani_df['ani_sum'] / ani_df['min_prots'], 3) # wgrr
    ani_df.to_csv(OUTPUT_PATH, index=False)
    print("Done!")
else:
    merged.to_csv(OUTPUT_PATH, index=False)









