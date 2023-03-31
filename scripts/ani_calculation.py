import pandas as pd
import gc
import numpy as np

try:
    INPUT_PATH = snakemake.input[0]
    LENGTHS_PATH = snakemake.input[1]
    BEST_HITS_PATH = snakemake.output[0]
    OUTPUT_PATH = snakemake.output[1]

    BBH = snakemake.params["CDS_BASED_BBH"]

except NameError:
    import argparse

    parser = argparse.ArgumentParser(description='Find best bidirectional hits from formatted MMSeqs2 results')
    parser.add_argument("input", help="Formatted MMSeqs2 results")
    parser.add_argument("fastalen", help="CSV file of fasta lengths")
    parser.add_argument("besthits", help="Output path for best hits")
    parser.add_argument("output", help="Output path for search results")
    parser.add_argument("--eval", default=10, type=float, help="E-value threshold (default: 10)")
    parser.add_argument("--coverage", default=0, type=float, help="Coverage threshold (default: 0)")
    parser.add_argument("--identity", default=0, type=float, help="Identity threshold (default: 0)")
    parser.add_argument("--bbh", default=False, action='store_true', help="Perform a CDS-based BBH ANI calculation")
    parser.set_defaults(bbh=False)

    args = parser.parse_args()

    INPUT_PATH = args.input
    LENGTHS_PATH = args.fastalen
    BEST_HITS_PATH = args.besthits
    OUTPUT_PATH = args.output

    EVAL_THR = args.eval
    COVERAGE_THR = args.coverage
    IDENTITY_THR = args.identity
    BBH = args.bbh


col_dtypes = {'length':'int32','mismatches':'int32','qlen':'int32','gaps':'int32','ani_alnlen':'int32','ani_ids':'int32'}
best_hits = pd.read_csv(INPUT_PATH,dtype=col_dtypes) 
print('dataframe loaded')


if BBH:
    print("Finding best bidirectional hits...")
    bbh = pd.merge(best_hits,best_hits, how="inner", \
      left_on=["query_fragment_id", "reference_fragment_id"], \
      right_on=["reference_fragment_id", "query_fragment_id"]) 
    del best_hits
    gc.collect()
    bbh.query("query_seq_x < reference_seq_x",inplace=True)  # Avoid duplication of each pair
    bbh["pident"] = (bbh.pident_x + bbh.pident_y) / 2
    bbh["ani_alnlen"] = (bbh.ani_alnlen_x + bbh.ani_alnlen_y) / 2
    best_hits_final = bbh
else:
    best_hits_final = best_hits
    del best_hits
    gc.collect()

best_hits_final.to_csv(BEST_HITS_PATH, index=False)
best_hits_final = best_hits_final.rename(columns = {"query_seq_x": "query_seq", "reference_seq_x": "reference_seq"})

best_hits_final.query_seq = best_hits_final.query_seq.str.split("_FRAGMENT_", expand = True)[0]
best_hits_final.reference_seq = best_hits_final.reference_seq.str.split("_FRAGMENT_", expand = True)[0]


ani = best_hits_final.groupby(["query_seq", "reference_seq"]).pident.mean().reset_index()

lengths_df=pd.read_csv(LENGTHS_PATH)  #read the fasta length csv
prot_length=lengths_df[['genome','length']]  #get a subset of the fasta_length for ani calculation
n_prot_df =lengths_df[['genome','n_prots']].copy()  #get a subset with the number of proteins for wgrr calculation 
n_prot_df.rename({'genome': 'Seq1'}, axis=1, inplace=True)  #renaming the genome column
fasta_lengths = prot_length.set_index('genome')



aligned_nucleotides = best_hits_final.groupby(["query_seq", "reference_seq"]).ani_alnlen.sum().reset_index()
aligned_nucleotides["len_1"] = fasta_lengths.loc[aligned_nucleotides.query_seq].reset_index(drop=True)
aligned_nucleotides["len_2"] = fasta_lengths.loc[aligned_nucleotides.reference_seq].reset_index(drop=True)
aligned_nucleotides["af_1"] = aligned_nucleotides.ani_alnlen / aligned_nucleotides.len_1
aligned_nucleotides["af_2"] = aligned_nucleotides.ani_alnlen / aligned_nucleotides.len_2
aligned_nucleotides["af_mean"] = 2* aligned_nucleotides.ani_alnlen / (aligned_nucleotides.len_1 + aligned_nucleotides.len_2)
aligned_nucleotides["af_min"] = aligned_nucleotides.ani_alnlen / aligned_nucleotides[["len_1", "len_2"]].min(axis=1)
aligned_nucleotides["af_jaccard"] = aligned_nucleotides.ani_alnlen / (aligned_nucleotides.len_1 + aligned_nucleotides.len_2 - aligned_nucleotides.ani_alnlen)

merged = pd.merge(ani, aligned_nucleotides, on = ["query_seq", "reference_seq"]) \
           .rename({"query_seq": "Seq1", "reference_seq": "Seq2", "pident": "ANI"}, axis=1)

if (BBH):
    print("wgrr stage")
    #####calculating wgrr
    #get number of bbh used to calculate mean ani
    bbh_cols = ['query_seq', 'reference_seq', 'query_fragment_id_x', 'reference_fragment_id_x']
    bbh_df = best_hits_final.groupby(bbh_cols).count().reset_index()[bbh_cols] # group by phages and orfs
    bbh_df['bbh_counts'] = 1
    bbh_df = bbh_df.groupby(['query_seq', 'reference_seq']).count()['bbh_counts'].reset_index()
    bbh_df.rename({'query_seq':'Seq1', 'reference_seq': 'Seq2'}, axis=1, inplace=True)

    #merge ani df to number of proteins 
    ani_df = pd.merge(merged, n_prot_df, how='left', on='Seq1')
    ani_df.rename({'n_prots': 'seq1_n_prots'}, axis=1, inplace=True)

    n_prot_df.rename({'Seq1': 'Seq2'}, axis=1, inplace=True)
    ani_df = pd.merge(ani_df, n_prot_df, how='left', on='Seq2')
    ani_df.rename({'n_prots': 'seq2_n_prots'}, axis=1, inplace=True)

    ##check point if its merged instead of ani_df
    ani_df = pd.merge(ani_df, bbh_df, on=['Seq1', 'Seq2'], how='left') # add number of bbh for each pair of genomes
    ani_df['min_prots'] = ani_df[['seq1_n_prots', 'seq2_n_prots']].min(axis=1) # min number of prots from pair of genomes

    ani_df['ani_sum'] = ani_df['ANI'] * ani_df['bbh_counts'] # sum of ani for pair of genomes from bbh
    ani_df['wgrr'] = np.round(ani_df['ani_sum'] / ani_df['min_prots'], 3) # wgrr

    ani_df.to_csv(OUTPUT_PATH, index=False)
else:
    merged.to_csv(OUTPUT_PATH, index=False)

print("Done!")







