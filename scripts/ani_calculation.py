# modules
from pathlib import Path
import pandas as pd
import gc
import numpy as np
import shutil
import os


# paths
SIGNIGICANT_HITS_PATH = snakemake.input[0]
PHAGE_LENGTHS_PATH = snakemake.input[1]

GENOME_ALIGNMENT = snakemake.output[0]
CDS_ALIGNMENT_RAW = snakemake.params.CDS_ALIGNMENT_RAW
CDS_ALIGNMENT_FILE = snakemake.params.CDS_ALIGNMENT_FILE
MMSEQS_TEMP_DIR = snakemake.params.MMSEQS_TEMP_DIR

# params
CDS_BASED = snakemake.params.CDS_BASED
CONSERVED_IDENTITY = snakemake.params.CONSERVED_IDENTITY # only used in CDS BASED
CONSERVED_COVERAGE = snakemake.params.CONSERVED_COVERAGE # only used in CDS BASED

DELETE_CDS_ALIGNMENT = snakemake.params.DELETE_CDS_ALIGNMENT
DELETE_INTERMEDIATE_FILES = snakemake.params.DELETE_INTERMEDIATE_FILES
INTERMEDIATE_FILES_DIR = Path(snakemake.input[0]).parent

col_dtypes = {'length':'int32','mismatches':'int32',
              'qlen':'int32','gaps':'int32',
              'ani_alnlen':'int32',
              'ani_ids':'int32'}


# load hits (RAM efficiently)
print('Load processed hits... ', end='')
significant_hits_df = pd.read_csv(SIGNIGICANT_HITS_PATH, dtype=col_dtypes) 
print('Done!')

if CDS_BASED:
    print("Finding best bidirectional hits (BBH)...", end='')
    bbh = pd.merge(significant_hits_df, significant_hits_df, how="inner", \
                                         left_on=["query_fragment_id", "reference_fragment_id"], \
                                         right_on=["reference_fragment_id", "query_fragment_id"]) 
    # release RAM
    del significant_hits_df
    gc.collect()

    # deduplicate
    bbh.query("query_seq_x < reference_seq_x", inplace=True)  
    bbh["pident"] = np.round((bbh.pident_x + bbh.pident_y) / 2, 6)
    bbh["ani_alnlen"] = np.round((bbh.ani_alnlen_x + bbh.ani_alnlen_y) / 2, 6)
    best_hits_final = bbh
    print("Done!")
else:
    best_hits_final = significant_hits_df

    # release RAM
    del significant_hits_df
    gc.collect()


print("Save best hits... ", end='')
best_hits_final.to_csv(CDS_ALIGNMENT_RAW, index=False)
print("Done!")

# clean best hits table 
best_hits_final = best_hits_final.rename(columns = {"query_seq_x": "query_seq", "reference_seq_x": "reference_seq"})
best_hits_final.query_seq = best_hits_final.query_seq.str.split("_FRAGMENT_", expand = True)[0]
best_hits_final.reference_seq = best_hits_final.reference_seq.str.split("_FRAGMENT_", expand = True)[0]

print("Calculating ANI... ", end='')
# load genome lengths | n_proteins per genome
lengths_df = pd.read_csv(PHAGE_LENGTHS_PATH)  

# genome lengths (ANImm calculation)
genome_length_df = lengths_df[['genome','length']]
genome_length_df = genome_length_df.set_index('genome')

# calculate ANI
ani = best_hits_final.groupby(["query_seq", "reference_seq"]).pident.mean().reset_index()

# calculate different measures
aligned_nucleotides = best_hits_final.groupby(["query_seq", "reference_seq"]).ani_alnlen.sum().reset_index()
aligned_nucleotides["len_1"] = genome_length_df.loc[aligned_nucleotides.query_seq].reset_index(drop=True).astype(int)
aligned_nucleotides["len_2"] = genome_length_df.loc[aligned_nucleotides.reference_seq].reset_index(drop=True).astype(int)
aligned_nucleotides["af_1"] = np.round(aligned_nucleotides.ani_alnlen / aligned_nucleotides.len_1, 6)
aligned_nucleotides["af_2"] = np.round(aligned_nucleotides.ani_alnlen / aligned_nucleotides.len_2, 6)
aligned_nucleotides["af_mean"] = np.round(2* aligned_nucleotides.ani_alnlen / (aligned_nucleotides.len_1 + aligned_nucleotides.len_2), 6)
aligned_nucleotides["af_min"] = np.round(aligned_nucleotides.ani_alnlen / aligned_nucleotides[["len_1", "len_2"]].min(axis=1), 6)
aligned_nucleotides["af_max"] = np.round(aligned_nucleotides.ani_alnlen / aligned_nucleotides[["len_1", "len_2"]].max(axis=1), 6)
aligned_nucleotides["af_jaccard"] = np.round(aligned_nucleotides.ani_alnlen / (aligned_nucleotides.len_1 + aligned_nucleotides.len_2 - aligned_nucleotides.ani_alnlen), 6)

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

    # get number of BBH used to calculate mean ani
    cds_alignment_cols = ['query_seq', 'reference_seq', 'query_fragment_id_x', 'reference_fragment_id_x', 'pident_x', 'ani_cov_x', 'pident_y', 'ani_cov_y']
    cds_alignment_df = best_hits_final.groupby(cds_alignment_cols).count().reset_index()[cds_alignment_cols] # group by phages and orfs
    
    # rename columns
    cds_alignment_df = cds_alignment_df.rename(columns = {'query_seq':'Seq1', 'reference_seq': 'Seq2',
                                      'query_fragment_id_x': 'seq1_fragment_id', 'reference_fragment_id_x': 'seq2_fragment_id',
                                      'pident_x': 'seq1_fragment_pident', 'ani_cov_x': 'seq1_fragment_cov',
                                      'pident_y': 'seq2_fragment_pident', 'ani_cov_y': 'seq2_fragment_cov'})
    

    # calculate homologous proteins per phage
    cds_alignment_df['cds_alignments_counts'] = 1

    ### seq1
    seq1_homologous_prots_df = cds_alignment_df.groupby(['Seq1', 'Seq2'])['seq1_fragment_id'].size() \
                                                                                             .reset_index() \
                                                                                             .rename({'seq1_fragment_id': 'seq1_n_prots_hom'}, axis=1)[['Seq1', 'Seq2', 'seq1_n_prots_hom']]
    
    ### seq2                                                                                   
    seq2_homologous_prots_df = cds_alignment_df.groupby(['Seq2', 'Seq1'])['seq2_fragment_id'].size() \
                                                                                             .reset_index() \
                                                                                             .rename({'seq2_fragment_id': 'seq2_n_prots_hom'}, axis=1)[['Seq2', 'Seq1', 'seq2_n_prots_hom']]

    # merge
    homologous_prots_df = pd.merge(seq1_homologous_prots_df, seq2_homologous_prots_df, how='outer', on=['Seq1', 'Seq2'])


    # calculate conserved proteins per phage
    filt_conserved_identity = (cds_alignment_df['seq1_fragment_pident'] >= CONSERVED_IDENTITY) & (cds_alignment_df['seq2_fragment_pident'] >= CONSERVED_IDENTITY)
    filt_conserved_coverage = (cds_alignment_df['seq1_fragment_cov'] >= CONSERVED_COVERAGE) & (cds_alignment_df['seq2_fragment_cov'] >= CONSERVED_COVERAGE)
    conserved_df = cds_alignment_df.loc[filt_conserved_identity & filt_conserved_coverage]

    seq1_conserved_prots_df = cds_alignment_df.groupby(['Seq1', 'Seq2'])['seq1_fragment_id'].size() \
                                                                                            .reset_index() \
                                                                                            .rename({'seq1_fragment_id': 'seq1_n_prots_cons'}, axis=1)[['Seq1', 'Seq2', 'seq1_n_prots_cons']]


    seq2_conserved_prots_df = cds_alignment_df.groupby(['Seq2', 'Seq1'])['seq2_fragment_id'].size() \
                                                                                            .reset_index() \
                                                                                            .rename({'seq2_fragment_id': 'seq2_n_prots_cons'}, axis=1)[['Seq2', 'Seq1', 'seq2_n_prots_cons']]
    
    # merge
    conserved_prots_df = pd.merge(seq1_conserved_prots_df, seq2_conserved_prots_df, how='outer', on=['Seq1', 'Seq2'])


    # save BBH
    cols2save = ['Seq1', 'Seq2', 'seq1_fragment_id', 'seq2_fragment_id', 'seq1_fragment_pident', 'seq2_fragment_pident', 'seq1_fragment_cov', 'seq2_fragment_cov']
    #'seq1_n_prots_hom', 'seq2_n_prots_hom', 'seq1_n_prots_cons', 'seq2_n_prots_cons']

    if not DELETE_CDS_ALIGNMENT: 
        print('Delete CDS alignment file... ', end='')
        cds_alignment_df[cols2save].to_csv(CDS_ALIGNMENT_FILE, index=False)
        print('Done!')

    # count BBH
    cds_alignment_df = cds_alignment_df.groupby(['Seq1', 'Seq2']).count()['cds_alignments_counts'].reset_index()

    ### add number protein types to genome alignments

    # total number of proteins
    genome_alignment_df = pd.merge(merged, n_prot_df, how='left', on='Seq1')
    genome_alignment_df.rename({'n_prots': 'seq1_n_prots'}, axis=1, inplace=True)

    n_prot_df.rename({'Seq1': 'Seq2'}, axis=1, inplace=True)
    genome_alignment_df = pd.merge(genome_alignment_df, n_prot_df, how='left', on='Seq2')
    genome_alignment_df.rename({'n_prots': 'seq2_n_prots'}, axis=1, inplace=True)

    # homologous proteins
    genome_alignment_df = genome_alignment_df.merge(homologous_prots_df, on=['Seq1', 'Seq2'], how='left')

    # conserved proteins
    genome_alignment_df = genome_alignment_df.merge(conserved_prots_df, on=['Seq1', 'Seq2'], how='left')

    # BBH hits number
    genome_alignment_df = pd.merge(genome_alignment_df, cds_alignment_df, on=['Seq1', 'Seq2'], how='left')

    # minimum number of proteins
    genome_alignment_df['min_prots'] = np.round(genome_alignment_df[['seq1_n_prots', 'seq2_n_prots']].min(axis=1), 6)

    # calculate wgrr
    genome_alignment_df['cds_alignments_ani_sum'] = np.round(genome_alignment_df['ANI'] * genome_alignment_df['cds_alignments_counts'], 6)
    genome_alignment_df['wgrr'] = np.round(genome_alignment_df['cds_alignments_ani_sum'] / genome_alignment_df['min_prots'], 3)

    # save
    genome_alignment_df.to_csv(GENOME_ALIGNMENT, index=False)
else:
    merged.to_csv(GENOME_ALIGNMENT, index=False)


# remove intermediate files
if DELETE_INTERMEDIATE_FILES: 
    print('Removing intermediate files... ', end='')
    shutil.rmtree(INTERMEDIATE_FILES_DIR)
    print('Done!')

# remove mmseqs temp dir
try: shutil.rmtree(MMSEQS_TEMP_DIR)
except FileNotFoundError: pass

# end
print("\n\nSuccess! Thank you for using MANIAC!\n\n")