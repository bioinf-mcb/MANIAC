# modules
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
import gc


# paths
INPUT_PATH = snakemake.input[0]
BEST_HITS_PATH = snakemake.output[0]

# params
EVAL_THR = snakemake.params.EVALUE
COVERAGE_THR = snakemake.params.COVERAGE
IDENTITY_THR = snakemake.params.IDENTITY
CDS_BASED = snakemake.params.CDS_BASED
MEMORY_EFFICIENT = snakemake.params.MEMORY_EFFICIENT


# loading params
header = ["query_seq", "query_fragment_id", "reference_seq", "reference_fragment_id",
          "matches", "length", "mismatches", "pident", "evalue", "qlen",
          "gapopen", "qstart", "qend", "rstart", "rend", "bitscore"]

usecols = list(range(0,16))

dtypes = {0:'object',1:'object', 2:'object',
          3:'object',4:'int32', 5:'int32',
          6:'int32', 7: 'float64', 8:'float64',
          9:'int32', 10: 'int32', 11: 'int32', 
          12: 'int32', 13: 'int32', 14: 'int32', 15: 'int32'}


if MEMORY_EFFICIENT:
    print("WARNING: Memory efficient mode ON! (dropping some columns)... ", end='')

    # load only nine first columns
    header = header[0:10]
    usecols = usecols[0:10]

print("Loading mmseqs results table... ", end='')
mmseqs_results = pd.read_csv(INPUT_PATH, sep = "\t", header = None, usecols=usecols, dtype=dtypes)
mmseqs_results.columns = header
print("Done!")


# curate phage IDs
print("Curate phage indetifiers... ", end='')

# queries phage identifiers [from fragments/ORFs/proteins]
query_fragment_id_list = mmseqs_results['query_fragment_id'].to_list()
query_seq = ['_'.join(query_fragment_id.split('_')[:-2]) for query_fragment_id in query_fragment_id_list]
mmseqs_results['query_seq'] = query_seq

# reference phage indetifiers
if CDS_BASED: 
    # fragments/ORFs/proteins
    reference_fragment_id_list = mmseqs_results['reference_fragment_id'].to_list()
    reference_seq = ['_'.join(reference_fragment_id.split('_')[:-2]) for reference_fragment_id in reference_fragment_id_list]
    mmseqs_results['reference_seq'] = reference_seq
else: 
    # map phageIDs
    mmseqs_results['reference_seq'] = mmseqs_results['reference_fragment_id']

print("Done!")

print("Filtering...", end='')

mmseqs_results = mmseqs_results[mmseqs_results.query_seq != mmseqs_results.reference_seq].copy() # We don't want ANI with itself

mmseqs_results["gaps"] = mmseqs_results.length - mmseqs_results.matches - mmseqs_results.mismatches
mmseqs_results["ani_alnlen"] = mmseqs_results.mismatches + mmseqs_results.matches
mmseqs_results["ani_ids"] = mmseqs_results.matches
mmseqs_results["ani_cov"] = mmseqs_results.ani_alnlen / mmseqs_results.qlen
mmseqs_results["ani_pid"] = mmseqs_results.ani_ids / mmseqs_results.qlen
mmseqs_results.pident = mmseqs_results.pident*0.01

mmseqs_results_filtered = mmseqs_results[
    (mmseqs_results.ani_cov > float(COVERAGE_THR)) &
    (mmseqs_results.ani_pid > float(IDENTITY_THR)) &
    (mmseqs_results.evalue < float(EVAL_THR))
].reset_index(drop=True)

# RAM release
del mmseqs_results
gc.collect()

print('grouping by query fragment ID... ', end='')
hits = mmseqs_results_filtered.groupby(['query_fragment_id', 'reference_seq'], as_index=False).evalue.min()
print('Done!')


# only take the best hit for each fragment-reference pair
print("Selecting best hits... ", end='')
best_hits = pd.merge(mmseqs_results_filtered, hits, how='inner', on=['query_fragment_id', 'reference_seq', 'evalue'])
best_hits = best_hits.groupby(['query_fragment_id', 'reference_seq'], as_index=False).first()
print("Done!")

# RAM release
del hits
del mmseqs_results_filtered
gc.collect()

print('Saving... ', end='')
best_hits.to_csv(BEST_HITS_PATH, index=False)
print('Done!')



