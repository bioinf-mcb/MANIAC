import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
import gc



try:
    INPUT_PATH = snakemake.input[0]
    BEST_HITS_PATH = snakemake.output[0]
    
    EVAL_THR = snakemake.params["eval_threshold"]
    COVERAGE_THR = snakemake.params["coverage_threshold"]
    IDENTITY_THR = snakemake.params["identity_threshold"]
    BBH = snakemake.params["bbh_calc"]
    INPUT_EXTENSION = snakemake.params["input_extension"]
    Low_memory_mode = snakemake.params["memory_mode"]
except NameError:
    import argparse

    print('Warning! Script modified and not tested in this mode.')
    parser = argparse.ArgumentParser(description='Find best bidirectional hits from formatted MMSeqs2 results')
    parser.add_argument("input", help="Formatted MMSeqs2 results")
    parser.add_argument("fastalen", help="CSV file of fasta lengths")
    parser.add_argument("besthits", help="Output path for best hits")
    parser.add_argument("output", help="Output path for search results")
    parser.add_argument("--eval", default=10, type=float, help="E-value threshold (default: 10)")
    parser.add_argument("--coverage", default=0, type=float, help="Coverage threshold (default: 0)")
    parser.add_argument("--identity", default=0, type=float, help="Identity threshold (default: 0)")
    parser.add_argument("--bbh", default=False, action='store_true', help="Perform a CDS-based BBH ANI calculation")
    parser.add_argument("--input_extension", default=False, action='store_true', help="Extesion of input files")
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
    Low_memory_mode =args.Low_memory_mode

if (Low_memory_mode):
    print('low memory mode,some columns will not be loaded to save memomory')
    col_dtypes = {0:'object',1:'object',2:'object',3:'object',4:'int32',5:'int32',6:'int32',7:'float64',8:'float64',9:'int32'}
    RESULTS_HEADER = ["query_seq", "query_fragment_id", "reference_seq", "reference_fragment_id",
                  "matches", "length", "mismatches", "pident", "evalue", "qlen"]
                  
    print("Loading input files...")
    mmseqs_results = pd.read_csv(INPUT_PATH, sep = "\t", header = None,usecols=[0,1,2,3,4,5,6,7,8,9],dtype=col_dtypes)
    mmseqs_results.columns=RESULTS_HEADER


else:
    RESULTS_HEADER = ["query_seq", "query_fragment_id", "reference_seq", "reference_fragment_id",
                  "matches", "length", "mismatches", "pident", "evalue", "qlen",
                 "gapopen", "qstart", "qend", "rstart", "rend", "bitscore"]

    print("Loading input files...")
    mmseqs_results = pd.read_csv(INPUT_PATH, sep = "\t", header = None, names = RESULTS_HEADER)

print('dataframe loaded')
mmseqs_results.iloc[:,0]=mmseqs_results.iloc[:,1].str.split('.').str[0]
mmseqs_results.iloc[:,2]=mmseqs_results.iloc[:,3].str.split('.').str[0]


print("Filtering...")

mmseqs_results = mmseqs_results[mmseqs_results.query_seq != mmseqs_results.reference_seq].copy() # We don't want ANI with itself

mmseqs_results['gaps'] = mmseqs_results.length - mmseqs_results.matches - mmseqs_results.mismatches
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

del mmseqs_results
gc.collect()

hits = mmseqs_results_filtered.groupby(['query_fragment_id', 'reference_seq'], as_index=False).evalue.min()

print('groupby done')

# Only take the best hit for each fragment-reference pair
print("Selecting best hits...")
best_hits = pd.merge(mmseqs_results_filtered, hits, how='inner', on=['query_fragment_id', 'reference_seq', 'evalue'])
best_hits = best_hits.groupby(['query_fragment_id', 'reference_seq'], as_index=False).first()

del hits
del mmseqs_results_filtered
gc.collect()

best_hits.to_csv(BEST_HITS_PATH, index=False)




