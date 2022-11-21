import pandas as pd

try:
    INPUT_PATH = snakemake.input[0]
    LENGTHS_PATH = snakemake.input[1]
    BEST_HITS_PATH = snakemake.output[0]
    OUTPUT_PATH = snakemake.output[1]

    EVAL_THR = snakemake.params["eval_threshold"]
    COVERAGE_THR = snakemake.params["coverage_threshold"]
    IDENTITY_THR = snakemake.params["identity_threshold"]
    BBH = snakemake.params["bbh_calc"]
    INPUT_EXTENSION = snakemake.params["input_extension"]

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


RESULTS_HEADER = ["query_seq", "query_fragment_id", "reference_seq", "reference_fragment_id",
                  "matches", "length", "mismatches", "pident", "evalue", "qlen",
                 "gapopen", "qstart", "qend", "rstart", "rend", "bitscore"]

print("Loading input files...")
mmseqs_results = pd.read_csv(INPUT_PATH, sep = "\t", header = None, names = RESULTS_HEADER)
mmseqs_results.iloc[:,0]=mmseqs_results.iloc[:,1].str.split('.').str[0]
mmseqs_results.iloc[:,2]=mmseqs_results.iloc[:,3].str.split('.').str[0]
fasta_lengths = pd.read_csv(LENGTHS_PATH, index_col=0)


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

# Only take the best hit for each fragment-reference pair
print("Selecting best hits...")
print("ORFs headers have to be unique for each phage!")

best_hits = mmseqs_results_filtered.iloc[
    mmseqs_results_filtered.groupby(['query_fragment_id', 'reference_seq']).evalue.idxmin()
].copy()

if (BBH):
    print("Finding best bidirectional hits...")
    bbh = pd.merge(best_hits, best_hits, how="inner", \
      left_on=["query_fragment_id", "reference_fragment_id"], \
      right_on=["reference_fragment_id", "query_fragment_id"])
    bbh = bbh[bbh.query_seq_x < bbh.reference_seq_x].copy() # Avoid duplication of each pair
    bbh["pident"] = (bbh.pident_x + bbh.pident_y) / 2
    bbh["ani_alnlen"] = (bbh.ani_alnlen_x + bbh.ani_alnlen_y) / 2
    best_hits_final = bbh
else:
    best_hits_final = best_hits

best_hits_final = best_hits_final.rename(columns = {"query_seq_x": "query_seq", "reference_seq_x": "reference_seq"})

try:
    best_hits_final.query_seq = best_hits_final.query_seq.str.split(f".{INPUT_EXTENSION}", expand = True)[0]
    best_hits_final.reference_seq = best_hits_final.reference_seq.str.split(f".{INPUT_EXTENSION}", expand = True)[0]
except:
    print('No bbh hits were found! Throws and error and stops snakemake :/ ')

best_hits_final.to_csv(BEST_HITS_PATH, index=False)

ani_mean = best_hits_final.groupby(["query_seq", "reference_seq"]).pident.mean().reset_index()
ani_mean.rename({'pident':'ani_mean'},axis=1,inplace=True)


aligned_nucleotides = best_hits_final.groupby(["query_seq", "reference_seq"]).ani_alnlen.sum().reset_index()
aligned_nucleotides["len_1"] = fasta_lengths.loc[aligned_nucleotides.query_seq].reset_index(drop=True)
aligned_nucleotides["len_2"] = fasta_lengths.loc[aligned_nucleotides.reference_seq].reset_index(drop=True)
aligned_nucleotides["af_1"] = aligned_nucleotides.ani_alnlen / aligned_nucleotides.len_1
aligned_nucleotides["af_2"] = aligned_nucleotides.ani_alnlen / aligned_nucleotides.len_2
aligned_nucleotides["af_mean"] = 2* aligned_nucleotides.ani_alnlen / (aligned_nucleotides.len_1 + aligned_nucleotides.len_2)
aligned_nucleotides["af_min"] = aligned_nucleotides.ani_alnlen / aligned_nucleotides[["len_1", "len_2"]].min(axis=1)
aligned_nucleotides["af_jaccard"] = aligned_nucleotides.ani_alnlen / (aligned_nucleotides.len_1 + aligned_nucleotides.len_2 - aligned_nucleotides.ani_alnlen)

merged = pd.merge(ani_mean, aligned_nucleotides, on = ["query_seq", "reference_seq"]) \
           .rename({"query_seq": "Seq1", "reference_seq": "Seq2"}, axis=1)

merged.to_csv(OUTPUT_PATH, index=False)

print("Done!")
