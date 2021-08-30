
import pandas as pd

try:
    INPUT_PATH = snakemake.input[0]
    BEST_HITS_PATH = snakemake.output[0]
    BBH_PATH = snakemake.output[1]

    EVAL_THR = snakemake.params["eval_threshold"]
    COVERAGE_THR = snakemake.params["coverage_threshold"]
    IDENTITY_THR = snakemake.params["identity_threshold"]

except NameError:
    import argparse

    parser = argparse.ArgumentParser(description='Find best bidirectional hits from formatted MMSeqs2 results')
    parser.add_argument("input", help="Formatted MMSeqs2 results")
    parser.add_argument("besthits", help="Output path for best hits")
    parser.add_argument("bbh", help="Output path for best bidirectional hits")
    parser.add_argument("--eval", default=10, type=float, help="E-value threshold (default: 10)")
    parser.add_argument("--coverage", default=0, type=float, help="Coverage threshold (default: 0)")
    parser.add_argument("--identity", default=0, type=float, help="Identity threshold (default: 0)")

    args = parser.parse_args()

    INPUT_PATH = args.input
    BEST_HITS_PATH = args.besthits
    BBH_PATH = args.bbh

    EVAL_THR = args.eval
    COVERAGE_THR = args.coverage
    IDENTITY_THR = args.identity


RESULTS_HEADER = ["query_seq", "query_fragment_id", "reference_seq", "reference_fragment_id", "n_ident", "aln_len",
                  "n_mismatches", "n_gaps", "evalue", "qlen"]


print("Loading results...")
mmseqs_results = pd.read_csv(INPUT_PATH, sep = "\t", header = None,
                             names = RESULTS_HEADER)

print("Filtering...")

mmseqs_results = mmseqs_results[mmseqs_results.query_seq != mmseqs_results.reference_seq] # We don't want ANI with itself

mmseqs_results["aligned_pairs"] = mmseqs_results.n_ident + mmseqs_results.n_mismatches
mmseqs_results["query_coverage"] = mmseqs_results.aligned_pairs / mmseqs_results.qlen
mmseqs_results["query_identity"] = mmseqs_results.n_ident / mmseqs_results.qlen
mmseqs_results["aln_identity"] = mmseqs_results.n_ident / mmseqs_results.aligned_pairs

mmseqs_results_filtered = mmseqs_results[
    (mmseqs_results.query_coverage > COVERAGE_THR) &
    (mmseqs_results.query_identity > IDENTITY_THR) &
    (mmseqs_results.evalue < EVAL_THR)
].reset_index(drop=True)

# Only take the best hit for each fragment-reference pair
print("Selecting best hits...")
best_hits = mmseqs_results_filtered.iloc[
    mmseqs_results_filtered.groupby(['query_fragment_id', 'reference_seq']).evalue.idxmin()
].copy()

best_hits.to_csv(BEST_HITS_PATH, index=False)

print("Finding best bidirectional hits...")
bbh = pd.merge(best_hits, best_hits, how="inner", \
      left_on=["query_fragment_id", "reference_fragment_id"], \
      right_on=["reference_fragment_id", "query_fragment_id"])

bbh = bbh[bbh.query_seq_x < bbh.reference_seq_x].copy() # Avoid duplication of each pair

bbh["aln_identity"] = (bbh.aln_identity_x + bbh.aln_identity_y) / 2
bbh["aligned_pairs"] = (bbh.aligned_pairs_x + bbh.aligned_pairs_y) / 2

bbh.to_csv(BBH_PATH, index=False)

print("Done!")
