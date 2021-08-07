
import pandas as pd

INPUT_PATH = snakemake.input[0]
OUTPUT_PATH = snakemake.output[0]

EVAL_THR = snakemake.params["eval_threshold"]
COVERAGE_THR = snakemake.params["coverage_threshold"]
IDENTITY_THR = snakemake.params["identity_threshold"]

# RESULTS_HEADER = ["query", "reference", "identity", "aln_len",
#                   "n_mismatches", "n_gaps", "q_start", "q_end",
#                   "r_start", "r_end", "eval", "bitscore"]

RESULTS_HEADER = ["query_seq", "fragment_id", "reference_seq", "n_ident", "aln_len",
                  "n_mismatches", "n_gaps", "evalue", "qlen"]

mmseqs_results = pd.read_csv(INPUT_PATH, sep = "\t", header = None,
                             names = RESULTS_HEADER)

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
# Very time-consuming --> Maybe split the table and process individually?
best_hits = mmseqs_results_filtered.iloc[
    mmseqs_results_filtered.groupby(['fragment_id', 'reference_seq']).aln_identity.idxmax()
].copy()

best_hits.query_seq = best_hits.query_seq.str.split(".", expand = True)[0]
best_hits.reference_seq = best_hits.reference_seq.str.split(".", expand = True)[0]

ani = best_hits.groupby(["query_seq", "reference_seq"]).aln_identity.mean().reset_index()
ani.to_csv(OUTPUT_PATH, index=False)