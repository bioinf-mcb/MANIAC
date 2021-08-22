import pandas as pd

BEST_HITS_PATH = snakemake.input[0]
BBH_PATH = snakemake.input[1]
LENGTHS_PATH = snakemake.input[2]

OUTPUT_PATH = snakemake.output[0]

print("Loading best bidirectional hits...")
bbh = pd.read_csv(BBH_PATH)

print("Formatting...")
bbh.query_seq_x = bbh.query_seq_x.str.split(".", expand = True)[0]
bbh.reference_seq_x = bbh.reference_seq_x.str.split(".", expand = True)[0]


print("Computing ANI...")
ani = bbh.groupby(["query_seq_x", "reference_seq_x"]).aln_identity.mean().reset_index()

print("Computing coverage...")

genome_lengths = pd.read_csv(LENGTHS_PATH, index_col=0)

aligned_nucleotides = bbh.groupby(["query_seq_x", "reference_seq_x"]).aligned_pairs.sum().reset_index()
aligned_nucleotides["len_1"] = genome_lengths.loc[aligned_nucleotides.query_seq_x].reset_index(drop=True)
aligned_nucleotides["len_2"] = genome_lengths.loc[aligned_nucleotides.reference_seq_x].reset_index(drop=True)
aligned_nucleotides["af_1"] = aligned_nucleotides.aligned_pairs / aligned_nucleotides.len_1
aligned_nucleotides["af_2"] = aligned_nucleotides.aligned_pairs / aligned_nucleotides.len_2
aligned_nucleotides["af_mean"] = 2* aligned_nucleotides.aligned_pairs / (aligned_nucleotides.len_1 + aligned_nucleotides.len_2)
aligned_nucleotides["af_min"] = aligned_nucleotides.aligned_pairs / aligned_nucleotides[["len_1", "len_2"]].min(axis=1)
aligned_nucleotides["af_jaccard"] = aligned_nucleotides.aligned_pairs / (aligned_nucleotides.len_1 + aligned_nucleotides.len_2 - aligned_nucleotides.aligned_pairs)

merged = pd.merge(ani, aligned_nucleotides, on = ["query_seq_x", "reference_seq_x"]) \
           .rename({"query_seq_x": "Seq1", "reference_seq_x": "Seq2", "aln_identity": "ANI"}, axis=1)

merged.to_csv(OUTPUT_PATH, index=False)