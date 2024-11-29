
##Get params
INPUT_PATH = snakemake@input[[1]]
PHAGE_LENGTHS_PATH = snakemake@input[[2]]

COVERAGE_THR = snakemake@params[["COVERAGE"]]
IDENTITY_THR = snakemake@params[["IDENTITY"]]
CDS_BASED = snakemake@params[["CDS_BASED"]]
SEPARATOR = snakemake@params[["SEPARATOR"]]
MODE = snakemake@params[["MODE"]]

GENOME_ALIGNMENT = snakemake@output[[1]]
CDS_ALIGNMENT_FILE = snakemake@params[["CDS_ALIGNMENT_FILE"]]
MMSEQS_TEMP_DIR = snakemake@params[["MMSEQS_TEMP_DIR"]]

CONSERVED_IDENTITY = snakemake@params[["CONSERVED_IDENTITY"]]
CONSERVED_COVERAGE = snakemake@params[["CONSERVED_COVERAGE"]]

DELETE_INTERMEDIATE_FILES = snakemake@params[["DELETE_INTERMEDIATE_FILES"]]
INTERMEDIATE_FILES_DIR = dirname(snakemake@input[[1]])

# Load necessary libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


# loading params
col.names <- c("query_fragment_id", "reference_fragment_id", "matches", "length", "mismatches", "pident", "evalue", "qlen")

# Two methods depending on the mode (fragment or CDS: fragment is processed by chunks while CDS have to be done all at once)

if(CDS_BASED){ #METHOD 1: CDS
	# Load mmseqs results
	writeLines("Loading mmseqs results table...")
	mmseqs_results <- fread(INPUT_PATH, header = F)
	setnames(mmseqs_results, col.names)

	# queries phage identifiers [from fragments/ORFs/proteins]
	writeLines("Curate phage identifiers...")
	mmseqs_results[, query_seq := sapply(strsplit(mmseqs_results[['query_fragment_id']], SEPARATOR, fixed = TRUE), `[`, 1)]
	mmseqs_results[, reference_seq := sapply(strsplit(mmseqs_results[['reference_fragment_id']], SEPARATOR, fixed = TRUE), `[`, 1)]

	writeLines("Filtering...")
	mmseqs_results[, gaps := length - matches - mismatches]
	mmseqs_results[, ani_alnlen := mismatches + matches]
	mmseqs_results[, ani_cov := ani_alnlen / qlen]
	mmseqs_results[, ani_pid := matches / qlen]
	mmseqs_results[, pident := pident * 0.01]

	mmseqs_results <- mmseqs_results[query_seq != reference_seq & ani_cov > COVERAGE_THR & ani_pid > IDENTITY_THR]

	writeLines("Finding best bidirectional hits (BBH)...")
	# Perform an inner join to find bidirectional hits
	mmseqs_results <- merge(mmseqs_results, mmseqs_results, 
	                         by.x = c("query_fragment_id", "reference_fragment_id"),
	                         by.y = c("reference_fragment_id", "query_fragment_id"))
	    
	# Deduplicate by filtering where query_seq_x < reference_seq_x
	mmseqs_results <- mmseqs_results[query_seq.x < reference_seq.x]
	  
	 # Calculate the average pident and ani_alnlen
	mmseqs_results[, `:=`(
	pident = round((pident.x + pident.y) / 2, 6),
	ani_alnlen = round((ani_alnlen.x + ani_alnlen.y) / 2, 6)
	)]

	# clean best hits table 
	setnames(mmseqs_results, old = c("query_seq.x", "reference_seq.x"), new = c("query_seq", "reference_seq"))

	writeLines("Calculating ANI...")
	# Load genome lengths
	lengths_df <- fread(PHAGE_LENGTHS_PATH)

	# genome lengths (ANImm calculation)
	genome_length_df <- lengths_df[, .(genome, length)]
	setkey(genome_length_df, genome)

	# calculate ANI
	ani <- mmseqs_results[, .(ANI = mean(pident)), by = .(query_seq, reference_seq)]

	# calculate different measures
	aligned_nucleotides <- mmseqs_results[, .(ani_alnlen = sum(ani_alnlen)), by = .(query_seq, reference_seq)]
	aligned_nucleotides[, len_1 := genome_length_df[query_seq, length]]
	aligned_nucleotides[, len_2 := genome_length_df[reference_seq, length]]
	aligned_nucleotides[, af_1 := round(ani_alnlen / len_1, 6)]
	aligned_nucleotides[, af_2 := round(ani_alnlen / len_2, 6)]
	aligned_nucleotides[, af_mean := round(2 * ani_alnlen / (len_1 + len_2), 6)]
	aligned_nucleotides[, af_min := round(ani_alnlen / pmin(len_1, len_2), 6)]
	aligned_nucleotides[, af_max := round(ani_alnlen / pmax(len_1, len_2), 6)]
	aligned_nucleotides[, af_jaccard := round(ani_alnlen / (len_1 + len_2 - ani_alnlen), 6)]
	# add measures and rename columns
	merged <- merge(ani, aligned_nucleotides, by = c("query_seq", "reference_seq"))
	merged[, wgANI := round(ANI * af_mean, 6)]
	setnames(merged, c("query_seq", "reference_seq", "ANI"), c("Seq1", "Seq2", "ANI"))

	writeLines("Calculating wGRR...")
	# Proteins number per genome (wGRR calculation)
	n_prot_df <- lengths_df[, .(Seq1 = genome, n_prots)]

	# Get number of BBH used to calculate mean ani
	cds_alignment_cols <- c("query_seq", "reference_seq", "query_fragment_id", "reference_fragment_id", "pident.x", "ani_cov.x", "pident.y", "ani_cov.y")
	cds_alignment_df <- mmseqs_results[, .N, by = cds_alignment_cols][, cds_alignment_cols, with=FALSE] # Group by phages and ORFs

	# Rename columns
	setnames(cds_alignment_df, old = c("query_seq", "reference_seq", "query_fragment_id", "reference_fragment_id", 
	                                   "pident.x", "ani_cov.x", "pident.y", "ani_cov.y"),
	                            new = c("Seq1", "Seq2", "seq1_fragment_id", "seq2_fragment_id", 
	                                    "seq1_fragment_pident", "seq1_fragment_cov", 
	                                    "seq2_fragment_pident", "seq2_fragment_cov"))

	# Calculate homologous proteins per phage
	cds_alignment_df[, cds_alignments_counts := 1]

	### seq1
	seq1_homologous_prots_df <- cds_alignment_df[, .(seq1_n_prots_hom = .N), by = .(Seq1, Seq2)]

	### seq2
	seq2_homologous_prots_df <- cds_alignment_df[, .(seq2_n_prots_hom = .N), by = .(Seq2, Seq1)]

	# Merge
	homologous_prots_df <- merge(seq1_homologous_prots_df, seq2_homologous_prots_df, by = c("Seq1", "Seq2"), all = TRUE)

	# Calculate conserved proteins per phage
	filt_conserved_identity <- (cds_alignment_df$seq1_fragment_pident >= CONSERVED_IDENTITY) & (cds_alignment_df$seq2_fragment_pident >= CONSERVED_IDENTITY)
	filt_conserved_coverage <- (cds_alignment_df$seq1_fragment_cov >= CONSERVED_COVERAGE) & (cds_alignment_df$seq2_fragment_cov >= CONSERVED_COVERAGE)
	conserved_df <- cds_alignment_df[filt_conserved_identity & filt_conserved_coverage]

	# Calculate conserved proteins per phage for Seq1
	seq1_conserved_prots_df <- conserved_df[, .(seq1_n_prots_cons = .N), by = .(Seq1, Seq2)]

	# Calculate conserved proteins per phage for Seq2
	seq2_conserved_prots_df <- conserved_df[, .(seq2_n_prots_cons = .N), by = .(Seq2, Seq1)]

	# Merge the results for Seq1 and Seq2
	conserved_prots_df <- merge(seq1_conserved_prots_df, seq2_conserved_prots_df, by = c("Seq1", "Seq2"), all = TRUE)

	# Specify columns to save
	cols2save <- c("Seq1", "Seq2", "seq1_fragment_id", "seq2_fragment_id", "seq1_fragment_pident", "seq2_fragment_pident", "seq1_fragment_cov", "seq2_fragment_cov")

	if (!DELETE_INTERMEDIATE_FILES) {
	    writeLines("Saving CDS alignment file... ")
	    fwrite(cds_alignment_df[, ..cols2save], CDS_ALIGNMENT_FILE, row.names = FALSE)
	}

	# Count BBH
	cds_alignment_df <- cds_alignment_df[, .N, by = .(Seq1, Seq2)]
	setnames(cds_alignment_df, "N", "cds_alignments_counts")

	### Add number of protein types to genome alignments

	# Total number of proteins
	genome_alignment_df <- merge(merged, n_prot_df, by = "Seq1", all.x = TRUE)
	setnames(genome_alignment_df, "n_prots", "seq1_n_prots")

	setnames(n_prot_df, "Seq1", "Seq2")
	genome_alignment_df <- merge(genome_alignment_df, n_prot_df, by = "Seq2", all.x = TRUE)
	setnames(genome_alignment_df, "n_prots", "seq2_n_prots")

	# Homologous proteins
	genome_alignment_df <- merge(genome_alignment_df, homologous_prots_df, by = c("Seq1", "Seq2"), all.x = TRUE)

	# Conserved proteins
	genome_alignment_df <- merge(genome_alignment_df, conserved_prots_df, by = c("Seq1", "Seq2"), all.x = TRUE)
	genome_alignment_df[, c("seq1_n_prots_cons", "seq2_n_prots_cons") := lapply(.SD, function(x) replace(x, is.na(x), 0)), .SDcols = c("seq1_n_prots_cons", "seq2_n_prots_cons")]
	# BBH hits number
	genome_alignment_df <- merge(genome_alignment_df, cds_alignment_df, by = c("Seq1", "Seq2"), all.x = TRUE)

	# Minimum number of proteins
	genome_alignment_df[, min_prots := pmin(seq1_n_prots, seq2_n_prots, na.rm = TRUE)]

	# Calculate wGRR
	genome_alignment_df[, cds_alignments_ani_sum := round(ANI * cds_alignments_counts, 6)]
	genome_alignment_df[, wgrr := round(cds_alignments_ani_sum / min_prots, 3)]

	if (MODE == "CDS_AA") {
		setnames(genome_alignment_df, old = c("ANI", "wgANI"), new = c("AAI", "wgAAI"))
	}
	# Save
	writeLines("Saving final table...")
	fwrite(genome_alignment_df, GENOME_ALIGNMENT, row.names = FALSE)

} else{ #METHOD 2: FRAGMENT
	writeLines("Starting post-processing...")
	files <- list.files(path=dirname(INPUT_PATH), pattern="chunk_.*.tsv", full.names=TRUE)
	SEPARATOR2 <- paste(SEPARATOR,".*",sep="")

	lengths_df <- fread(PHAGE_LENGTHS_PATH)
	genome_length_df <- lengths_df[, .(genome, length)]
	setkey(genome_length_df, genome)

	invisible(lapply(files, function(x) {
		mmseqs_results_csv <- open_dataset(sources = x,format='tsv', read_options = CsvReadOptions$create(column_names = col.names))
		mmseqs_results <- mmseqs_results_csv |> 
		  filter(sub(SEPARATOR2,"", query_fragment_id) != reference_fragment_id, (mismatches + matches) / qlen > COVERAGE_THR, matches / qlen > IDENTITY_THR) |>
		  group_by(Seq1=sub(SEPARATOR2,"", query_fragment_id),Seq2=reference_fragment_id) |>
		  summarize(ANI = mean(pident)* 0.01, ani_alnlen = sum(mismatches + matches)) |>
		  collect()
		mmseqs_results <- data.table(mmseqs_results)
		mmseqs_results[, len_1 := genome_length_df[Seq1, length]]
		mmseqs_results[, len_2 := genome_length_df[Seq2, length]]
		mmseqs_results[, af_1 := round(ani_alnlen / len_1, 6)]
		mmseqs_results[, af_2 := round(ani_alnlen / len_2, 6)]
		mmseqs_results[, af_mean := round(2 * ani_alnlen / (len_1 + len_2), 6)]
		mmseqs_results[, af_min := round(ani_alnlen / pmin(len_1, len_2), 6)]
		mmseqs_results[, af_max := round(ani_alnlen / pmax(len_1, len_2), 6)]
		mmseqs_results[, af_jaccard := round(ani_alnlen / (len_1 + len_2 - ani_alnlen), 6)]
		mmseqs_results[, wgANI := round(ANI * af_mean, 6)]
		mmseqs_results <- mmseqs_results[order(Seq1,Seq2)]
		fwrite(mmseqs_results, GENOME_ALIGNMENT, row.names = FALSE, append = TRUE)
		unlink(x)
	}))
}

# Remove intermediate files
if (DELETE_INTERMEDIATE_FILES) {
  writeLines('Removing intermediate files...')
  unlink(INTERMEDIATE_FILES_DIR, recursive = TRUE)
}

# Remove mmseqs temp dir
tryCatch({
  unlink(MMSEQS_TEMP_DIR, recursive = TRUE)
}, error = function(e) {
  if (!grepl("no temporary mmseqs dir found", e$message, ignore.case = TRUE)) stop(e)
})

# End
writeLines("\n\nSuccess! Thank you for using MANIAC!\n\n")
