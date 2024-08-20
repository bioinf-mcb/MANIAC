##Add the conda install r-base and r-essentials with versions (4.4.1 and 4.4)

##Get params
INPUT_PATH = snakemake@input[[1]]
PHAGE_LENGTHS_PATH = snakemake@input[[2]]

BEST_HITS_PATH = snakemake@params[["BEST_HITS_PATH"]]
COVERAGE_THR = snakemake@params[["COVERAGE"]]
IDENTITY_THR = snakemake@params[["IDENTITY"]]
CDS_BASED = snakemake@params[["CDS_BASED"]]
SEPARATOR = snakemake@params[["SEPARATOR"]]

MEMORY_EFFICIENT = snakemake@params[["MEMORY_EFFICIENT"]]
GENOME_ALIGNMENT = snakemake@output[[1]]
CDS_ALIGNMENT_FILE = snakemake@params[["CDS_ALIGNMENT_FILE"]]
MMSEQS_TEMP_DIR = snakemake@params[["MMSEQS_TEMP_DIR"]]

CONSERVED_IDENTITY = snakemake@params[["CONSERVED_IDENTITY"]]
CONSERVED_COVERAGE = snakemake@params[["CONSERVED_COVERAGE"]]

DELETE_CDS_ALIGNMENT = snakemake@params[["DELETE_CDS_ALIGNMENT"]]
DELETE_INTERMEDIATE_FILES = snakemake@params[["DELETE_INTERMEDIATE_FILES"]]
INTERMEDIATE_FILES_DIR = dirname(snakemake@input[[1]])

# Load necessary libraries
suppressPackageStartupMessages(library(data.table))

# loading params
col.names <- c("query_seq", "query_fragment_id", "reference_seq", "reference_fragment_id",
               "matches", "length", "mismatches", "pident", "evalue", "qlen",
               "gapopen", "qstart", "qend", "rstart", "rend", "bitscore")
cols <- 1:16

if(MEMORY_EFFICIENT){
  writeLines("WARNING: Memory efficient mode ON! (dropping some columns)...")
  cols <- cols[1:10]
  col.names <- col.names[1:10]
}

# Load mmseqs results
writeLines("Loading mmseqs results table...")
mmseqs_results <- fread(INPUT_PATH, select = cols, header = F)
setnames(mmseqs_results, col.names)

# queries phage identifiers [from fragments/ORFs/proteins]
writeLines("Curate phage identifiers...")
mmseqs_results[, separator_pos := regexpr(SEPARATOR, query_fragment_id, fixed = TRUE)]
mmseqs_results[, query_seq := ifelse(separator_pos > 0, substr(query_fragment_id, 1, separator_pos - 1), query_fragment_id)]
mmseqs_results[, separator_pos := NULL]

# reference phage identifiers
if(CDS_BASED){
	# fragments/ORFs/proteins
    mmseqs_results[, reference_seq := sapply(strsplit(mmseqs_results[['reference_fragment_id']], SEPARATOR, fixed = TRUE), `[`, 1)]
} else{
	# map phageIDs
	mmseqs_results[, reference_seq := reference_fragment_id]
}

# Filtering
writeLines("Filtering...")
mmseqs_results[, gaps := length - matches - mismatches]
mmseqs_results[, ani_alnlen := mismatches + matches]
mmseqs_results[, ani_ids := matches]
mmseqs_results[, ani_cov := ani_alnlen / qlen]
mmseqs_results[, ani_pid := ani_ids / qlen]
mmseqs_results[, pident := pident * 0.01]

mmseqs_results <- mmseqs_results[query_seq != reference_seq & ani_cov > COVERAGE_THR & ani_pid > IDENTITY_THR]

# grouping by query fragment ID
writeLines("Grouping by query fragment ID...")
mmseqs_results <- mmseqs_results[order(evalue)]

# only take the best hit for each fragment-reference pair
writeLines("Selecting best hits...")
best_hits <- mmseqs_results[, .SD[1], by = .(query_fragment_id, reference_seq)]

# Saving
writeLines("Saving best hits...")
fwrite(best_hits, BEST_HITS_PATH, row.names = FALSE)

# Clean up
rm(mmseqs_results)

# If CDS calculate BBH
if (CDS_BASED) {
  writeLines("Finding best bidirectional hits (BBH)...")
  
  # Perform an inner join to find bidirectional hits
  best_hits_final <- merge(best_hits, best_hits, 
                           by.x = c("query_fragment_id", "reference_fragment_id"),
                           by.y = c("reference_fragment_id", "query_fragment_id"))
    
  # Deduplicate by filtering where query_seq_x < reference_seq_x
  best_hits_final <- best_hits_final[query_seq.x < reference_seq.x]
  
  # Calculate the average pident and ani_alnlen
  best_hits_final[, `:=`(
    pident = round((pident.x + pident.y) / 2, 6),
    ani_alnlen = round((ani_alnlen.x + ani_alnlen.y) / 2, 6)
  )]

} else {
  best_hits_final <- best_hits
}

# Remove the original dataframe to free memory
rm(best_hits)

# clean best hits table 
setnames(best_hits_final, old = c("query_seq.x", "reference_seq.x"), new = c("query_seq", "reference_seq"))

writeLines("Calculating ANI...")
# Load genome lengths
lengths_df <- fread(PHAGE_LENGTHS_PATH)

# genome lengths (ANImm calculation)
genome_length_df <- lengths_df[, .(genome, length)]
setkey(genome_length_df, genome)

# calculate ANI
ani <- best_hits_final[, .(ANI = mean(pident)), by = .(query_seq, reference_seq)]

# calculate different measures
aligned_nucleotides <- best_hits_final[, .(ani_alnlen = sum(ani_alnlen)), by = .(query_seq, reference_seq)]
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
setnames(merged, c("query_seq", "reference_seq", "ANI"), c("Seq1", "Seq2", "ANI"))

# save
if (CDS_BASED) {
	writeLines("Calculating wGRR...")
	# Proteins number per genome (wGRR calculation)
	n_prot_df <- lengths_df[, .(Seq1 = genome, n_prots)]

	# Get number of BBH used to calculate mean ani
	cds_alignment_cols <- c("query_seq", "reference_seq", "query_fragment_id", "reference_fragment_id", "pident.x", "ani_cov.x", "pident.y", "ani_cov.y")
	cds_alignment_df <- best_hits_final[, .N, by = cds_alignment_cols][, cds_alignment_cols, with=FALSE] # Group by phages and ORFs

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

	if (!DELETE_CDS_ALIGNMENT) {
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

	# Save
	writeLines("Saving final table...")
	fwrite(genome_alignment_df, GENOME_ALIGNMENT, row.names = FALSE)
	
} else{
	writeLines("Saving final table...")
	fwrite(merged, GENOME_ALIGNMENT, row.names = FALSE)
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
