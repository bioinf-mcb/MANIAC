#!/bin/bash

MEMORY_THIRD=$((snakemake_params[MEMORY_GB] * 1 / 3))
if [ "${snakemake_params[MMSEQS_THREADS]}" -le "$MEMORY_THIRD" ]; then
  JOBNB=${snakemake_params[MMSEQS_THREADS]}
else
  JOBNB=$MEMORY_THIRD
fi

split -l ${snakemake_params[CHUNKSIZE]} "${snakemake_input[0]}" "${snakemake_params[CHUNK]}"

ls "${snakemake_params[CHUNK]}"* | parallel --will-cite --silent -j$JOBNB Rscript ${snakemake_params[SORTPATH]} ${snakemake_params[SEPARATOR]} {} 1> /dev/null

sort --parallel=${snakemake_params[MMSEQS_THREADS]} -m -k1,1 -k9,9 -k7,7g "${snakemake_params[INTERMEDIATE_FILES_DIR]}"/sorted_* > "${snakemake_params[MERGED]}"

rm "${snakemake_params[CHUNK]}"* "${snakemake_params[INTERMEDIATE_FILES_DIR]}"/sorted_*

datamash -g 1,9 first 2 first 3 first 4 first 5 first 6 first 7 first 8 < "${snakemake_params[MERGED]}" | cut -f1,3-9 > "${snakemake_output[0]}"

rm "${snakemake_params[MERGED]}"

awk -v chunk_size=${snakemake_params[CHUNKSIZE]} -v prefix="${snakemake_params[CHUNK]}" -v separator="${snakemake_params[SEPARATOR]}" '
BEGIN { 
	current_id = ""; 
	chunk_count = 0; 
	line_count = 0; 
}
{
	split($1, id_parts, separator);
	id = id_parts[1];
	if (line_count >= chunk_size && id != current_id) {
		close(output_file);
		chunk_count++;
		line_count = 0;
	}
	if (line_count == 0) {
		output_file = sprintf("%s%04d.tsv", prefix, chunk_count);
	}
	print $0 >> output_file;
	line_count++;
	current_id = id;
}' "${snakemake_output[0]}"
