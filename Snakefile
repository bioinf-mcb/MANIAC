
import itertools
import os

from helpers import format_mmseqs_params

IN_DIR = config["input_dir"]
OUT_DIR = config["output_dir"]
TMP_DIR = config.get("tmp_dir", os.path.join(OUT_DIR, "tmp"))

FRAGMENT_SIZE = config.get("fragment_size", 1020)
INPUT_EXTENSION = config.get("input_extension", "fna")

MMSEQS_PARAMS = config.get("mmseqs_params",
    {}
)
FORMATTED_MMSEQS_PARAMS = format_mmseqs_params(MMSEQS_PARAMS)
MMSEQS_THREADS = config.get("threads", 8)


genomes, = glob_wildcards(os.path.join(IN_DIR, "{genome}."+INPUT_EXTENSION))
print("Running in all vs. all mode.\nGenomes:")
print(genomes)


rule target:
    input: os.path.join(OUT_DIR, "anib.csv")


rule make_reference_db:
    input: expand(os.path.join(IN_DIR, "{genome}."+INPUT_EXTENSION), genome=genomes)
    output: os.path.join(OUT_DIR, "reference-db/reference-db")
    shell: "mmseqs createdb {input} {output}" 


rule split_query_genome:
    input: os.path.join(IN_DIR, "{query_genome}."+INPUT_EXTENSION)
    output: os.path.join(OUT_DIR, "split/{query_genome}.fasta")
    params:
        fragment_size = FRAGMENT_SIZE
    script: "scripts/split_fasta.py"


rule get_query_lengths:
    input: expand(os.path.join(IN_DIR, "{genome}."+INPUT_EXTENSION), genome = genomes)
    output: os.path.join(OUT_DIR, "query_lengths.csv")
    script: "scripts/get_fasta_lengths.py"


rule make_query_db:
    input: expand(os.path.join(OUT_DIR, "split/{genome}.fasta"), genome=genomes)
    output: os.path.join(OUT_DIR, "query-db/query-db")
    shell: "mmseqs createdb {input} {output}" 


rule mmseqs_search:
    input:
        os.path.join(OUT_DIR, "query-db/query-db"),
        os.path.join(OUT_DIR, "reference-db/reference-db")
    output:
        os.path.join(OUT_DIR, "results-db/results-db.index")
    params:
        results_basename = os.path.join(OUT_DIR, "results-db/results-db"),
        tmp_dir = TMP_DIR,
        mmseqs_params = FORMATTED_MMSEQS_PARAMS
    threads: 
        MMSEQS_THREADS
    shell:
        """
        mmseqs search {input} {params.results_basename} {params.tmp_dir} \
        --search-type 3 -a --max-seqs 10000 --threads {threads} \
        {params.mmseqs_params}
        """


rule mmseqs_convert_results:
    input: 
        os.path.join(OUT_DIR, "query-db/query-db"),
        os.path.join(OUT_DIR, "reference-db/reference-db"),
        os.path.join(OUT_DIR, "results-db/results-db.index")
    output: 
        os.path.join(OUT_DIR, "search_results.tsv")
    params:
        results_basename = os.path.join(OUT_DIR, "results-db/results-db")
    threads:
        MMSEQS_THREADS
    shell:
        """
        mmseqs convertalis {input[0]} {input[1]} \
        {params.results_basename} {output} --threads {threads} \
        --format-output 'qset,query,tset,nident,alnlen,mismatch,gapopen,evalue,qlen'
        """


rule process_results:
    input: os.path.join(OUT_DIR, "search_results.tsv")
    output: os.path.join(OUT_DIR, "anib.csv")
    params: eval_threshold = config.get("eval_filter", 1.0E-15),
            coverage_threshold = config.get("coverage_filter", 0.7),
            identity_threshold = config.get("identity_filter", 0.3)
    script: "scripts/process_results.py"

