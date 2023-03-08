
import itertools
import os

from helpers import format_mmseqs_params


# Parameters
IN_DIR = config["input_dir"]
OUT_DIR = config["output_dir"]
TMP_DIR = config.get("tmp_dir", os.path.join(OUT_DIR, "tmp"))

FRAGMENT_SIZE = config.get("fragment_size", 1020)
INPUT_EXTENSION = config.get("input_extension", "fna")

MMSEQS_PARAMS = config.get("mmseqs_params",
    "--search-type 3 -a --max-seqs 10000"
)
MMSEQS_THREADS = config.get("threads", 8)

CDS_BASED_BBH = config.get("bbh", False)

# Find genomes

genomes, = glob_wildcards(os.path.join(IN_DIR, "{genome}."+INPUT_EXTENSION))

print("Running in all vs. all mode.\nGenomes:")
print(genomes)

print('------checkpoint------')


rule target:
    input: os.path.join(OUT_DIR, "search_results.tsv"),
           os.path.join(OUT_DIR, "best_hits.csv")
           

rule split_genomes:
    input: os.path.join(IN_DIR)
    output: os.path.join(OUT_DIR, "split/IN_DIR")
    params:
        fragment_size = FRAGMENT_SIZE
    script: "scripts/split_fasta.py"


rule make_query_db:
    input: os.path.join(OUT_DIR, "split/IN_DIR")
    output: os.path.join(OUT_DIR, "query-db/query-db")
    shell: "mmseqs createdb {input} {output}"

rule make_db:
    input: os.path.join(IN_DIR)
    output: os.path.join(OUT_DIR, "cds-db/cds-db") if CDS_BASED_BBH else \
			os.path.join(OUT_DIR, "reference-db/reference-db")
    shell: "mmseqs createdb {input} {output}"

rule get_total_cds_lengths:
    input: os.path.join(IN_DIR)
    output: os.path.join(OUT_DIR, "fasta_lengths.csv")
    script: "scripts/get_fasta_lengths.py"

# the following rule is carried out for Gorie et al ANI calculation
# as it searches fragmented pieces against the reference (whole genome) db
rule mmseqs_qr_search:
    input:
        os.path.join(OUT_DIR, "query-db/query-db"),
        os.path.join(OUT_DIR, "reference-db/reference-db")
    output:
        os.path.join(OUT_DIR, "results-qr-db/results-qr-db.index")
    params:
        results_basename = os.path.join(OUT_DIR, "results-qr-db/results-qr-db"),
        tmp_dir = TMP_DIR,
        mmseqs_params = MMSEQS_PARAMS
    threads:
        MMSEQS_THREADS
    shell:
        """
        mmseqs search {input} {params.results_basename} {params.tmp_dir} \
        --threads {threads} {params.mmseqs_params}
        """

# the following rule is carried out for the CDS-based BBH ANI calculation
# and it carries out an all-by-all CDS search
rule mmseqs_cds_search:
    input:
        os.path.join(OUT_DIR, "cds-db/cds-db")
    output:
        os.path.join(OUT_DIR, "results-cds-db/results-cds-db.index")
    params:
        results_basename = os.path.join(OUT_DIR, "results-cds-db/results-cds-db"),
        tmp_dir = TMP_DIR,
        mmseqs_params = MMSEQS_PARAMS
    threads:
        MMSEQS_THREADS
    shell:
        """
        mmseqs search {input} {input} {params.results_basename} {params.tmp_dir} \
        --threads {threads} {params.mmseqs_params}
        """

rule mmseqs_cds_convert:
    input:
        os.path.join(OUT_DIR, "cds-db/cds-db"),
        os.path.join(OUT_DIR, "results-cds-db/results-cds-db.index")
    output:
        os.path.join(OUT_DIR, "search_results.tsv")
    params:
        results_basename = os.path.join(OUT_DIR, "results-cds-db/results-cds-db")
    threads:
        MMSEQS_THREADS
    shell:
        """
        mmseqs convertalis {input[0]} {input[0]} \
        {params.results_basename} {output} --threads {threads} \
        --format-output 'qset,query,tset,target,nident,alnlen,mismatch,pident,evalue,qlen,gapopen,qstart,qend,tstart,tend,bits'
        """

rule mmseqs_qr_convert:
    input:
        os.path.join(OUT_DIR, "query-db/query-db"),
        os.path.join(OUT_DIR, "reference-db/reference-db"),
        os.path.join(OUT_DIR, "results-qr-db/results-qr-db.index")
    output:
        os.path.join(OUT_DIR, "search_results.tsv")
    params:
        results_basename = os.path.join(OUT_DIR, "results-qr-db/results-qr-db")
    threads:
        MMSEQS_THREADS
    shell:
        """
        mmseqs convertalis {input[0]} {input[1]} \
        {params.results_basename} {output} --threads {threads} \
        --format-output 'qset,query,tset,target,nident,alnlen,mismatch,pident,evalue,qlen,gapopen,qstart,qend,tstart,tend,bits'
        """

rule process_results:
    input:
        os.path.join(OUT_DIR, "search_results.tsv"),
        os.path.join(OUT_DIR, "fasta_lengths.csv")
    output:
        os.path.join(OUT_DIR, "best_hits.csv"),
        os.path.join(OUT_DIR, "ani.csv")
    params: eval_threshold = config.get("eval_filter", 1.0E-15),
            coverage_threshold = config.get("coverage_filter", 0.7),
            identity_threshold = config.get("identity_filter", 0.3),
            bbh_calc = config.get("bbh", False),
            memory_mode = config.get("Low_memory_mode", False),
            input_extension=INPUT_EXTENSION
    script: "scripts/process_results.py"

