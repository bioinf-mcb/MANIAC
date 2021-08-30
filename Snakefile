
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

SKIP_FRAGMENTATION = config.get("prefragmented", False)

# Find genomes

genomes, = glob_wildcards(os.path.join(IN_DIR, "{genome}."+INPUT_EXTENSION))

print("Running in all vs. all mode.\nGenomes:")
print(genomes)


rule target:
    input: os.path.join(OUT_DIR, "ani.csv")


rule split_genomes:
    input: os.path.join(IN_DIR, "{genome}."+INPUT_EXTENSION)
    output: os.path.join(OUT_DIR, "split/{genome}.fasta")
    params:
        fragment_size = FRAGMENT_SIZE
    script: "scripts/split_fasta.py"


rule get_genome_lengths:
    input: expand(os.path.join(IN_DIR, "{genome}."+INPUT_EXTENSION), genome = genomes)
    output: os.path.join(OUT_DIR, "genome_lengths.csv")
    script: "scripts/get_fasta_lengths.py"


rule make_db:
    input: 
        expand(os.path.join(IN_DIR, "{genome}."+INPUT_EXTENSION), genome=genomes) \
        if SKIP_FRAGMENTATION else \
        expand(os.path.join(OUT_DIR, "split/{genome}.fasta"), genome=genomes)
    output: os.path.join(OUT_DIR, "genome-db/genome-db")
    shell: "mmseqs createdb {input} {output}" 


rule mmseqs_search:
    input:
        os.path.join(OUT_DIR, "genome-db/genome-db"),
    output:
        os.path.join(OUT_DIR, "results-db/results-db.index")
    params:
        results_basename = os.path.join(OUT_DIR, "results-db/results-db"),
        tmp_dir = TMP_DIR,
        mmseqs_params = MMSEQS_PARAMS
    threads: 
        MMSEQS_THREADS
    shell:
        """
        mmseqs search {input} {input} {params.results_basename} {params.tmp_dir} \
        --threads {threads} {params.mmseqs_params}
        """


rule mmseqs_convert_results:
    input: 
        os.path.join(OUT_DIR, "genome-db/genome-db"),
        os.path.join(OUT_DIR, "results-db/results-db.index")
    output: 
        os.path.join(OUT_DIR, "search_results.tsv")
    params:
        results_basename = os.path.join(OUT_DIR, "results-db/results-db")
    threads:
        MMSEQS_THREADS
    shell:
        """
        mmseqs convertalis {input[0]} {input[0]} \
        {params.results_basename} {output} --threads {threads} \
        --format-output 'qset,query,tset,target,nident,alnlen,mismatch,gapopen,evalue,qlen'
        """


rule process_results:
    input: 
        os.path.join(OUT_DIR, "search_results.tsv")
    output: 
        os.path.join(OUT_DIR, "best_hits.csv"),
        os.path.join(OUT_DIR, "best_bidirectional_hits.csv")
    params: eval_threshold = config.get("eval_filter", 1.0E-15),
            coverage_threshold = config.get("coverage_filter", 0.7),
            identity_threshold = config.get("identity_filter", 0.3)
    script: "scripts/process_results.py"


rule compute_ani_af:
    input:
        os.path.join(OUT_DIR, "best_hits.csv"),
        os.path.join(OUT_DIR, "best_bidirectional_hits.csv"),
        os.path.join(OUT_DIR, "genome_lengths.csv")
    output:
        os.path.join(OUT_DIR, "ani.csv")
    script: "scripts/compute_ani_af.py"        