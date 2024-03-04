# MANIAC
<p align="center"><img src="https://github.com/bioinf-mcb/MANIAC/blob/main/extras/maniac-logo.png" alt="MANIAC" width="500"></p>

## 1. What is MANIAC?
MANIAC stands for **M**Mseqs2-based **A**verage **N**ucleotide **I**dentity **A**ccurate **C**alculator. It is a bioinformatic pipeline, written using SnakeMake, for rapid and accurate calculation of average nucleotide identity (ANI) and Alignment Fraction (AF) between viral genomes. The goal of MANIAC is to provide a user-friendly and efficient tool for researchers in genomics, bioinformatics, and virology. MANIAC has been developed and optimised for bacteriophages but in principle can be used on any microbial genomes.

## 2. Features
- High throughput: MANIAC can efficiently process large datasets (thousands) of viral genomes.
- Accurate: Uses MMseqs2 to ensure accurate calculation of average nucleotide identity (ANI) and alignment fraction (AF).
- Comprehensive: Provides analysis at both nucleotide and amino-acid level.
- User-friendly: Easy-to-use Snakemake workflow.
- Reproducible: Conda-based installation support ensures reproducibility.

## 3. ANI calculation
### Fragment mode
<p align="center"><img src="https://github.com/bioinf-mcb/MANIAC/blob/main/extras/ani-fragment.png" alt="fragment" width="600"></p>

The standard and quickest way of ANI calculation is based on the approach proposed by Goris et al. for bacterial genomes [1]. Specifically, each query is chopped into short fragments of pre-defined length (by default 1020 nt). Then, each fragment is aligned with the subject and the best hit is found – but only if the query coverage is at least 70% and the sequence identity is 30% across the entire query length. ANI is then taken as the mean percentage identity of all aligned fragments and query AF is calculated as the length of the aligned query genome (i.e., the summed length of all aligned fragments) to the full query length.


### Best-bidirectional hits mode
<p align="center"><img src="https://github.com/bioinf-mcb/MANIAC/blob/main/extras/ani-bbh.png" alt="BBH" width="600"></p>

In addition to the standard, fragment-based ANI calculation, MANIAC carries out the calculation using best-bidirectional hits approach should the user provide coding sequences (CDSs) for input genomes, either in nucleotide or amino-acid. The calculation is then carried out analogously as in the fragment mode with the following differences:

1. CDS are being used instead of fragments
2. To calculate ANI and AF, in both query and subject only CDSs which are each others best hits are considered.


## 4. Installation

### Clone repository
First clone the GitHub directory
```
git clone https://github.com/bioinf-mcb/MANIAC
```

### Install dependencies (conda)
#### Linux & MacOS

```
conda create -n maniac -c conda-forge mamba python=3.9
conda activate maniac
mamba install -c conda-forge -c bioconda snakemake pandas biopython=1.79 mmseqs2
```

#### Test
```
cd MANIAC
snakemake --cores 8 --snakefile MANIAC --configfile test/configs/fragment-based.yml
snakemake --cores 8 --snakefile MANIAC --configfile test/configs/cds-aa.yml
snakemake --cores 8 --snakefile MANIAC --configfile test/configs/cds-nt.yml
```


#### MANIAC dependencies:

- python=3.9
- snakemake=8.5
- pandas=2.2
- biopython=1.79
- python=3.11
- mmseqs2=15.6



## 5. Running MANIAC
This section will guide you on how to prepare your input files, create a yaml configuration file, and run the MANIAC software. We'll also cover the types of output files you can expect from MANIAC.

### Input files
MANIAC requires one of two types of input files: 

1. Full genome files (for the fragment calculation), 
2. Nucleotide or amino-acid coding-sequences (for the BBH calculation).

Each file should be in FASTA format. The header convention for CDS input should be the genome name, followed by a `_CDS` sting, followed by its unique suffix. For example, if genome named **XYZ_phageVp123** has three coding sequences, the input file headers could be  

`>XYZ_phageVp123_CDS1`, `>XYZ_phageVp123_CDS2` and `>XYZ_phageVp123_CDS5`

Examples of input files are located in `test/data`.

### Configuration file
MANIAC uses a yaml configuration file to set the workflow parameters. Here's an example of what this file might look like:

```
INPUT_FILE: "test/data/fragment-based.fasta"
OUTPUT_DIR: "test/output/FRAGMENT-BASED"

THREADS: 16
CDS_BASED: False
MEMORY_EFFICIENT: True
DELETE_CDS_ALIGNMENT: False
DELETE_INTERMEDIATE_FILES: True

EVALUE: 10
IDENTITY: 0.3
COVERAGE: 0.7

MMSEQS_PARAMS:
  '--search-type 3 -a --max-seqs 10000 --max-seq-len 100000 -s 7.5 --mask 0 -e 1e-15 -k 11 --zdrop 150 
  --seed-sub-mat "scoring/blastn-scoring.out" 
  --sub-mat "scoring/blastn-scoring.out"'
```

Here are details of various parameters.

#### Parameters: required
* `INPUT_FILE`: full genome or CDS files
* `OUTPUT_DIR`: where the output should be written
* `CDS_BASED`: use BBH to calculate ANI (only when ORFs/CDSs are provided instead of full genomes) [True/False]

#### Parameters: filtering (optional)
* `EVALUE`: maximal allowed E-value (default: `1.0E-15`)
* `COVERAGE`: minimal query coverage (default: `0.7`)
* `IDENTITY`: minimal query identity (default: `0.3`)

#### Parameters: other optional
* `FRAGMENT_SIZE`: length of the genome fragments to be used in search (default: `1020`)
* `MMSEQS_PARAMS`: any additional parameters to be passed to MMseqs2 search (default: `--search-type 3 -a --max-seqs 10000 --max-seq-len 100000 -s 7.5 --mask 0 -e 1e-15 -k 13 --zdrop 150`; calibrated with Pyani)
* `THREADS`: number of threads to be used for computationally intensive steps. Will be reduced if is greater than the `cores` parameter of Snakemake (default: `8`)
* `MEMORY_EFFICIENT` : mode used to run ANImm in a memory stringent manner. Only loads table columns that are important for the analysis and drops all columns that are not used for ANI calculation.

#### Parameters: recommendations

For full genome and nucleotide CDS mode, the alignment scoring matrix should be provided. Matrices for the blastn and unit-scoring modes are provided in the repository. Please note that the sensitivity parameter will not matter for nucleotide-based calculations, only k-mer size will.

For amino-acid calculations, no scoring matrix has to be provided but a more sensitive search is recommended (such as `-s 7.5` or higher). Please refer to the original mmseqs publication [2].

Examples of input files for different calculation modes are located in `test/configs`. **We strongly recommend against changing the mmseqs input parameters** as they have been optimised for different calculation modes.


### Running MANIAC
After your input files are ready and your configuration file is set, you can run MANIAC as follows:
```
snakemake --cores 8 --snakefile MANIAC --configfile your-path-to-configuration-file.yml
```
where `your-path-to-configuration-file.yml` is the full path to your configuration file. The type of the configuration file will determine whether MANIAC runs in the fragment mode or the BBH mode.

### Output files
(TO BE WRITTEN)


## 6. References
1. Goris, J. et al. DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. Int. J. Syst. Evol. Microbiol. 57, 81–91 (2007).
2. Steinegger, M. & Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat. Biotechnol. 35, 1026–1028 (2017).
  


## NOTES

* Fragment-based calculation has duplicate entries (a-b & b-a) (all modes do?)
* ANI for proteins (CDS) is AAI
* wGRR for ORFs is not *sensu stricto* wGRR

