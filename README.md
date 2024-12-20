# MANIAC
<p align="center"><img src="https://github.com/bioinf-mcb/MANIAC/blob/main/extras/maniac-logo.png" alt="MANIAC" width="500"></p>

## 1. What is MANIAC?
MANIAC stands for **M**Mseqs2-based **A**verage **N**ucleotide **I**dentity **A**ccurate **C**alculator. It is a bioinformatic pipeline, written using SnakeMake, for rapid and accurate calculation of average nucleotide identity (ANI) and Alignment Fraction (AF) between viral genomes. The goal of MANIAC is to provide a user-friendly and efficient tool for researchers in genomics, bioinformatics, and virology. MANIAC has been developed and optimised for bacteriophages but in principle can be used on any microbial genomes.

Publication: https://doi.org/10.1101/2024.04.23.590796<br>
Contact: jade.leconte@uj.edu.pl

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


### Best-bidirectional hits or CDS mode
<p align="center"><img src="https://github.com/bioinf-mcb/MANIAC/blob/main/extras/ani-bbh.png" alt="BBH" width="600"></p>

In addition to the standard, fragment-based ANI calculation, MANIAC carries out the calculation using best-bidirectional hits approach should the user provide coding sequences (CDSs) for input genomes, either in nucleotide or amino-acid. The calculation is then carried out analogously as in the fragment mode with the following differences:

1. CDS are being used instead of fragments
2. To calculate ANI and AF, in both query and subject only CDSs which are each others best hits are considered.


## 4. Install

### Linux

Create and activate a conda environment.
```
conda create -n maniac -c conda-forge mamba python=3.9
conda activate maniac
mamba install -c conda-forge -c bioconda snakemake pandas biopython=1.79 mmseqs2 r-base r-essentials r-arrow datamash pyopenssl=24.2 parallel=20240922
```

Clone MANIAC repository. 
```
git clone https://github.com/bioinf-mcb/MANIAC
```

Test using example input data and configuration files in the `test` folder.
```
snakemake --cores 8 --quiet --snakefile MANIAC --configfile test/configs/easy-fragment-based.yml
snakemake --cores 8 --quiet --snakefile MANIAC --configfile test/configs/easy-cds-aa.yml
snakemake --cores 8 --quiet --snakefile MANIAC --configfile test/configs/easy-cds-nt.yml
```

### macOS

Install dependencies using [homebrew](https://brew.sh/)
```
brew update
brew install bash
brew install coreutils
brew install gnu-sed
brew install gawk
brew install parallel
brew install datamash
brew install mmseqs2

export PATH="$(brew --prefix coreutils)/libexec/gnubin:$PATH"
export PATH="$(brew --prefix gnu-sed)/libexec/gnubin:$PATH"
export PATH="$(brew --prefix gawk)/libexec/gnubin:$PATH"
export PATH="/usr/local/opt/mmseqs2/bin:$PATH"
```

Create and activate a conda environment.
```
conda create -n maniac -c conda-forge mamba python=3.9
conda activate maniac
mamba install -c conda-forge -c bioconda snakemake pandas biopython=1.79 r-base r-essentials r-arrow datamash pyopenssl=24.2
```

Clone MANIAC repository. 
```
git clone https://github.com/bioinf-mcb/MANIAC
```

Test using example input data and configuration files in the `test` folder.
```
snakemake --cores 8 --quiet --snakefile MANIAC --configfile test/configs/easy-fragment-based.yml
snakemake --cores 8 --quiet --snakefile MANIAC --configfile test/configs/easy-cds-aa.yml
snakemake --cores 8 --quiet --snakefile MANIAC --configfile test/configs/easy-cds-nt.yml
```

### Windows

To install MANIAC on Windows, you first need to install Windows Subsystem for Linux (WSL) and set it up. Once WSL is installed, follow the instructions for installing MANIAC on Linux.

1. Click the Start menu, type "PowerShell," right-click on Windows PowerShell, and select Run as administrator.
2. In the PowerShell window, enter the following command `wsl --install` to install WSL.
3. Restart Your Computer, choose Linux to launch and follow the on-screen instructions.
4. Once your Linux environment is ready, follow the [Linux](#linux) Debian-Based installation steps to install MANIAC.


### Dependencies details:
MANIAC was successfully tested on Linux, macOS and Windows with the following dependencies versions. Make sure the bash version used is 5.0 or above by running `bash --version` and that all dependencies were successfully installed with homebrew and/or conda.
- python=3.9
- bash=5.2.21
- r-base=4.4.1
- r-essentials=4.4
- r-arrows=17.0.0
- snakemake=7.32.4
- pandas=2.2
- biopython=1.79
- mmseqs2=15.6
- datamash=1.8
- pyopenssl=24.2
- parallel=20240922

## 5. Running MANIAC
This section will guide you on how to prepare your input files, create a yaml configuration file, and run the MANIAC software. We'll also cover the types of output files you can expect from MANIAC.

### Input files
MANIAC requires one of two types of input files: 

1. Full genome files (for the fragment calculation), 
2. Nucleotide or amino-acid coding-sequences (for the BBH calculation).

Each file should be in FASTA format. The header convention for CDS input should be the genome name, followed by a `_CDS` string, followed by its unique suffix. For example, if genome named **XYZ_phageVp123** has three coding sequences, the input file headers could be  

`>XYZ_phageVp123_CDS1`, `>XYZ_phageVp123_CDS2` and `>XYZ_phageVp123_CDS5`

Examples of input files are located in `test/data`.

### Configuration file
MANIAC uses a yaml configuration file to set the workflow parameters. Here's an example of what a simple configuration file might look like:

```
INPUT_FILE: "test/data/fragment-based.fasta"
OUTPUT_DIR: "test/output/FRAGMENT-BASED"
MODE: FRAGMENTS_NT
FAST: False
MEMORY_GB: 16
```
Here are details of various parameters.

#### Parameters: required
* `INPUT_FILE`: full genome or CDS file
* `OUTPUT_DIR`: directory where the output should be written
* `MODE`: FRAGMENTS_NT requires full genomes as an input, while CDS_NT and CDS_AA use BBH to calculate ANI and require the input to be CDS (nucleotide or protein respectively) [FRAGMENTS_NT | CDS_NT | CDS_AA]
* `FAST`: Enable Fast mode. Fast mode will automatically overwrite some parameters to prioritize speed over accuracy (KMER: 15, FRAGMENT_SIZE: 1020) [True/False]
* `MEMORY_GB`: Declare the memory available for MANIAC in GB. Note: This will not actually limit the memory and is only used to optimize post-processing run speed. (default: `16`)

#### Parameters: specific to fragment mode (optional)
* `COVERAGE`: minimal query coverage used for filtering (default: `0.7`)
* `IDENTITY`: minimal query identity used for filtering (default: `0.3`)
* `FRAGMENT_SIZE`: length of the genome fragments to be used in search (default: `500`)

#### Parameters: specific to BBH mode (optional)
* `HOMOLOGS:` BBH & homologous CDS definition
  * `IDENTITY`: (default: `0.3`)
  * `COVERAGE`: (default: `0.7`)
* `CONSERVED`: conservative CDS definition
  * `IDENTITY`: (default: `0.8`)
  * `COVERAGE`: (default: `0.5`)

#### Parameters: others (optional)
* `DELETE_INTERMEDIATE_FILES`: [True/False] (default: `True`)
* `MEMORY_EFFICIENT`: mode used to run in a memory stringent manner. Only loads table columns that are important for the analysis and drops all columns that are not used for ANI calculation [True/False] (default: `True`)
* `MMSEQS_PARAMS`: any additional parameters to be passed to MMseqs2 search, default values calibrated with Pyani
  * `EVALUE`: (default: `1e-15`)
  * `SENSITIVITY`: (default: `7.5`)
  * `ZDROP`: (default: `40`)
  * `MAX_SEQS`: (default: `10000`)
  * `MAX_SEQ_LEN`: (default: `65000`)
  * `KMER`: (default: `11`)
  * `SEED_SUB_MATRIX`: (default: `scoring/blastn-scoring.out`)
  * `SUB_MATRIX`: (default: `scoring/blastn-scoring.out`)

#### Parameters: recommendations
For full genome and nucleotide CDS mode, the alignment scoring matrix should be provided. Matrices for the blastn and unit-scoring modes are provided in the repository. Please note that the sensitivity parameter will not matter for nucleotide-based calculations, only k-mer size will. If FAST is enabled, k-mer size will be forced to 15.

For amino-acid calculations, no scoring matrix has to be provided but a more sensitive search is recommended (such as `-s 7.5` or higher). Please refer to the original mmseqs publication [2].

Examples of input files for different calculation modes are located in `test/configs`. A minumum working example is provided, as well as different examples with more complete sets of parameters for advanced users. **We strongly recommend against changing the mmseqs input parameters** as they have been optimised for different calculation modes.


### Running MANIAC
After your input files are ready and your configuration file is set, you can run MANIAC as follows:
```
snakemake --cores 8 --quiet --snakefile MANIAC --configfile your-path-to-configuration-file.yml
```
where `your-path-to-configuration-file.yml` is the full path to your configuration file. The type of the configuration file will determine whether MANIAC runs in the fragment mode or the BBH mode. `cores` should be adapted to the machine you are using to run MANIAC.

If MANIAC ran into an error because of an issue (such as wrong installation) that you were able to correct, please delete the output directory before trying to run the pipeline again.

### Maniac Output Description

Maniac generates output files in the user-defined output directory. The `genome-alignment.csv` file contains the ANI results along with associated metrics. The file is a table with fields detailed below:

| Metrics        | Description                                                                                                                                                                         |
|----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **ANI/AAI**    | Average nucleotide or aminoacid identity between the query and reference sequences                                                                                                  |
| **len_1**      | The length of the query sequence                                                                                                                                                    |
| **len_2**      | The length of the reference sequence                                                                                                                                                |
| **ani_alnlen** | The total length of aligned nucleotides between the query and reference sequences                                                                                                    |
| **af_1**       | Alignment fraction of the query sequence calculated by dividing the aligned length by the total length of the query sequence                                                                                      |
| **af_2**       | Alignment fraction of the reference sequence calculated by dividing the aligned length by the total length of the reference sequence                                                                                  |
| **af_min**     | The minimum alignment fraction between the query and reference sequence calculated by dividing the aligned nucleotide length by the shorter sequence between the query and reference sequence                 |
| **af_max**     | The maximum alignment fraction between the query and reference sequence calculated by dividing the aligned nucleotide length by the longer sequence between the query and reference sequence                   |
| **af_mean**    | Mean alignment fraction between the query and reference sequences. It is calculated by averaging the alignment fraction of both query and reference sequences weighted by their length. Users can also calculate `af_mean` by considering the alignment fraction between pairs since the results of MANIAC are asymmetrical i.e (af_1 + af_2)/2                                                                                                                                 |
| **af_jaccard** | The jaccard index of the alignment fraction calculated as the ratio of the aligned length to the total length of the union of the query and reference sequences                                                              |
| **wgANI/wgAAI**        | wgANI/AAI is the whole genome ANI/AAI. It is calculated by multiplying ANI or AAI by the mean AF |
| **seq1_n_prots** | Number of proteins or CDS in the query sequence |
| **seq2_n_prots** | Number of proteins or CDS in the reference sequence |
| **seq1_n_prots_hom** | Number of homologous proteins or CDS in the query sequence |
| **seq2_n_prots_hom** | Number of homologous proteins or CDS in the reference sequence |
| **seq1_n_prots_cons** | Number of conservative proteins or CDS in the query sequence |
| **seq2_n_prots_cons** | Number of conservative proteins or CDS in the reference sequence |
| **cds_alignments_counts** | Number of alignments proteins or CDS between query and reference sequences |
| **cds_alignments_ani_sum** | Sum of nucleotide or aminoacid identities of aligned proteins or CDS between query and reference sequences |
| **min_prots**    | The minimum number of proteins or CDS between the query and reference sequences |
| **wgrr**         | wGRR is the weighted gene repertoire relatedness. It is calculated as the ratio of bi-directional best hits between the query and reference genomes weighted by the sequence identity of homologs (CDS or protein homologs for the CDS or protein mode respectively) |


## 6. References
1. Goris, J. et al. DNA-DNA hybridization values and their relationship to whole-genome sequence similarities. Int. J. Syst. Evol. Microbiol. 57, 81–91 (2007).
2. Steinegger, M. & Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat. Biotechnol. 35, 1026–1028 (2017).
3. Johannes Köster, Sven Rahmann, Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, Volume 28, Issue 19, Pages 2520–2522 (2012).
4. Richardson N, Cook I, Crane N, Dunnington D, François R, Keane J, Moldovan-Grünfeld D, Ooms J, Wujciak-Jens J, Apache Arrow. arrow: Integration to 'Apache' 'Arrow'. R package version 17.0.0, https://github.com/apache/arrow (2024).
5. Tange, O. GNU Parallel 20240922 ('Gold Apollo AR924'). Zenodo. doi: 10.5281/zenodo.13826092 (2024) 
