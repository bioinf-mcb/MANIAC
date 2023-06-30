# MANIAC

### What is MANIAC?
MANIAC stands for **M**Mseqs2-based **A**verage **N**ucleotide **I**dentity **A**ccurate **C**alculator. It is a bioinformatic pipeline, written using SnakeMake, for rapid and accurate calculation of average nucleotide identity (ANI) and Alignment Fraction (AF) between viral genomes. The goal of MANIAC is to provide a user-friendly and efficient tool for researchers in genomics, bioinformatics, and virology.

### Features
- High throughput: MANIAC can efficiently process large datasets (thousands) of viral genomes.
- Accurate: Uses MMseqs2 to ensure accurate calculation of average nucleotide identity (ANI) and alignment fraction (AF).
- Comprehensive: Provides analysis at both nucleotide and amino-acid level.
- User-friendly: Easy-to-use Snakemake workflow.
- Reproducible: Conda-based installation support ensures reproducibility.

### Prerequisities (to be verified)
- mmseqs2
- snakemake
- biopython=1.79
- pathlib=1.0.1
- pandas

### Installation (Mac OS X & Linux)

#### Clone repository and install dependencies **(not tested)**.
First clone the Github directory
```
git clone https://github.com/bioinf-mcb/MANIAC
```
Then create a conda environment with required dependencies, activate it and install the dependencies.
```
conda create --name maniac
conda activate maniac
conda install -c conda-forge -c bioconda snakemake mamba biopython=1.79 pathlib=1.0.1 pandas mmseqs2
```
NOTE: If you're using Apple M1/M2 computer, you may get a `PackagesNotFoundError`. If this should happen, run the following command:
```
conda config --add subdirs osx-64
```
and then run `conda install ...`. Hopefully you're good to go.

### Execution
MANIAC operates in three modes:

#### Fragment-based ANI
In this mode, MANIAC calculates ANI following the approach of Goris et al. [1]. [DETAILS]. To run the fragment-based mode, type 
```
snakemake --use-conda --cores 8 --snakefile MANIAC --configfile test/configs/fragment-based.yml
```

#### BBH using nucleotide-based CDS
In this mode, MANIAC calculates ANI using best-bidirectional hits based on the user-provided CDS (nucleotide level). 
```
snakemake --use-conda --cores 8 --snakefile MANIAC --configfile test/configs/orf-based.yml
```

#### BBH using amino acid-based CDS
In this mode, MANIAC calculated AAI (average amino-acid identity) using best-bidirectional hits based on the user-provided CDS (amino-acid level). 
```
snakemake --use-conda --cores 8 --snakefile MANIAC --configfile test/configs/cds-based.yml
```

### Configuration file
The configuration file is expected to be a yaml file, in which the various options can be specified. Each record header in input file has to be unique and follow a convenction. 

Examples of [config files](./test/configs) and [headers formatting](./test/data).

* fragment-based: fragment-based ANI calculation, based on Goris et al (PMID 17220447) (headers: any unique set)
* orf-based: best bidirectional hit calculation open reading frames (ORFs). (headers: PHAGEID_ORF_NUMBER)
* cds-based: best bidirectional hit calculation using protein sequences (CDSs). (headers: PHAGEID_PROTEIN_NUMBER)


### Details

#### Configuration files

Required:
* `INPUT_FILE`: one file with sequences (genomes/ORFs/CDSs)
* `OUTPUT_DIR`: where the output should be written
* `CDS_BASED`: use best bidirectional hits to calculate ANI (only when ORFs/CDSs are provided instead of full genomes) [True/False]

Filtering (optional):
* `EVALUE`: maximal allowed E-value (default: `1.0E-15`)
* `COVERAGE`: minimal query coverage (default: `0.7`)
* `IDENTITY`: minimal query identity (default: `0.3`)

Other optional:
* `FRAGMENT_SIZE`: length of the genome fragments to be used in search (default: `1020`)
* `MMSEQS_PARAMS`: any additional parameters to be passed to MMseqs2 search (default: `--search-type 3 -a --max-seqs 10000 --max-seq-len 100000 -s 7.5 --mask 0 -e 1e-15 -k 13 --zdrop 150`; calibrated with Pyani)
* `THREADS`: number of threads to be used for computationally intensive steps. Will be reduced if is greater than the `cores` parameter of Snakemake (default: `8`)
* `MEMORY_EFFICIENT` : mode used to run ANImm in a memory stringent manner. Only loads table columns that are important for the analysis and drops all columns that are not used for ANI calculation.

For more sensitive search it is recommended to use higher sensitivity settings than default (such as `-s 7.5`) as well as the blastn scoring matrix (provided in this repository).

#### Run time

  <table>
    <thead>
      <tr>
        <th>data</th>
        <th>type</th>
        <th>RAM</th>
        <th>time</th>
      </tr>
    </thead>
    <tbody>
        <tr>
            <td>"Inphared(~23k)"</td>
            <td>diverse</td>
            <td>~270GB</td>
            <td>~30HRS</td>
        </tr>
        <tr>
            <td>"Klebsiella prophages(~8K)"</td>
            <td>similar</td>
            <td>~270GB</td>
            <td>~7HRS to parse results table</td>
        </tr>
    </tbody>
  </table>

#### Results notes

* Fragment-based calculation has duplicate entries (a-b & b-a) (all modes do?)
* ANI for proteins (CDS) is AAI
* wGRR for ORFs is not *sensu stricto* wGRR

