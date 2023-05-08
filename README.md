# ANImm
Snakemake workflow for computation of average nucleotide identity with the use of MMseqs2 optimized for phages. 


### Installation and execution (Linux)

#### Clone repository and install dependencies **(not tested)**.

```
git clone https://github.com/bioinf-mcb/ANImm
conda install -c conda-forge -c bioconda snakemake mamba biopython=1.79 pathlib=1.0.1 pandas
```

#### Test workflow

```
# fragment based
snakemake --use-conda --cores 8 --snakefile MANIAC --configfile test/configs/fragment-based.yml

# orf based
snakemake --use-conda --cores 8 --snakefile MANIAC --configfile test/configs/orf-based.yml

# cds based
snakemake --use-conda --cores 8 --snakefile MANIAC --configfile test/configs/cds-based.yml
```


### Configuration file
The configuration file is expected to be a yaml file, in which the following options can be specified. A [set of examples](./test/configs) covering the two basic types of calculation:

* fragment-based.yml: fragment-based ANI calculation, based on Goris et al (PMID 17220447)
* orf-based.yml: best bidirectional hit calculation open reading frames (ORFs).
* cds-based.yml: best bidirectional hit calculation using protein sequences (CDSs).


*Header format*
Each record header in input file has to be unique and follow a convenction (examples in test/data). 

* fragment-based.yml: any set of unique headers
* orf-based.yml: {PHAGEID}_ORF_{NUMBER}
* cds-based.yml: {PHAGEID}_PROTEIN_{NUMBER}


##### Configuration files details

Required:
* `INPUT_FILE`: one file with sequences (genomes/ORFs/CDSs) *[WARNING! Look at header formatting paragraph]*
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

##### Run time details
* Data type and runtime

  <table>
    <thead>
      <tr>
        <th>data</th>
        <th>type</th>
        <th>RAM</th>
        <th>run_time</th>
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