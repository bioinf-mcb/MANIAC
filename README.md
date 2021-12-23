# ANImm
Computation of average nucleotide identity with the use of MMseqs2.

## How to run
This pipeline is built with the use of [Snakemake](https://snakemake.github.io/) and can be run (for example) with the following command:

```
snakemake --cores 8 --snakefile <path to Snakefile in this repository> --configfile <path to configuration file>
```

## Configuration file
The configuration file is expected to be a yaml file, in which the following options can be specified. A [set of examples](./sample-configs) covering three basic scenarios (fragment-based ANI, CDS-based ANI and AAI) is provided.

Required:
*  `input_dir`: directory in which the genomes are present (one genome per file)
* `output_dir`: where the output should be written
* `bbh`: use best bidirectional hits to calculate ANI (only when CDS are provided instead of full genomes)

Filtering (optional):
* `eval_filter`: maximal allowed E-value (default: `1.0E-15`)
* `coverage_filter`: minimal query coverage (default: `0.7`)
* `identity_filter`: minimal query identity (default: `0.3`)

Other optional:
* `input_extension`: file extension of the genome files (default: `fna`)
* `fragment_size`: length of the genome fragments to be used in search (default: `1020`)
* `mmseqs_params`: any additional parameters to be passed to MMseqs2 search (default: `--search-type 3 -a --max-seqs 10000 --max-seq-len 100000 -s 7.5 --mask 0`)
* `threads`: number of threads to be used for computationally intensive steps. Will be reduced if is greater than the `cores` parameter of Snakemake (default: `8`)

For more sensitive search it is recommended to use higher sensitivity settings than default (such as `-s 7.5`) as well as unit match/mismatch scoring matrix (provided in this repository).