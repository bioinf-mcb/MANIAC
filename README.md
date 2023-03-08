# ANImm
Computation of average nucleotide identity with the use of MMseqs2.

**wGRR comment**
**To run wGRR put your protein sequences into one file per phage.**
**mmseqs --search-flag in cds-based.yml had to be set to one.**


## How to run
This pipeline is built with the use of [Snakemake](https://snakemake.github.io/) and can be run (for example) with the following command:

```
snakemake --cores 8 --snakefile <path to Snakefile in this repository> --configfile <path to configuration file>
```

## Configuration file
The configuration file is expected to be a yaml file, in which the following options can be specified. A [set of examples](./sample-configs) covering the two basic types of calculation:
* fragment-based.yml: fragment-based ANI calculation, based on Goris et al (PMID 17220447),
* cds-based.yml: best bidirectional hit calculation using predicted CDS as fragments.

Required:
*  `input_dir`: directory in which the genomes are present (one genome per file)
* `output_dir`: where the output should be written
* `bbh`: use best bidirectional hits to calculate ANI (only when CDS are provided instead of full genomes; True or False)

Filtering (optional):
* `eval_filter`: maximal allowed E-value (default: `1.0E-15`)
* `coverage_filter`: minimal query coverage (default: `0.7`)
* `identity_filter`: minimal query identity (default: `0.3`)

Other optional:
* `input_extension`: file extension of the genome files (default: `fna`)
* `fragment_size`: length of the genome fragments to be used in search (default: `1020`)
* `mmseqs_params`: any additional parameters to be passed to MMseqs2 search (default: `--search-type 3 -a --max-seqs 10000 --max-seq-len 100000 -s 7.5 --mask 0 -e 1e-15 -k 13 --zdrop 150`; calibrated with Pyani)
* `threads`: number of threads to be used for computationally intensive steps. Will be reduced if is greater than the `cores` parameter of Snakemake (default: `8`)
* `Low_memory_mode` : mode used to run ANImm in a memory stringent manner.Only loads table columns that are important for the analysis and drops all columns that are not used for ANI calculation(default:`false` for fragment_based and `true` for tr and cds)
For more sensitive search it is recommended to use higher sensitivity settings than default (such as `-s 7.5`) as well as the blastn scoring matrix (provided in this repository).

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
            <td>~2HRS to parse results table</td>
        </tr>
    </tbody>
  </table>