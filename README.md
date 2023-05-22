<center><img src="Logo.png" alt="Taxonomic classification with XOR filters" width="500"/></center>

# Taxor: Fast and space-efficient taxonomic classification of long reads with hierarchical interleaved XOR filters

## Citation

## Table of contents

* [Description](#description)
* [Installation](#installation)
* [Commands](#commands)
* [Quickstart](#quickstart)

## Installation

The easiest way is to simply download [executable binaries](https://github.com/shenwei356/kmcp/releases) for Linux x86_64. 

## Commands

|Subcommand                                                                |Function                                                        |
|:-------------------------------------------------------------------------|:---------------------------------------------------------------|
|[**build**](#build)                                                       | Construct HIXF index from fasta reference files                |
|[**search**](#search)                                                     | Search sequences against a database index                      |
|[**profile**](#profile)                                                   | Generate the taxonomic profile from search results             |

```
genome_updater.sh \
    -d "refseq"\
    -g "archaea,bacteria,fungi,viral" \
    -c "all" \
    -l "complete genome,chromosome" \
    -f "genomic.fna.gz" \
    -o "refseq-abv" \
    -t 12 \
    -A "species:1" \
    -m -a -p

# cd to 2021-09-30_19-35-19

# taxdump
mkdir -p taxdump
tar -zxvf taxdump.tar.gz -C taxdump

cut -f 1,7,20 ../assembly_summary.txt \
| taxonkit lineage -i 2 -r -n -L --data-dir taxdump \
| taxonkit reformat -I 2 -P -t --data-dir taxdump \
| cut -f 1,2,3,4,6,7 > refseq_accessions_taxonomy.csv

```
