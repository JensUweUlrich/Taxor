<center><img src="img/Logo.png" alt="Taxonomic classification with XOR filters" width="300"/></center>

# Taxor: Fast and space-efficient taxonomic classification of long reads with hierarchical interleaved XOR filters

## Citation

## Table of contents

* [Description](#description)
* [Installation](#installation)
* [Commands](#commands)
* [Usage](#usage)

## <a name="installation"></a>Installation

The easiest way is to simply download [executable binaries](https://github.com/shenwei356/kmcp/releases) for Linux x86_64. 

## <a name="commands"></a>Commands

|Subcommand                                                                |Function                                                        |
|:-------------------------------------------------------------------------|:---------------------------------------------------------------|
|[**build**](#build)                                                       | Construct HIXF index from fasta reference files                |
|[**search**](#search)                                                     | Search sequences against a database index                      |
|[**profile**](#profile)                                                   | Generate the taxonomic profile from search results             |


### <a name="build"></a>Taxor build

```
taxor-build - Creates and HIXF index of a given set of fasta files
==================================================================

DESCRIPTION
    Creates an HIXF index using either k-mers or syncmers

OPTIONS

  Basic options:
    -h, --help
          Prints the help page.
    -hh, --advanced-help
          Prints the help page including advanced options.
    --version
          Prints the version information.
    --copyright
          Prints the copyright/license information.
    --export-help (std::string)
          Export the help page information. Value must be one of [html, man].

  Main options:
    --input-file (std::string)
          tab-separated-value file containing taxonomy information and reference file names
    --input-sequence-dir (std::string)
          directory containing the fasta reference files Default: .
    --output-filename (std::string)
          A file name for the resulting index. Default: .
    --kmer-size (signed 32 bit integer)
          size of kmers used for index construction Default: 20. Value must be in range [1,30].
    --syncmer-size (signed 32 bit integer)
          size of syncmer used for index construction Default: 10. Value must be in range [1,26].
    --threads (signed 32 bit integer)
          The number of threads to use. Default: 1. Value must be in range [1,32].
    --use-syncmer
          enable using syncmers for smaller index size
```

<b> input-file</b><br>
This file contains all relevant information about the organisms in the database, which will be indexed. All values are tab-separated and the file should have following columns:
* Column 1: Assembly accession: the assembly accession.version reported in this field is 
   a unique identifier for the set of sequences in this particular version of 
   the genome assembly.
* Column 2: Taxonomy ID: the NCBI taxonomy identifier for the organism from which the 
   genome assembly was derived. The NCBI Taxonomy Database is a curated 
   classification and nomenclature for all of the organisms in the public 
   sequence databases. The taxonomy record can be retrieved from the NCBI 
   Taxonomy resource:
   https://www.ncbi.nlm.nih.gov/taxonomy/
 * Column 3: FTP path: the path to the directory on the NCBI genomes FTP site from which 
   data for this genome assembly can be downloaded
 * Column 4: Organism name
 * Column 5: Taxonomy string
 * Column 6: Taxonomy ID string

A two-line example of such a file is provided below. You can easily create such a file by following the preprocessing steps described in the [Usage](#usage) section.

```
GCF_000002495.2	318829	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/495/GCF_000002495.2_MG8	Pyricularia oryzae	k__Eukaryota;p__Ascomycota;c__Sordariomycetes;o__Magnaporthales;f__Pyriculariaceae;g__Pyricularia;s__Pyricularia oryzae	2759;4890;147550;639021;2528436;48558;318829
GCF_000002515.2	28985	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/515/GCF_000002515.2_ASM251v1	Kluyveromyces lactis	k__Eukaryota;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycetaceae;g__Kluyveromyces;s__Kluyveromyces lactis	2759;4890;4891;4892;4893;4910;28985
```

<b> input-sequence-dir</b><br>
Path to the directory containing fasta files (compressed) of organisms listed in the tab-separated file explained above. The file stem of the fasta files needs to match the last directory path string of the FTP path in column 3 of the input file (e.g. GCF_000002495.2_MG8)

<b> output-filename</b><br>
Path to the output file containing the hierarchical interleaved XOR filter index of the reference sequences and taxonomy information for the profiling step.

<b> kmer-size</b><br>
Size of k-length-substrings used for pseudo-mapping. When using syncmers for downsampling, the kmer-size has to be even-numbered because of using open canonical syncmers. The maximum supported k-mer size is 30.

<b> syncmer-size</b><br>
Size of the substrings used for selecting a k-mer for pseudo-mapping. The syncmer-size also has to be even-numbered because of the usage of open canonical syncmers. This number needs to be smaller than the k-mer size and the maximum supported size is 26.

<b> use-syncmer</b><br>
Switch that enables the usage of syncmers for downsampling of k-mers.

<b> threads</b><br>
Number of threads used for computing the hierarchical structure and building the HIXF index.

### <a name="search"></a>Taxor search



## <a name="usage"></a>Usage

![](img/Workflow.png)

First download the reference sequences and taxonomy dump of the sequences from the NCBI using [genome_updater](https://github.com/pirovc/genome_updater). 
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
```
Then, unpack the taxonomy dump and create a tab-separated-values file using the Linux command [cut](https://man7.org/linux/man-pages/man1/cut.1.html) and [taxonkit](https://github.com/shenwei356/taxonkit).

```
# cd to 2023-03-15_12-56-12

# taxdump
mkdir -p taxdump
tar -zxvf taxdump.tar.gz -C taxdump

cut -f 1,7,20 assembly_summary.txt \
| taxonkit lineage -i 2 -r -n -L --data-dir taxdump \
| taxonkit reformat -I 2 -P -t --data-dir taxdump \
| cut -f 1,2,3,4,6,7 > refseq_accessions_taxonomy.csv

```
Now we can build the hierarchical interleaved XOR filter (HIXF) index of the reference sequences and the NCBI taxonomy.
```
taxor build --input-file refseq_accessions_taxonomy.csv --input-sequence_dir refseq/2023-03-15_12-56-12/files \
--output-filename refseq-abfv-k22-s12.hixf --threads 6 --kmer-size 22 --syncmer-size 12 --use-syncmer
```
Then, we query the sample fastq file against the index allowing in this case a sequencing error rate of 15%. 
```
taxor search --index-file refseq-abfv-k22-s12.hixf --query-file SAMPLE.fq.gz \
--output-file SAMPLE.search.txt --error-rate 0.15 --threads 6 
```
Finally, the query result file is used as input for taxonomic profiling, which has three output files containing taxonomic abundances and sequence abundances in CAMI report format as well as a binning file with final read to reference assignments.
```
taxor profile --search-file SAMPLE.search.txt --cami-report-file SAMPLE.report \
--seq-abundance-file SAMPLE.abundance --binning-file SAMPLE.binning --sample-id SAMPLE
```
