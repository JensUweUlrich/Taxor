# PoreClassify

```
genome_updater.sh \
    -d "refseq"\
    -g "viral" \
    -c "all" \
    -l "all" \
    -f "genomic.fna.gz" \
    -o "refseq-viral" \
    -t 12 \
    -m -a -p

# cd to 2021-09-30_19-35-19

# taxdump
mkdir -p taxdump
tar -zxvf taxdump.tar.gz -C taxdump

cut -f 1,6 ../assembly_summary.txt \
| taxonkit lineage -i 2 -r -n -L --data-dir taxdump \
| taxonkit reformat -I 2 -P --data-dir taxdump \
| cut -f 1,3,4,5 > refseq_accessions_taxonomy.csv

mkdir -p files.renamed

cd files.renamed
find ../files -name "*.fna.gz" \
    | rush 'ln -s {}'
    
cd ..    
# rename
brename -R -p '^(\w{3}_\d{9}\.\d+).+' -r '${1}_genomic.fna.gz' files.renamed
```
