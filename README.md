*In development*
# CELLO_pipeline
Nextflow based CELLO-seq pipeline optimised for SLURM clusters. 
If steps have same number, they are synchronous. 

## Step I: pre demultiplexing 
### Description
Input: reads (*fastq.gz)
1. Concatenates them into one file
1. Removes reads >20kb
2. Yields % contamination in .txt
2. Runs dT and TSO adaptor qc (see .html)
2. Demultiplexes plate into barcodes
Output: barcode_\*.fastq and barcode_\*.fastq.rds

### How to run
```
echo hi
```

## Step II: post demultiplexing
Input: barcode_\*.fastq, barcode_\*.fastq.rds, tso.rds, dT.rds. 
- Analysis is repeated per barcode
1. dT adaptor filer
1. TSO adaptor filter
1. Align to reference genome
2. Grouping of reads by UMI
3. Error-correction of reads by UMI
Output: error-corrected and demultiplexed fastq reads. 


### Da fare 
1. Strategia d'errore - limitata
4. Porechop
5. Rallentare per non overwelmare il cluster
6. Aggiungere FLARE?
7. Nascondi piu cose come i barcodes?
8. Sputa fuori i html utili
9. crea folder per output?
10. pachettizza
