*In development*
# CELLO_pipeline
Nextflow based CELLO-seq pipeline optimised for SLURM clusters. 
If steps have same number, they are synchronous. 

## Step I: pre demultiplexing 
### How to run
1. Edit parameters.json
```
{
    "singularity": "/path/to/sarlacc.img", 
    "path": "/path/to/your/folder",

    "queue": "test/short/long",
    "cpus": "#",
    "time": "#day ##hours ##minutes ##seconds",
    "memory": "### GB",
    "minimap_reference_index": "/path/to/genome.fa.mmi",
    "experiment_name": "experiment_name",
    "input_files" : "Common ending of input files: *.fastq.gz"

}
```
2. Run nextflow
```
module load nextflow
nextflow -bg run step_1.nf -params-file parameters.json > step_1.log
```
- bg: background, enables run to continue even if you log out from cluster
- stdout is saved into step_1.log
- If you see barcode_\*.fastq and barcode_\*.fastq.rds , process is done.
- Check run_report_YYYY-MM-DD_hh-ss.html for recap
- Only remove work/ dir once you are happy with the outcome
```
rm -rf work/
```

### Description
Input: reads (*fastq.gz)
1. Concatenates them into one file
1. Removes reads >20kb
2. Yields % contamination in .txt
2. Runs dT and TSO adaptor qc (see .html)
2. Demultiplexes plate into barcodes
Output: barcode_\*.fastq and barcode_\*.fastq.rds

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
