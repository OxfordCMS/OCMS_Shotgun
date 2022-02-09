# OCMS pipelines for shotgun metagenomic data analysis

## Pipelines

* pipeline_classifyreads.py - Classify reads from paired-end fastq files using kraken2
* pipeline_humann.py - functional profiling of fastq file using humann3

### Running pipeline_humann.py on BMRC
Set up environment on BMRC with

```
module load Python/3.8.2-GCCcore-9.3.0
source ~/devel/venv/Python-3.8.2-GCCcore-9.3.0/${MODULE_CPU_TYPE}/bin/activate
module load Bowtie2/2.4.1-GCC-9.3.0
module load DIAMOND/0.9.36-GCC-9.3.0
module load R/4.0.0-foss-2020a
module load R-bundle-Bioconductor/3.11-foss-2020a-R-4.0.0
module load Pandoc/2.13
module load X11/20200222-GCCcore-9.3.0
```

If need to concatenate paired-end reads, can run concatenate_fastq.py in OCMS_Sandbox/Py_utility
```
# will concatenate paired-end fastqs based on file name
## example: sample1.fastq.1.gz sample1.fastq.2.gz
python /path/to/OCMS_Sandbox/Py_utility/concatenate_fastq.py make full -p10 -v5
```

Run pipeline_humann.py in working directory containing fastq files
```
python /path/to/OCMS_Shotgun/pipelines/pipeline_humann.py make full -p10 -v5
```

Build report
```
python /path/to/OCMS_Shotgun/pipelines/pipeline_humann.py make build_report
```
