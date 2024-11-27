# OCMS pipelines for shotgun metagenomic data analysis
This repository contains a series of pipelines used for processing shotgun metagenomic data. Pipelines are written within the [CGAT framework](https://github.com/cgat-developers/cgat-core). OMCS_Shotgun has a command line interface, and can be installed and executed as a stand alone command line tool. OCMS_Shotgun is primarily written for usage within the OCMS on our HPC system, however, it can be used on other HPC systems, or used locally.

## Install
Clone the OCMS_Shotgun repository and install using pip, ideally within a python virtual environment.

```
# Download the repo
git clone https://github.com/OxfordCMS/OCMS_Shotgun.git

# Activate python virtual environment (if applicable) and install OCMS_Shotgun
cd OCMS_Shotgun
pip install .
```

## Pipeline Environments
Each pipeline has it's own set of dependencies. It is recommended that you only load the tools necessary for the pipeline being used. If you are working within the BMRC HPC, you can load the pipeline modulefile. Please see the [OCMS_Modulefiles](https://github.com/OxfordCMS/OCMS_Modulefiles) repository for details and to download the modulefiles and read the OCMS modulefiles SOP for more details. If you are not working within the BMRC, please ensure that you have the appropriate software installed and in your PATH.

## Quick Start
All pipelines are written to be used within a HPC system, but can be run using the `--local` flag to run locally. 

Set up the pipeline configuration file within your working directory.

```
ocms_shotgun preprocess config
```

You can see the pipeline tasks with `show full`. 

```
ocms_shotgun preprocess show full
```

Run pipeline individual pipeline tasks with `make` followed by the pipeline task or run all pipeline tasks with `make full`

```
ocms_shotgun preprocess make full -p 20 -v 5
```

# Pipelines

<details>
  <summary>Pipeline Databases</summary>

Each pipeline requires certain databases and indexes to run. For consistency between members of the group and to ensure compatibility with the different tool chains on the BMRC we have developed pipeline_databases. The pipeline is designed to either download pre-built index files (e.g. kraken2) or to download flat files and index them with a particular tool.

### Dependencies

The software dependencies are the same as for the pipelines below. If you are working on BMRC then please also see [OCMS_Modulefiles](https://github.com/OxfordCMS/OCMS_Modulefiles) for an explanation of how to set up your environment.

### Configuration

To create a pipeline.yml file for configuration type

```
ocms_shotgun databases config
```

Edit this file with versions of software for database creation.

### Run pipeline_databases

There are no inout files specified as the task dependencies are created on the fly. If you only want to create databases for a single pipeline then you can use, for example,

```
ocms_shotgun databases make buildPreprocessDatabases
```

which will build the databases for the preprocess pipeline. If you want to produce databases for every pipeline in OCMS_Shotgun repository then you can use

```
ocms_shotgun databases make full
```

Note that some database files that are downloaded from various repositories are very large and so may take some hours to download.

### Output

The structure of the output for the "full" pipeline is as follows:
```
microbiome/
|--genomes/
|  |-human/
|    |-build/
|      |-build.fa.gz
|   |-mouse/
|     |-build
|       |-build.fa.gz
|
|--SRPRISM/
|  |-human/
|    |-build/
|      |-GCC-version/
|        |-srprism-version/
|          |-build.srprism*
|  |-mouse 
|    |-build/
|      |-GCC-version/
|        |-srprism-version/
|          |-build.srprism*
|
|--bmtool/
|  |-human/
|    |-build/
|      |-GCC-version/
|        |-bmtool-version/
|          |-build.bitmask
|  |-mouse 
|    |-build/
|      |-GCC-version/
|        |-bmtool-version/
|          |-build.bitmask
|
|--kraken2/
|   |-kraken2-version/
|     |-database_version/
|       |-database_files
|
|--sortmerna/
|  |-GCC-version/
|    |-sortmerna-version/
|      |-rrna/
|        |-smr_*.fasta
|      |-index/
|        |-smr_*.db
```
</details>

<details>
  <summary>Pipeline Preprocess</summary>

## Pipeline Preprocess
This pipeline pre-processes shotgun metagenome or metatranscriptome data. It performes the following:

* summarise raw input read counts
* remove duplicate sequences with Cdhit
* removeAdapters with Trimmomatic
* remove rRNA with SortMeRNA
* remove host reads with SortMeRNA
* mask low complexity reads with bmtagger
* summrise preprocessed read counts

### Dependencies

Software requirements:

| Software      |
|---------------|
| CDHIT         |
| CDHITauxtools |
| SortMeRNA     |
| bmtagger      |
| BBMap         |
| SAMtools      |
| SRPRISM       |
| Trimmomatic   |

If using OCMS_Modulefiles you can simply load the modules:

```
module load pipelines/preprocess
```

### Configuration
Initiate the configuration file.

```
ocms_shotgun preprocess config
```

### Input files
Pipeline preprocess takes in single or paired end reads. Input files should use the notation `fastq.1.gz`, `fastq.2.gz`. Input files should be located in the working directory, alternatively, an input directory called `input.dir` can be specified in the yml with:

```
# pipeline.yml
location_fastq: 1
```

### Pipeline tasks

```
Task = "mkdir('read_count_summary.dir')   before pipeline_preprocess.countInputReads "
Task = 'pipeline_preprocess.countInputReads'
Task = "mkdir('reads_deduped.dir')   before pipeline_preprocess.removeDuplicates "
Task = 'pipeline_preprocess.removeDuplicates'
Task = "mkdir('reads_adaptersRemoved.dir')   before pipeline_preprocess.removeAdapters "
Task = 'pipeline_preprocess.removeAdapters'
Task = "mkdir('reads_rrnaRemoved.dir')   before pipeline_preprocess.removeRibosomalRNA "
Task = 'pipeline_preprocess.removeRibosomalRNA'
Task = "mkdir('reads_hostRemoved.dir')   before pipeline_preprocess.removeHost "
Task = 'pipeline_preprocess.removeHost'                                                 
Task = "mkdir('reads_dusted.dir')   before pipeline_preprocess.maskLowComplexity "
Task = 'pipeline_preprocess.maskLowComplexity'
Task = 'pipeline_preprocess.countOutputReads'
Task = 'pipeline_preprocess.collateReadCounts'
Task = 'pipeline_preprocess.summarizeReadCounts'
Task = 'pipeline_preprocess.full'         
```

### Run pipeline_preprocess
The pipeline must have input fastq files with the notation `.fastq.1.gz` and `pipeline.yml` in working directory. Set the number of jobs `-p` equal to the number of samples.

```
ocms_shotgun preprocess make full -p 20 -v 5
```

### Output
```
# Summary of reads remaining after each task
processing_summary.tsv

# output reads after serial filtering steps
reads_dusted.dir/

```

</details>

<details>
  <summary>Pipeline Kraken2</summary>

## Pipeline Kraken2
Uses Kraken2 to classify paired-end reads
Uses Bracken to estimate abundances at every taxonomic level
Uses Taxonkit to generate a taxonomy file listing taxonomic lineage in mpa style

### Dependencies
Taxonkit requires NCBI taxonomy files, which can be downloaded from the [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz). Path to directory of taxonomy files is specified in the `taxdump` parameter in the yml. 

Software requirements:

| Software	|
|---------------|
| Kraken2	|
| Bracken       |
| taxonkit	|


If using OCMS_Modulefiles you can simply load the modules:

```
module load pipelines/kraken2

```

### Configuration
Initiate the configuration file.

```
ocms_shotgun kraken2 config
```

### Input files
Pipeline kraken2 takes in single or paired end reads. Input files should use the notation `fastq.1.gz`, `fastq.2.gz`. Input files should be located in the working directory.

### Pipeline tasks

```
Task = "mkdir('taxonomy.dir')   before pipeline_kraken2.translateTaxonomy "
Task = "mkdir('bracken.dir')   before pipeline_kraken2.runBracken "
Task = 'pipeline_kraken2.runBracken'
Task = 'pipeline_kraken2.checkBrackenLevels'
Task = 'pipeline_kraken2.mergeBracken'
Task = 'pipeline_kraken2.translateTaxonomy'
Task = 'pipeline_kraken2.full'
```

### Run pipeline_kraken2
The pipeline must have input fastq files with the notation `.fastq.1.gz` and `pipeline.yml` in working directory. Set the number of jobs `-p` to 7 times the number of samples (so Bracken can be run on all taxonomic levels in parallel), however please be mindful of the number of jobs.

```
ocms_shotgun kraken2 make full -p 140 -v 5
```

### Output
```
# classified reads
kraken.dir/

# estimated abundances
bracken.dir/

# showing taxonomy as mpa-styled lineages
taxonomy.dir/

# counts with taxonomy information added as feature names
counts.dir/
```

</details>


<details>
  <summary>Pipeline Concatfastq</summary>

## Pipeline Concatfastq
This pipelines concatenates paired-end reads into one file. This is helpful when running Humann3.

### Dependencies
No dependencies

### Configuration
No configuration file needed

### Input files
Paired end reads should end in the notation `fastq.1.gz` and `fastq.2.gz`. Input files located in working directory.

### Run pipeline_concatfastq
Set number of jobs `-p` to the number of samples

```
ocms_shotgun concatfastq make full -p 20 -v 5
```

### Output
Concatenated fastq files located in `concat.dir/`

</details>

<details>
  <summary>Pipeline Humann3</summary>

## Pipeline Humann3
This pipeline performs functional profiling of fastq files using Humann3.

### Dependencies
This pipeline was written for Humann3 v3.8 and Metaphlan 3.1. If you're not working within BMRC, Humann3 and Metaphlan3 need to be installed according to their developers' instructions. 

Software requirements:

| Software	|
|---------------|
| Bowtie2       |
| MetaPhlAn     |
| humann        |
| DIAMOND       |
| R             |
| Pandoc        |

If using OCMS_Modulefiles you can simply load the modules:

```
module load pipelines/humann3
```

### Configuration
Initiate configuration file

```
ocms_shotgun humann3 config
```

### Input files
Humann3 takes in single end reads. If you have paired-end reads, paired-ends need to be concatenated into one file. Concatenating paired-end fastqs can be done with ` pipeline_concatfastq`. Input files should end in the notation `fastq.gz`, located in the working directory.

### Pipeline tasks

```
Task = "mkdir('humann3.dir')   before pipeline_humann3.runHumann3 "
Task = 'pipeline_humann3.runHumann3'
Task = 'pipeline_humann3.mergePathCoverage'
Task = 'pipeline_humann3.mergePathAbundance'
Task = 'pipeline_humann3.mergeGeneFamilies'
Task = 'pipeline_humann3.mergeMetaphlan'
Task = 'pipeline_humann3.splitMetaphlan'
```

### Run pipeline_humann3
Set number of jobs `-p` to number of samples.

```
ocms_shotgun humann3 make full -p 20 -v 5
```

### Output
Humann3 outputs for each sample are in their respective sample directories under `humann.dir`.
Humann3 outputs are automatically compressed once they are created. Metaphlan taxa abundances (`<sample>_metaphlan_bugs_list.tsv.gz` are moved out of the temporary direcory created by Humann3 and compressed. Metaphlan taxa abundances are split according by taxonomic levels. Each of the Humann3 outputs for all samples are merged into their respective files `merged_genefamilies.tsv`, `merged_pathabundance.tsv`, `merged_pathcoverage.tsv`, `merged_metaphlan.tsv`.

```
humann.dir/
    |- sample1/
    |- sample2/
    ...
    |- samplen/
        |- samplen_genefamilies.tsv.gz
	|- samplen_pathabundance.tsv.gz
	|- samplen_pathcoverage.tsv.gz
	|- samplen_metaphlan_bugs_list.tsv.gz
	|- samplen_humann_temp.tar.gz
    |- merged_genefamilies.tsv
    |- merged_metaphlan.tsv
    |- merged_metaphlan_class.tsv
    |- merged_metaphlan_family.tsv
    |- merged_metaphlan_genus.tsv
    |- merged_metaphlan_order.tsv
    |- merged_metaphlan_phylum.tsv
    |- merged_metaphlan_species.tsv
    |- merged_pathabundance.tsv
    |- merged_pathcoverage.tsv
```

### Report
Generate a report on humann3 results

```
ocms_shotgun humann3 make build_report
```

</details>