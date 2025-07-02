'''
====================
pipeline_meta_bin.py
====================

:Authors: Jethro Johnson, Uzma Basit Khan

A pipeline that takes contigs from metagenome assemblies and bins them
to generate metagenome-assembled genomes (MAGs).


Overview
========

Multiple tools exist that will take contiguous sequences from metagenome
assemblies and bin them based on coverage and other features.

Typically, tools will take as input either a pooled assembly (where
sequences from multiple samples have been pooled prior to assembly), or
a individual assemblies (where each sample is assembled separately).
While the latter situation isn't optimal, it's sometimes preferred due
to memory restrictions on large metagenome assemblies.

This pipeline is intended to cope with either scenario. It will perform
the following steps:

i.    Create a Bowtie2 index for metagenome assemblies in the
      input_metagenome.dir directory
ii.   Map fastq files to the reference index. If individual assemblies
      have been performed, then this will be a 1:1 mapping, based on
      matching sample names. If a pooled assembly has been performed,
      then input fastq files will be individually mapped to the same
      reference.
iii.  Mapped bamfiles + assembled metagenome files will then be passed
      to different binning tools to generate MAGs


Requirements
============

Input sequence data is expected to be in gzip compressed fastq format.
Samples can be single-end, paired-end, or paired-end plus singletons.
However, the naming convention is expected to be <sample-name>.fastq.1.gz,
<sample-name>.fastq.2.gz, <sample-name>.fastq.3.gz for read1, read2, and
singletons, respectively.

Mapping of input fastq files to their respective metagenome assembly is
done on the basis of reqular expressions specified in the pipeline.yml
config file. Input metagenome files must therefore be named according
to the regular expression provided. 


'''

import re
import os
import glob
import sys
from cgatcore import pipeline as P
import ocmstoolkit.modules.Utility as Utility
from ruffus import *

PARAMS = P.get_parameters("pipeline.yml")

# Access nested YAML entries under 'general'
fasta_indir = PARAMS["general"]["fasta.dir"]
fastq_indir = PARAMS["general"]["fastq.dir"]

#check all files to be processed
FASTQs = Utility.get_fastns(fastq_indir)

# List all FASTA files in the input directory
fasta_files = glob.glob(os.path.join(fasta_indir, "*.fasta"))

###############################################################################
# Index Fasta files (individual or pooled)
###############################################################################

@transform(fasta_files,
           regex(r'(.+)/(.+).fasta'),
           r'\1/\2_index/\2.1.bt2')
def indexfasta(infile, outfile):
    """
    Index metagenome assembly FASTA files using Bowtie2.
    Outputs index files to <sample>_index/ within the same base directory.
    """
    out_dir = re.sub(r'\.fasta$', '_index', infile)
    os.makedirs(out_dir, exist_ok=True)

    threads = PARAMS.get("bowtie2_indexing", {}).get("threads", 1)
    out_prefix = os.path.join(out_dir, re.sub(r'\.fasta$', '', os.path.basename(infile)))

    statement = (
        "bowtie2-build "
        " --threads %(threads)s "
        " %(infile)s "
        " %(out_prefix)s > %(out_prefix)s.log"
    )

    P.run(statement,
          infile=infile,
          out_prefix=out_prefix,
          threads=threads)

@follows(indexfasta)
@transform(FASTQs,
           regex(r'.*/(.+)\.fastq\.1\.gz'),
           r'mapping.dir/\1_depth.txt')
def mapfastq2fasta(infile, outfile):
    """
    Map FASTQ files to a reference: individual or pooled.
    - If 'pooled_fasta' is False → use sample-specific FASTA
    - If 'pooled_fasta' is True:
        - If only one *.fasta in fasta.dir → map all samples to it
        - If multiple *.fasta → map to group-specific FASTA based on sample prefix
    Each mapping is performed individually per sample.
    """
    os.makedirs("mapping.dir", exist_ok=True)

    sample = os.path.basename(infile).split(".fastq")[0]
    fasta_dir = PARAMS["general"]["fasta.dir"]
    use_pooled = PARAMS.get("bowtie2_mapping", {}).get("pooled_fasta", False)

    if use_pooled:
        # List all FASTA files to determine if this is single or grouped pooled mapping
        fasta_files = glob.glob(os.path.join(fasta_dir, "*.fasta"))

        if len(fasta_files) == 1:
            # Single pooled FASTA for all samples (no grouping)
            fasta = fasta_files[0]
        else:
            # Grouped pooled mapping — requires user-defined regex
            grouping_regex = PARAMS["bowtie2_mapping"].get("grouping_regex")
            if not grouping_regex:
                raise ValueError(
                    "Multiple FASTA files found in pooled mode, but no 'grouping_regex' defined in the YAML. "
                    "Either provide a regex or reduce to a single pooled FASTA."
                )
            match = re.search(grouping_regex, sample)
            if not match:
                raise ValueError(
                    f"Could not extract pooling group from sample '{sample}' using regex '{grouping_regex}'."
                )
            group = match.group(1)
            matches = glob.glob(os.path.join(fasta_dir, f"*{group}*.fasta"))  # match like pooled_chow.fasta or chow_pooled.fasta
            if len(matches) != 1:
                raise FileNotFoundError(
                    f"Could not uniquely find group-level pooled FASTA for group '{group}' in {fasta_dir}. "
                    f"Expected 1 match for '*{group}*.fasta', found {len(matches)}: {matches}"
                )
            fasta = matches[0]
    else:
        # Individual mapping — match based on sample name prefix
        matches = glob.glob(os.path.join(fasta_dir, f"{sample}*.fasta"))
        if len(matches) != 1:
            raise FileNotFoundError(
                f"Could not uniquely find FASTA for sample '{sample}' in {fasta_dir}. "
                f"Expected 1 match, found {len(matches)}: {matches}"
            )
        fasta = matches[0]

    # Determine correct Bowtie2 index path
    fasta_base = os.path.splitext(os.path.basename(fasta))[0]
    index_dir = os.path.join(fasta_dir, f"{fasta_base}_index")
    index = os.path.join(index_dir, fasta_base)

    # Input FASTQ pairs
    fastq_1 = infile
    fastq_2 = infile.replace(".1.gz", ".2.gz")

    # Output files
    bam = f"mapping.dir/{sample}.bam"
    sorted_bam = f"mapping.dir/{sample}_sorted.bam"
    depth = outfile
    threads = PARAMS.get("bowtie2_indexing", {}).get("threads", 1)

    # Run mapping
    statement = (
        "bowtie2 --threads %(threads)s -x %(index)s "
        "-1 %(fastq_1)s -2 %(fastq_2)s | "
        "samtools view -bS - > %(bam)s && "
        "samtools sort -o %(sorted_bam)s %(bam)s && "
        "jgi_summarize_bam_contig_depths --outputDepth %(depth)s %(sorted_bam)s"
    )

    P.run(statement,
          fastq_1=fastq_1,
          fastq_2=fastq_2,
          index=index,
          bam=bam,
          sorted_bam=sorted_bam,
          depth=depth,
          threads=threads)
    
def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    import sys
    sys.exit(main())

