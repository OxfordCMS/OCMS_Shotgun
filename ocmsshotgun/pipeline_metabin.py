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

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    import sys
    sys.exit(main())

