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
import cgatcore.pipeline as P
from ruffus import *

###############################################################################
# Setting up pipeline parameters
###############################################################################
# Load options from the config file
PARAMS = P.get_parameters("pipeline.yml")

FASTAS= glob.glob(os.path.join(PARAMS['fasta_input_dir'], "*.fasta"))

@collate(
    FASTAS,
    regex(r".*/(.+)_corrected\.(megahit|spades)\.contigs\.fasta$"),
    r"metawrap_binning.dir/\1"
)

def runMetawrapBinning(infiles, outfile):
    sample = re.sub(r"_corrected\.(megahit|spades)\.contigs\.fasta$", "", os.path.basename(infiles[0]))

    if PARAMS["is_pooled"]:
        fasta_file = infiles[0]  # only one pooled fasta file expected
        # Correctly ordered fastqs
        fastq_1_files = sorted(glob.glob(os.path.join(PARAMS["fastq_input_dir"], "*_1.fastq")))
        fastq_args_list = []

        for fwd in fastq_1_files:
            rev = re.sub(r'_1\.fastq$', '_2.fastq', fwd)
            if os.path.exists(rev):
                fastq_args_list.extend([fwd, rev])
            else:
                raise FileNotFoundError(f"Missing reverse pair for {fwd}")

        fastq_args = " ".join(fastq_args_list)
        outdir = "metawrap_binning.dir/pooled"
    else:
    
        fasta_file = infiles[0]
        fastq_1 = os.path.join(PARAMS["fastq_input_dir"], f"{sample}_1.fastq")
        fastq_2 = os.path.join(PARAMS["fastq_input_dir"], f"{sample}_2.fastq")
        fastq_args = f"{fastq_1} {fastq_2}"
        outdir = outfile

    
    threads = PARAMS["threads"]
    binning_tools = PARAMS["binning_tools"]
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    out_log = os.path.join(outdir, "metawrap_binning.log")

    print("Creating output dir:", outdir)
    print("Output log path:", out_log)

    statement = ("module purge && "
                  "module load metaWRAP/1.4-20230728-foss-2023a-Python-2.7.18 && "
                  "mkdir -p %(outdir)s && "
                  "metawrap binning "
                  " -o %(outdir)s "
                  " -t %(threads)s "
                  " -a %(fasta_file)s "
                  " %(binning_tools)s "
                  " %(fastq_args)s "
                  " --run-checkm "
                  " &> %(out_log)s "
                  )

    P.run(statement,
          job_threads=threads,
          job_memory=PARAMS["job_memory"],
          fasta_file=fasta_file,
          fastq_args=fastq_args,
          outdir=outdir,
          binning_tools= binning_tools,
          threads=threads,
          out_log=out_log)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    import sys
    sys.exit(main())

