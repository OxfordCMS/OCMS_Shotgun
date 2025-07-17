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
import ocmsshotgun.modules.MetaBin as MB
from pathlib import Path

PARAMS = P.get_parameters("pipeline.yml")

# Access nested YAML entries under 'general'
fasta_indir = PARAMS["general"]["fasta_dir"]
fastq_indir = PARAMS["general"]["fastq_dir"]

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

    threads = PARAMS["mapfastq2fasta"]["threads"]
    out_prefix = os.path.join(out_dir,
                              re.sub(r'\.fasta$', '', os.path.basename(infile)))

    statement = ("bowtie2-build "
                 " --threads %(threads)s "
                 " %(infile)s "
                 " %(out_prefix)s > %(out_prefix)s.log"
                )

    P.run(statement,
          infile=infile,
          out_prefix=out_prefix,
          threads=threads)

@follows(indexfasta)
@collate(FASTQs,
         # Regular expression for the query fastq
         regex('.+/' + PARAMS['mapfastq2fasta_fastq_regex']),
         # Regular expression for the reference fasta
         add_inputs(os.path.join(PARAMS['general_fasta_dir'],
                                 PARAMS['mapfastq2fasta_fasta_regex'])),
         # Regular expression for the output file
         os.path.join('01_mapping.dir',
                      PARAMS['mapfastq2fasta_output_regex']))
def mapfastq2fasta(infiles, outfile):
    """
    Map FASTQ files to reference FASTA based on regex-defined logic.
    The mapping can be:
      - 1-to-1: sample to its own FASTA
      - many-to-1: multiple samples to a shared pooled FASTA
    The FASTA to map to is resolved using `fastq_regex` and `fasta_regex` defined in the YAML.
    """
    
    # Ensure output directory exists
    os.makedirs("01_mapping.dir", exist_ok=True)

    # Unpack the tuple: (fastq_path, index_path)
    fastq_path, index_path = infiles[0]

    # Derive second read
    fastq_1 = fastq_path
    fastq_2 = fastq_path.replace(".1.gz", ".2.gz")
    sample = os.path.basename(fastq_path).split(".fastq")[0]

    # Bowtie2 index prefix (remove .1.bt2)
    index_prefix = index_path.replace(".1.bt2", "")

    # Output files
    bam = f"01_mapping.dir/{sample}.bam"
    sorted_bam = f"01_mapping.dir/{sample}_sorted.bam"
    depth = outfile
    threads = PARAMS["mapfastq2fasta"]["threads"]
    logfile = f"01_mapping.dir/{sample}_mapping.log"
    
    statement = (
                 "(bowtie2 --threads %(threads)s -x %(index)s "
                 "-1 %(fastq_1)s -2 %(fastq_2)s | "
                 "samtools view -bS - > %(bam)s && "
                 "samtools sort -o %(sorted_bam)s %(bam)s && "
                 "jgi_summarize_bam_contig_depths --outputDepth %(depth)s %(sorted_bam)s) "
                 " &> %(logfile)s"
                )

    P.run(statement,
          fastq_1=fastq_1,
          fastq_2=fastq_2,
          index=index_prefix,
          bam=bam,
          sorted_bam=sorted_bam,
          depth=depth,
          threads=threads,
          logfile = logfile)

@follows(mapfastq2fasta)
@active_if(PARAMS["mapfastq2fasta"]["mapping_mode"] == "many2one")
@originate("01_mapping.dir/cumulative_depth.txt")
def generate_cumulative_depth(outfile):
    bam_files = glob.glob("01_mapping.dir/*_sorted.bam")
    assert len(bam_files) > 1, "Expected multiple BAM files for pooled samples."

    bam_inputs = " ".join(bam_files)
    statement = f'''
    jgi_summarize_bam_contig_depths --outputDepth {outfile} {bam_inputs};
    touch {outfile}
    '''
    P.run(statement)

@transform(mapfastq2fasta,
           regex(r"01_mapping.dir/(.+)_depth.txt"),
           add_inputs(os.path.join(PARAMS['general']['fasta_dir'], r"\1.fasta")),
           r"03_bins.dir/\1_bins")
@active_if(PARAMS["mapfastq2fasta"]["mapping_mode"] == "one2one")
def run_metabat2_unpooled(infiles, outfile):
    depth_file, assembly = infiles
    sample = os.path.basename(depth_file).replace("_depth.txt", "")
    
    # Define output directory and bin file prefix
    output_dir = outfile  # e.g., '03_bins/<sample>_bins.dir'
    os.makedirs(output_dir, exist_ok=True)

    prefix = os.path.basename(outfile).replace("_bins", "")
    output_prefix = os.path.join(output_dir, prefix)

    binner = MB.MetaBAT2Runner(
        assembly=assembly,
        depth_file=depth_file,
        prefix=prefix,
        output_dir=output_dir,
        **PARAMS["metabat2"]
    )

    statement = binner.build_command()
    P.run(statement)
    Path(outfile).touch()

@merge(os.path.join("01_mapping.dir", "cumulative_depth.txt"),
       "03_bins.dir/pooled_bins_done.txt")
@active_if(PARAMS["mapfastq2fasta"]["mapping_mode"] == "many2one")
def run_metabat2_pooled(depth_file, outfile):

    # Get pooled fasta
    fasta_dir = PARAMS["general"]["fasta_dir"]
    fasta_pattern = glob.glob(os.path.join(fasta_dir, "*.fasta"))
    assert len(fasta_pattern) == 1, f"Expected one pooled fasta, found: {fasta_pattern}"
    pooled_fasta = fasta_pattern[0]

    # Setup output
    prefix = os.path.splitext(os.path.basename(pooled_fasta))[0]
    output_dir = "03_bins.dir"
    os.makedirs(output_dir, exist_ok=True)

    # Run MetaBAT2
    binner = MB.MetaBAT2Runner(
        assembly=pooled_fasta,
        depth_file=depth_file,
        prefix=prefix,
        output_dir=output_dir,
        **PARAMS["metabat2"]
    )

    statement = binner.build_command()
    P.run(statement)

    # Mark completion
    Path(outfile).touch()

def main(argv=None):

    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    import sys
    sys.exit(main())

