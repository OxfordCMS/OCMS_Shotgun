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
import sys
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

# ---------------------------------------------------------
# 1. Index Fasta files (individual or pooled)
# ---------------------------------------------------------
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

###############################################################################
# 2. Prepare depth/coverage files for all binning tools
###############################################################################

# ---------------------------------------------------------
# a. Map FASTQs to FASTAs and generate MetaBAT2 depth
# ---------------------------------------------------------
@follows(indexfasta)
@collate(FASTQs,
         # Regular expression for the query fastq
         regex('.+/' + PARAMS['mapfastq2fasta_fastq_regex']),
         # Regular expression for the reference fasta
         add_inputs(os.path.join(PARAMS['general_fasta_dir'],
                                 PARAMS['mapfastq2fasta_fasta_regex'])),
         # Regular expression for the output file
#         os.path.join('01_mapping.dir',
 #                     PARAMS['mapfastq2fasta_output_regex']+ ".done"))
         r"01_mapping.dir/" + PARAMS["mapfastq2fasta"]["output_regex"] + ".done")
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

    #Output files
    bam = os.path.join("01_mapping.dir", sample + ".bam")
    sorted_bam = os.path.join("01_mapping.dir", sample + "_sorted.bam")
    depth = os.path.join("01_mapping.dir", sample + "_metabat2_depth.txt")
    threads = PARAMS["mapfastq2fasta"]["job_threads"]
    logfile = os.path.join("01_mapping.dir", sample + "_mapping.log")
    
    statement = (
        "(bowtie2 --threads %(threads)s -x %(index_prefix)s "
        "-1 %(fastq_1)s -2 %(fastq_2)s | "
        "samtools view -bS - > %(bam)s && "
        "samtools sort -o %(sorted_bam)s %(bam)s && "
        "samtools index %(sorted_bam)s && "
        "jgi_summarize_bam_contig_depths --outputDepth %(depth)s %(sorted_bam)s) "
        " &> %(logfile)s"
    )

    P.run(statement,
          job_memory=PARAMS["mapfastq2fasta_job_memory"],
          job_threads=PARAMS["mapfastq2fasta_job_threads"])
    
    # Mark task complete
    Path(outfile).touch()

@follows(mapfastq2fasta)
@originate("01_mapping.dir/cumulative_metabat2_depth.txt")
def generate_cumulative_metabat2_depth(outfile):
    """
    Generate cumulative MetaBAT2 depth file (pooled mode only).
    """
    if PARAMS["mapfastq2fasta"]["mapping_mode"] != "many2one":
        Path(outfile).touch()
        return

    bam_files = glob.glob("01_mapping.dir/*_sorted.bam")
    if len(bam_files) <= 1:
        raise ValueError("Expected multiple BAM files for pooled samples.")

    bam_inputs = " ".join(bam_files)
    statement = "jgi_summarize_bam_contig_depths --outputDepth %(outfile)s %(bam_inputs)s"
    P.run(statement)
    Path(outfile).touch()

# ---------------------------------------------------------
# b. Generate MaxBin2 depth (per-sample OR pooled)
# ---------------------------------------------------------
@follows(mapfastq2fasta)
@transform("01_mapping.dir/*_metabat2_depth.txt",
           regex(r"01_mapping.dir/(?!cumulative)(.+)_metabat2_depth.txt"),
           r"01_mapping.dir/\1_maxbin2_depth.txt")
def generate_maxbin2_depth(infile, outfile):
    print("DEBUG >>> infiles:", infile)
    """
    Create MaxBin2-compatible depth file from MetaBAT2 depth,
    keeping only (contigName, contigLen, totalAvgDepth).
    """
    statement = "cut -f 1,2,3 %(infile)s > %(outfile)s"
    P.run(statement)

@follows(generate_maxbin2_depth)
@collate("01_mapping.dir/*_maxbin2_depth.txt",
         regex(r"01_mapping.dir/.+_maxbin2_depth.txt"),
         "01_mapping.dir/cumulative_maxbin2_depth.txt")
def generate_cumulative_maxbin2_depth(infiles, outfile):
    """
    Merge per-sample MaxBin2 depth files into a single cumulative file.
    Columns: contigName, contigLen, depth_sample1, depth_sample2, ...
    """
    if PARAMS["mapfastq2fasta"]["mapping_mode"] != "many2one":
        Path(outfile).touch()
        return

    input_files = sorted(infiles)
    
    sample_names = [
        os.path.basename(f).replace("_maxbin2_depth.txt", "")
        for f in input_files
    ]

    first = f"<(cut -f1,2,3 {input_files[0]})"
    others = " ".join([f"<(cut -f3 {f})" for f in input_files[1:]])

    if others:
        paste_cmd = f"paste {first} {others}"
    else:
        paste_cmd = f"cat {input_files[0]}"

    header = "contigName\tcontigLen\t" + "\t".join(sample_names)
    statement = f"(echo -e '{header}' && {paste_cmd} | tail -n +2) > {outfile}"

    P.run(f"bash -c \"{statement}\"")

# ---------------------------------------------------------
# c. Prepare all binning input files
# ---------------------------------------------------------
@follows(
    generate_maxbin2_depth,
    generate_cumulative_maxbin2_depth,
    generate_cumulative_metabat2_depth,
    mapfastq2fasta
)
def prepare_binning_inputs():
    """
    Prepares all depth/coverage files required for MetaBAT2, MaxBin2, and CONCOCT.
    """
    pass

# -------------------------------------------------------------------------
# User-selected binning tools
# -------------------------------------------------------------------------
selected_tools = [
    t.strip() for t in PARAMS.get("binning_tools", "metabat2").split(",")
]

# ---------------------------------------------------------
# 3. Run Binning Tools
# ---------------------------------------------------------
# -------------------------------------------------------------------------
# MetaBAT2 execution via DepthFileManager
# -------------------------------------------------------------------------

mapping_mode = PARAMS["mapfastq2fasta"]["mapping_mode"]

if mapping_mode == "many2one" and "metabat2" in selected_tools:

    @follows(prepare_binning_inputs)
    @files("01_mapping.dir/cumulative_metabat2_depth.txt",
           "02_bins.dir/pooled/pooled_metabat2_done.txt")
    def run_metabat2(infile, outfile):
        pooled_dir = "02_bins.dir/pooled/metabat2_bins"
        os.makedirs(pooled_dir, exist_ok=True)

        # Fetch depth files for MetaBAT2
        depth_manager = MB.DepthFileManager("01_mapping.dir")
        depth_files = depth_manager.get_depth_files("many2one", tool="metabat2")
        assert len(depth_files) == 1, f"Expected 1 pooled depth file, got {len(depth_files)}"

        commands = MB.MetaBAT2Runner.run_all("02_bins.dir", PARAMS, tool="metabat2")
        assert len(commands) == 1, f"Expected 1 pooled command, got {len(commands)}"
        _, statement = commands[0]

        P.run(statement,
              job_memory=PARAMS["binners_job_memory"],
              job_threads=PARAMS["binners_job_threads"])
        Path(outfile).touch()

elif mapping_mode == "one2one" and "metabat2" in selected_tools:
    
    @follows(prepare_binning_inputs)
    @subdivide(f"01_mapping.dir/*_metabat2_depth.txt",
               regex(rf"01_mapping\.dir/(?!cumulative)(.+)_metabat2_depth\.txt"),
               r"02_bins.dir/\1/\1_metabat2_done.txt")
    def run_metabat2(infile, outfile, tool="metabat2"):
        sample = os.path.basename(infile).replace("_metabat2_depth.txt", "")
        sample_dir = f"02_bins.dir/{sample}/{tool}_bins"
        os.makedirs(sample_dir, exist_ok=True)

        # Fetch depth files (per-sample)
        depth_manager = MB.DepthFileManager("01_mapping.dir")
        depth_files = depth_manager.get_depth_files("one2one", tool=tool)

        # Build all commands
        commands = MB.MetaBAT2Runner.run_all("02_bins.dir", PARAMS, tool=tool)
        command_map = dict(commands)

        if sample not in command_map:
            raise RuntimeError(f"No {tool} command built for {sample}")

        statement = command_map[sample]

        P.run(statement,
              job_memory=PARAMS["binners_job_memory"],
              job_threads=PARAMS["binners_job_threads"])
        Path(outfile).touch()

# -------------------------------------------------------------------------
# MaxBin2 execution via DepthFileManager
# -------------------------------------------------------------------------

mapping_mode = PARAMS["mapfastq2fasta"]["mapping_mode"]

if mapping_mode == "many2one" and "maxbin2" in selected_tools:

    @follows(prepare_binning_inputs)
    @files("01_mapping.dir/cumulative_maxbin2_depth.txt",
           "02_bins.dir/pooled/pooled_maxbin2_done.txt")
    def run_maxbin2(infile, outfile):
        pooled_dir = "02_bins.dir/pooled/maxbin2_bins"
        os.makedirs(pooled_dir, exist_ok=True)

        # Fetch depth files for MaxBin2
        depth_manager = MB.DepthFileManager("01_mapping.dir")
        depth_files = depth_manager.get_depth_files("many2one", tool="maxbin2")
        assert len(depth_files) == 1, f"Expected 1 pooled depth file, got {len(depth_files)}"

        # Build pooled command
        commands = MB.MaxBin2Runner.run_all("02_bins.dir", PARAMS, tool="maxbin2")
        assert len(commands) == 1, f"Expected 1 pooled command, got {len(commands)}"
        _, statement = commands[0]

        P.run(statement,
              job_memory=PARAMS["binners_job_memory"],
              job_threads=PARAMS["binners_job_threads"])
        Path(outfile).touch()

elif mapping_mode == "one2one" and "maxbin2" in selected_tools:
    @follows(prepare_binning_inputs)
    @subdivide(f"01_mapping.dir/*_maxbin2_depth.txt",
               regex(rf"01_mapping\.dir/(?!cumulative)(.+)_maxbin2_depth\.txt"),
               r"02_bins.dir/\1/\1_maxbin2_done.txt")
    def run_maxbin2(infile, outfile, tool ="maxbin2"):
        sample = os.path.basename(infile).replace("_maxbin2_depth.txt", "")
        sample_dir = f"02_bins.dir/{sample}/{tool}_bins"
        os.makedirs(sample_dir, exist_ok=True)

        # Fetch depth files (per-sample)
        depth_manager = MB.DepthFileManager("01_mapping.dir")
        depth_files = depth_manager.get_depth_files("one2one", tool=tool)

        # Build all commands
        commands = MB.MaxBin2Runner.run_all("02_bins.dir", PARAMS, tool=tool)
        command_map = dict(commands)

        if sample not in command_map:
            raise RuntimeError(f"No {tool} command built for {sample}")

        statement = command_map[sample]

        P.run(statement,
              job_memory=PARAMS["binners_job_memory"],
              job_threads=PARAMS["binners_job_threads"])
        Path(outfile).touch()

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(main())

