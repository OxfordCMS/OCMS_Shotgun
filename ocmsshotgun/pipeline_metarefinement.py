import sys
import os
from pathlib import Path
from cgatcore import pipeline as P
from ruffus import *
import subprocess
from cgatcore import iotools as IOTools
import glob
import csv
import json
import pysam
import gzip
import re
from Bio import SeqIO
import ocmsshotgun.modules.MetaRefinement as MR
import ocmsshotgun.modules.MetaAssembly as PMA
import ocmsshotgun.modules.BinAssembly as BA

# load options from the config file
PARAMS = P.get_parameters(["pipeline.yml"])
try:
    IOTools.open_file("pipeline.yml")
except FileNotFoundError as e:
    indir = "."
else:
    indir = PARAMS.get('general_input.dir','input.dir')

print("indir =", indir)

###############################################################################
# Run binning_refiner per sample
###############################################################################

@follows(mkdir("01_refined_bins.dir"))
@transform(str(indir) + "/*",
           regex(r".*/([^/]+)$"),
           r"01_refined_bins.dir/\1/done.txt")
def run_binning_refiner(infile, outfile):
    outdir = Path(outfile).parent
    outdir.mkdir(parents=True, exist_ok=True)

    # make absolute paths
    infile = Path(infile).resolve()
    outdir = outdir.resolve()

    sample_name = infile.name

    # run inside outdir using absolute input path
    statement = f"cd {outdir} && Binning_refiner -i {infile} -p {sample_name} -plot"
    P.run(statement)

    Path(outfile).write_text("done\n")

###############################################################################
# Create contig-to-bin mapping dictionaries
###############################################################################

@follows(run_binning_refiner, mkdir("02_contig_bin_maps.dir"))
@transform("01_refined_bins.dir/*/*_Binning_refiner_outputs/*_contigs.txt",
           regex(r".*/([^/]+)_Binning_refiner_outputs/([^/]+)_contigs.txt$"),
           r"02_contig_bin_maps.dir/\2_contig_to_bin.json")
def create_contig_to_bin_map(infile, outfile):
    infile = Path(infile).resolve()
    outfile = Path(outfile).resolve()
    outdir = outfile.parent
    outdir.mkdir(parents=True, exist_ok=True)

    contig_to_bin = {}

    # Read mapping file
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            refined_bin = row["Refined_bin"]
            contigs = [c.strip() for c in row["Contigs"].split(",") if c.strip()]
            for contig in contigs:
                contig_to_bin[contig] = refined_bin

    # Saving as JSON for easy loading later
    with open(outfile, "w") as out_json:
        json.dump(contig_to_bin, out_json, indent=2)

    print(f"[INFO] Created contig→bin map for {len(contig_to_bin)} contigs → {outfile}")

###############################################################################
# Extract read IDs per refined bin
###############################################################################

@follows(create_contig_to_bin_map, mkdir("03_refined_bin_reads.dir"))
@subdivide(PARAMS["general_input_bams_dir"] + "/*_sorted.bam",
           regex(r".*/([^/]+)_sorted\.bam$"),
           add_inputs(r"02_contig_bin_maps.dir/\1_contig_to_bin.json"),
           r"03_refined_bin_reads.dir/\1_*.*")
def extract_refined_bin_reads(inputs, outfiles):
    """Extract read IDs for contigs belonging to refined bins."""
    tool = MR.ExtractRefinedBinReads(inputs, **PARAMS)
    statement = tool.build_statement()
    P.run(
        statement,
        job_memory=PARAMS["extract_refined_bin_reads_memory"],
        job_threads=PARAMS["extract_refined_bin_reads_threads"]
    )

###############################################################################
# Extract fastqs for each bin using reads IDs
###############################################################################

#@follows(extract_refined_bin_reads, mkdir("04_refined_bin_fastqs.dir"))
@transform("03_refined_bin_reads.dir/*/*_read_ids.txt",
           regex(r"03_refined_bin_reads\.dir/([^/]+)/(.+)_([0-9]+)_read_ids\.txt$"),
           [r"04_refined_bin_fastqs.dir/\1/bin\3.fastq.1.gz",
            r"04_refined_bin_fastqs.dir/\1/bin\3.fastq.2.gz"])
def extract_fastqs_by_bin(infile, outfiles):
    tool = MR.FilterFastqByIds(infile, outfiles, **PARAMS)
    statement = tool.build_statement()
    P.run(
        statement,
        job_memory=PARAMS["extract_fastqs_by_bin_memory"],
        job_threads=PARAMS["extract_fastqs_by_bin_threads"]
    )

##############################################################################
# Assemble refined bins
##############################################################################

@follows(extract_fastqs_by_bin, mkdir("05_refined_bin_assemblies.dir"))
@transform(
    "04_refined_bin_fastqs.dir/*/bin*.fastq.1.gz",
    regex(r"04_refined_bin_fastqs\.dir/([^/]+)/bin([0-9]+)\.fastq\.1\.gz$"),
    r"05_refined_bin_assemblies.dir/\1/bin\2.spades.contigs.fasta"
)
def assemble_refined_bins(infile, outfile):

    # Create output directory
    outdir = os.path.dirname(outfile)
    os.makedirs(outdir, exist_ok=True)

    # Import and run bin-level SPAdes wrapper
    import ocmsshotgun.modules.BinAssembly as BA
    assembler = BA.runBinSpades()

    # Build SPAdes run command (auto-detects R2 and singletons)
    statement = assembler.build(infile, outfile, **PARAMS)

    # Execute
    P.run(
        statement,
        job_threads=PARAMS.get("spades_threads", 8),
        job_memory=PARAMS.get("spades_memory", "32G"),
        job_options=PARAMS.get("spades_job_options", "")
    )

@follows(assemble_refined_bins)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

