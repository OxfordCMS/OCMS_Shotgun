import sys
import os
from pathlib import Path
from cgatcore import pipeline as P
from ruffus import *
import subprocess
from cgatcore import iotools as IOTools
from cgatcore import pipeline as P
import glob
import csv
import json
import pysam
import gzip
from Bio import SeqIO

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

@follows(create_contig_to_bin_map, mkdir("03_refined_bin_reads.dir"))
@subdivide("input_bams.dir/*_sorted.bam",
           regex(r".*/([^/]+)_sorted\.bam$"),
           add_inputs(r"02_contig_bin_maps.dir/\1_contig_to_bin.json"),
           r"03_refined_bin_reads.dir/\1_*.txt")
def extract_refined_bin_reads(inputs, outfiles):
    """
    Extract read IDs for contigs belonging to refined bins.
    Inputs:
      - Sorted BAM file (with alignments)
      - JSON mapping file (contig → bin)
    Outputs:
      - One text file per refined bin containing read IDs
    """

    infile, map_file = map(Path, inputs)
    sample_id = infile.stem.replace("_sorted", "")
    outdir = Path("03_refined_bin_reads.dir") / sample_id
    outdir.mkdir(parents=True, exist_ok=True)

    # --- Load contig → bin mapping ---
    if not map_file.exists():
        raise FileNotFoundError(f"Mapping file not found for {sample_id}: {map_file}")

    with open(map_file) as f:
        contig_to_bin = json.load(f)

    # --- Group contigs by bin ---
    bin_to_contigs = {}
    for contig, bin_name in contig_to_bin.items():
        bin_to_contigs.setdefault(bin_name, []).append(contig)

    # --- Prepare output files ---
    bin_to_handles = {}
    for bin_name in bin_to_contigs:
        out_path = outdir / f"{bin_name}_read_ids.txt"
        bin_to_handles[bin_name] = open(out_path, "w")

    valid_contigs = set(contig_to_bin.keys())

    # --- Iterate over BAM and write read IDs ---
    with pysam.AlignmentFile(infile, "rb") as bam_in:
        for read in bam_in.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            ref_name = bam_in.get_reference_name(read.reference_id)
            if ref_name not in valid_contigs:
                continue
            bin_name = contig_to_bin[ref_name]
            handle = bin_to_handles.get(bin_name)
            if handle:
                handle.write(read.query_name + "\n")

    # --- Close output handles ---
    for handle in bin_to_handles.values():
        handle.close()

    # --- Sanity check for missing contigs ---
    with pysam.AlignmentFile(infile, "rb") as bam_in:
        bam_contigs = set(bam_in.references)
    missing_contigs = [c for c in contig_to_bin if c not in bam_contigs]

    if missing_contigs:
        print(f"[WARN] {sample_id}: {len(missing_contigs)} refined-bin contigs not present in BAM header")
    else:
        print(f"[OK] {sample_id}: all refined-bin contigs found in BAM")

    print(f"[{sample_id}] Extracted read IDs for {len(bin_to_contigs)} refined bins.")

#Filter original FASTQ files using read ID lists
@follows(extract_refined_bin_reads, mkdir("04_refined_bin_fastqs.dir"))
@transform("03_refined_bin_reads.dir/*/*_read_ids.txt",
           regex(r"03_refined_bin_reads\.dir/(.+)/(.+)_read_ids\.txt"),
           [r"04_refined_bin_fastqs.dir/\1/\2_R1.fastq.gz",
            r"04_refined_bin_fastqs.dir/\1/\2_R2.fastq.gz"])
def extract_fastqs_by_bin(read_id_file, outfiles):
    """
    Filter original FASTQ files using read ID lists from the previous step.
    Produces paired-end FASTQs for each refined bin.
    """

    out_r1, out_r2 = map(Path, outfiles)
    out_r1.parent.mkdir(parents=True, exist_ok=True)

    # Derive sample ID
    sample_id = out_r1.parent.name
    print(f"[{sample_id}] Processing refined bin: {out_r1.stem.replace('_R1', '')}")

    # --- Load read IDs ---
    with open(read_id_file) as f:
        read_ids = set(line.strip() for line in f if line.strip())
    print(f"  Loaded {len(read_ids)} read IDs")

    if not read_ids:
        print(f"  [WARN] No read IDs found, skipping bin")
        return

    # --- Locate original FASTQs ---
    fq_dir = Path("input_fastqs.dir")
    r1_path = fq_dir / f"{sample_id}.fastq.1.gz"
    r2_path = fq_dir / f"{sample_id}.fastq.2.gz"

    if not r1_path.exists() or not r2_path.exists():
        raise FileNotFoundError(f"FASTQ files for {sample_id} not found in {fq_dir}")

    # --- Define helper to filter FASTQ ---
    def filter_fastq(in_path, out_path):
        count_in, count_out = 0, 0
        with gzip.open(in_path, "rt") as in_fq, gzip.open(out_path, "wt") as out_fq:
            for rec in SeqIO.parse(in_fq, "fastq"):
                count_in += 1
                if rec.id in read_ids:
                    SeqIO.write(rec, out_fq, "fastq")
                    count_out += 1
        return count_in, count_out

    # --- Filter R1 and R2 ---
    c1_in, c1_out = filter_fastq(r1_path, out_r1)
    c2_in, c2_out = filter_fastq(r2_path, out_r2)

    print(f"  R1: {c1_out}/{c1_in} reads written")
    print(f"  R2: {c2_out}/{c2_in} reads written")
    print(f"  [OK] Wrote FASTQs: {out_r1.name}, {out_r2.name}")



@follows(extract_fastqs_by_bin)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

