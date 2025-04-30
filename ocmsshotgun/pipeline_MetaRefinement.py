import os
import glob
import cgatcore.pipeline as P
from ruffus import *
import random

# Load config
PARAMS = P.get_parameters("pipeline.yml")

# Define binning inputs from previous step
BINNING_OUTPUTS = glob.glob(os.path.join(PARAMS['bins_dir'],'*'))

@collate(
    BINNING_OUTPUTS,
    regex(rf"{PARAMS['bins_dir']}/(.+)"),
    r"metawrap_bin_refinement.dir/\1"
)
def runMetawrapBinRefinement(infiles, outfile):
    sample = os.path.basename(outfile)

    bins_base = os.path.join(PARAMS["bins_dir"], sample)
    outdir = f"metawrap_bin_refinement.dir/{sample}"
        
    if PARAMS["is_pooled"]:
        # Run once for pooled case only
        first_sample = os.path.basename(infiles[0])
        if sample != first_sample:
            return

    metabat_dir = os.path.join(bins_base, "metabat2_bins")
    maxbin_dir = os.path.join(bins_base, "maxbin2_bins")
    concoct_dir = os.path.join(bins_base, "concoct_bins")

    threads = PARAMS["binrefinement"]["threads"]
    job_memory = PARAMS["binrefinement"]["job_memory"]
    completeness = PARAMS["binrefinement"]["completeness"]
    contamination = PARAMS["binrefinement"]["contamination"]
   
    log_file = os.path.join(outdir, "bin_refinement.log")

    statement = (
        "module purge && "
        "module load metaWRAP/1.4-20230728-foss-2023a-Python-2.7.18 && "
        "mkdir -p {outdir} && "
        "metawrap bin_refinement "
        "-o {outdir} "
        "-t {threads} "
        "-A {metabat_dir} "
        "-B {maxbin_dir} "
        "-C {concoct_dir} "
        "-c {completeness} -x {contamination} "
        "> {log_file} 2>&1"
        ).format(**locals())

    P.run(statement,
          job_threads=threads,
          job_memory=job_memory)

# Reassemble bins after refinement
@transform(
    runMetawrapBinRefinement,
    regex(r"metawrap_bin_refinement.dir/(.+)"),
    r"metawrap_reassembled_bins.dir/\1"
)
def runMetawrapReassembleBins(infiles,outfile):
    sample = os.path.basename(outfile)

    # Dynamically find directory like metawrap_50_10_bins (thats is one of the outputs from BinRefienemnt)
    refinement_base = os.path.join("metawrap_bin_refinement.dir", sample)
    refinement_dirs = glob.glob(os.path.join(refinement_base,"metawrap_*_bins"))
    
    if not refinement_dirs:
        raise FileNotFoundError(f"No metawrap_*_bins directory found in {refinement_base}")

    refinement_dir = refinement_dirs[0]  # Take the first (assumed only) match
    reassembly_outdir = os.path.join("metawrap_reassembled_bins.dir", sample)
   
    fastq_dir = PARAMS["binreassembly"]["fastqs"]

    if PARAMS.get("is_pooled", False):
        left_candidates = sorted(glob.glob(os.path.join(fastq_dir, "*.fastq.1.gz")))
        right_candidates = sorted(glob.glob(os.path.join(fastq_dir, "*.fastq.2.gz")))

        if not left_candidates or not right_candidates:
            raise FileNotFoundError(f"Could not find *.fastq.1.gz or *.fastq.2.gz in pooled input directory: {fastq_dir}")

        tmp_fastq_dir = os.path.join("tmp_pooled_fastqs", sample)
        os.makedirs(tmp_fastq_dir, exist_ok=True)

        # Process left and right reads
        def decompress_and_rename(source, suffix):
            basename = os.path.basename(source).replace(".fastq." + suffix + ".gz", "")
            target = os.path.abspath(os.path.join(tmp_fastq_dir, f"{basename}_{suffix}.fastq"))

            print(f"Decompressing: {source} → {target}")
            if os.path.exists(target):
                os.remove(target)

            ret = os.system(f"zcat {source} > {target}")
            if ret != 0:
                raise RuntimeError(f"Failed to decompress {source}")
            return target

        left_reads = decompress_and_rename(left_candidates[0], "1")
        right_reads = decompress_and_rename(right_candidates[0], "2")

    else:
        # Use sample-named FASTQs for unpooled data
        left_reads = os.path.join(PARAMS["binreassembly"]["fastqs"], f"{sample}_1.fastq")
        right_reads = os.path.join(PARAMS["binreassembly"]["fastqs"], f"{sample}_2.fastq")

    threads = PARAMS["binreassembly"]["threads"]
    job_memory_raw = PARAMS["binreassembly"]["job_memory"]

    # Extract numeric part of memory (e.g., 20)
    memory_per_cpu = int(''.join(filter(str.isdigit, job_memory_raw)))

    # Now calculate total memory for SPAdes
    total_memory_spades = memory_per_cpu * threads  # e.g., 20 * 4 = 80

    # SLURM still uses raw memory per CPU (like "20G")
    job_memory_slurm = job_memory_raw
    
    log_file = os.path.join(reassembly_outdir, "reassemble_bins.log")
    
    statement = (
        "module purge && "
        "module load metaWRAP/1.4-20230728-foss-2023a-Python-2.7.18 && "
        "mkdir -p {reassembly_outdir} && "
        "metawrap reassemble_bins "
        "-o {reassembly_outdir} "
        "-1 {left_reads} "
        "-2 {right_reads} "
        "-t {threads} "
        "-m {total_memory_spades} "
        "-b {refinement_dir} "
        "> {log_file} 2>&1"
    ).format(**locals())

    P.run(statement,
          job_threads=threads,
          job_memory=job_memory_slurm)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    import sys
    sys.exit(main())
