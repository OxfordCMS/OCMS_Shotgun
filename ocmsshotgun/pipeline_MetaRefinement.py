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

    threads = PARAMS["binrefinement"]["threads"]
    job_memory = PARAMS["binrefinement"]["job_memory"]
    completeness = PARAMS["binrefinement"]["completeness"]
    contamination = PARAMS["binrefinement"]["contamination"]
   
    log_file = os.path.join(outdir, "bin_refinement.log")

    # User-specified binning tools
    binning_tools = PARAMS.get("binning_tools", [])

    # Map tool names to subdirectory names
    tool_subdirs = {
        "metabat2": "metabat2_bins",
        "maxbin2": "maxbin2_bins",
        "concoct": "concoct_bins",
        "vamb": "vamb_bins"
        # Add more tools and subdirs as needed
    }

    # Validate number of tools
    if not (2 <= len(binning_tools) <= 3):
        raise ValueError("MetaWRAP bin_refinement requires 2 or 3 binning tool outputs.")

    # Fixed flag order
    metawrap_flags = ["-A", "-B", "-C"]
    bin_args = []

    for flag, tool in zip(metawrap_flags, binning_tools):
        if tool not in tool_subdirs:
            raise ValueError(f"Unknown binning tool: {tool}")
        subdir = tool_subdirs[tool]
        bin_dir = os.path.join(bins_base, subdir)
        if not os.path.isdir(bin_dir):
            raise FileNotFoundError(f"Bin directory for tool '{tool}' not found: {bin_dir}")
        bin_args.append(f"{flag} {bin_dir}")

    bin_args_str = " ".join(bin_args)

    statement = (
        "module purge && "
        "module load metaWRAP/1.4-20230728-foss-2023a-Python-2.7.18 && "
        "mkdir -p {outdir} && "
        "metawrap bin_refinement "
        "-o {outdir} "
        "-t {threads} "
        "{bin_args_str} "
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
def runMetawrapReassembleBins(infiles, outfile):
    sample = os.path.basename(outfile)

    # Get parameters from YAML
    completeness = PARAMS["bin_refinement"]["completeness"]
    contamination = PARAMS["bin_refinement"]["contamination"]

    # Build expected directory name
    expected_dirname = f"metawrap_{completeness}_{contamination}_bins"

    # Full path to the output directory from bin_refinement
    refinement_dir = os.path.join("metawrap_bin_refinement.dir", sample, expected_dirname)

    # Check if the directory exists
    if not os.path.exists(refinement_dir):
        raise FileNotFoundError(f"Expected refinement directory not found: {refinement_dir}")

    # Define reassembly output directory
    reassembly_outdir = os.path.join("metawrap_reassembled_bins.dir", sample)

    # Get fastq directory from config
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
