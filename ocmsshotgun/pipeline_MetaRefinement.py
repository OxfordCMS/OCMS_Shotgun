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
   
    left_reads = os.path.join(PARAMS["binreassembly"]["fastqs"], f"{sample}_1.fastq")
    right_reads = os.path.join(PARAMS["binreassembly"]["fastqs"], f"{sample}_2.fastq")

    threads = PARAMS["binreassembly"]["threads"]
    job_memory = PARAMS["binreassembly"]["job_memory"]
    
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
        "-m {job_memory} "
        "-c {refinement_dir} "
        "> {log_file} 2>&1"
    ).format(**locals())

    P.run(statement,
          job_threads=threads,
          job_memory=job_memory)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    import sys
    sys.exit(main())
