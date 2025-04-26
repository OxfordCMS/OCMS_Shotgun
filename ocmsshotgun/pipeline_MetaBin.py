#import modules
import gzip
import shutil
import os
import glob
import ocmsshotgun.modules.Utility as utility
import cgatcore.pipeline as P
from ruffus import *

# Load options from the config file
PARAMS = P.get_parameters("pipeline.yml")

FASTQS = glob.glob(os.path.join(PARAMS['fastq_input_dir'], "*.fastq.*.gz"))

# Decompress and rename to *_1.fastq or *_2.fastq
@transform(
    FASTQS,
    regex(r".*/(.+)\.fastq\.(\d)\.gz"),
    r"fastq_prepped_for_metawrap.dir/\1_\2.fastq"
)
def decompress_and_rename_fastq(infile, outfile):
    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))  
    print(f"Decompressing {infile} to {outfile}")
    
    with gzip.open(infile, "rt") as fin, open(outfile, "wt") as fout:
        shutil.copyfileobj(fin, fout)


@collate(
    decompress_and_rename_fastq,
    regex(r".*/(.+)_\d\.fastq"),
    r"metawrap_binning.dir/\1"
)
def runMetawrapBinning(infiles, outfile):
    sample = os.path.basename(infiles[0]).rsplit("_", 1)[0]

    # Handle pooled or unpooled
    if PARAMS["is_pooled"]:
        # Pooled mode: run once using all FASTQs - skip remaining task calls
        first_sample = os.path.basename(infiles[0]).rsplit("_", 1)[0]
        if sample != first_sample:
            return

        fasta_file = f"{PARAMS['fasta_input_dir']}/*.fasta"
        fastq_args = f"{PARAMS['fastq_prepped_dir']}/*_[12].fastq"
        outdir = "metawrap_binning.dir/pooled"
    else:
        fasta_file = f"{PARAMS['fasta_input_dir']}/{sample}_corrected.spades.contigs.fasta"
        fastq_args = " ".join(infiles)  
        outdir = outfile

    threads = PARAMS["threads"]

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    out_log = os.path.join(outdir, "metawrap_binning.log")

    print("Creating output dir:", outdir)
    print("Output log path:", out_log)

    statement = ("module purge && "
                  "module load metaWRAP/1.4-20230728-foss-2023a-Python-2.7.18 && "
                  "mkdir -p %(outdir)s && "
                  "metawrap binning"
                  " -o %(outdir)s"
                  " -t %(threads)s"
                  " -a %(fasta_file)s"
                  " --metabat2 --maxbin2 --concoct"
                  " %(fastq_args)s"
                  " --run-checkm"
                  " &> %(out_log)s"
                  )

    P.run(statement,
          job_threads=threads,
          job_memory=PARAMS["job_memory"])


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    import sys
    sys.exit(main())

