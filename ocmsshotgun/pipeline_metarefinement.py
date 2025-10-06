import sys
import os
from pathlib import Path
from cgatcore import pipeline as P
from ruffus import *
import subprocess
from cgatcore import iotools as IOTools
from cgatcore import pipeline as P
import glob

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
@follows(mkdir("refined_bins.dir"))
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

@follows(run_binning_refiner)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

