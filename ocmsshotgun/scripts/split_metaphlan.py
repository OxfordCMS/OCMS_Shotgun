'''
split_metaphlan.py
==================

:Author: Sandi Yen
:Tags: Python

Purpose
-------
Split merged metaphlan output by taxonomy levels

Usage
-----
.. Example use case

Example::
    ocms_shotgun split_metaphlan

Type::
    ocms_shotgun split_metaphlan --help
    
for command line help.

Command line options
--------------------
'''

# load modules
import sys
import os
import cgatcore.experiment as E
import cgatcore.iotools as IOTools

def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    Expects merged_metaphlan.tsv as input
    """

    if argv is None:
        argv = sys.argv

    # set up command line parser
    parser = E.ArgumentParser(description = __doc__)

    parser.add_argument("-i", "--infile", dest="infile", type=str,
                        help="merged metaphlan file"),
    parser.add_argument("-o", "--outdir", dest="outdir", type=str,
                        default=os.getcwd(),
                        help="output directory; defaults to current directory")

    # add common options(-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    # create output directory if it doesn't exist
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # initialise file for each tax level
    kingdom_file = IOTools.open_file(os.path.join(args.outdir, "metaphlan_kingdom.tsv.gz"), "w")
    phylum_file = IOTools.open_file(os.path.join(args.outdir, "metaphlan_phylum.tsv.gz"), "w")
    class_file = IOTools.open_file(os.path.join(args.outdir, "metaphlan_class.tsv.gz"), "w")
    order_file = IOTools.open_file(os.path.join(args.outdir, "metaphlan_order.tsv.gz"), "w")
    family_file = IOTools.open_file(os.path.join(args.outdir, "metaphlan_family.tsv.gz"), "w")
    genus_file = IOTools.open_file(os.path.join(args.outdir, "metaphlan_genus.tsv.gz"), "w")
    species_file = IOTools.open_file(os.path.join(args.outdir, "metaphlan_species.tsv.gz"), "w")
    unknown_file = IOTools.open_file(os.path.join(args.outdir, "metaphlan_unknown.tsv.gz"), "w")

    out_list = [kingdom_file, phylum_file, class_file, order_file, family_file, genus_file, species_file, unknown_file]

    # read in header
    with IOTools.open_file(args.infile) as f:
        header = f.readline()
        header = f.readline()
    for of in out_list:
        of.write(header)
    
    # read in file by line
    with IOTools.open_file(args.infile) as f:
        # skip first line
        next(f)
        next(f)
        for line in f:
            if "s__" in line:
                species_file.write(line)
            elif "g__" in line:
                genus_file.write(line)
            elif "f__" in line:
                family_file.write(line)
            elif "o__" in line:
                order_file.write(line)
            elif "c__" in line:
                class_file.write(line)
            elif "p__" in line:
                phylum_file.write(line)
            elif "k__" in line:
                kingdom_file.write(line)
            elif line.startswith("unclassified"):
                unknown_file.write(line)
            else:
                raise ValueError("Unexpected line format in merged input: %s" % line)

    # close files
    for of in out_list:
        of.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
