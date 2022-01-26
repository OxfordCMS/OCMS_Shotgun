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
    python cgat_script_template.py

Type::
    python cgat_script_tempalte.py --help
    
for command line help.

Command line options
--------------------
'''

# load modules
import sys
import os
import cgatcore.experiment as E

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
    kingdom_file = open(os.path.join(args.outdir, "merged_metaphlan_kingdom.tsv"), "w")
    phylum_file = open(os.path.join(args.outdir, "merged_metaphlan_phylum.tsv"), "w")
    class_file = open(os.path.join(args.outdir, "merged_metaphlan_class.tsv"), "w")
    order_file = open(os.path.join(args.outdir, "merged_metaphlan_order.tsv"), "w")
    family_file = open(os.path.join(args.outdir, "merged_metaphlan_family.tsv"), "w")
    genus_file = open(os.path.join(args.outdir, "merged_metaphlan_genus.tsv"), "w")
    species_file = open(os.path.join(args.outdir, "merged_metaphlan_species.tsv"), "w")

    out_list = [kingdom_file, phylum_file, class_file, order_file, family_file, genus_file, species_file]

    # read in header
    with open(args.infile) as f:
        header = f.readline()
    for of in out_list:
        of.write(header)
    
    # read in file by line
    with open(args.infile) as f:
        # skip first line
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

    # close files
    for of in out_list:
        of.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
