'''
split_taxonomy_abundances.py
=============================

:Author: Nick Ilott
:Tags: Python

Purpose
-------

Split kraken2 abundance output into different taxonomic levels.

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python cgat_script_template.py --help

for command line help.

Command line options
--------------------

'''

import sys
import os
import cgatcore.experiment as E


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description=__doc__)

    parser.add_argument("-o", "--outdir", dest="outdir", type=str,
                        help="supply output directory")
    parser.add_argument("-p", "--prefix", dest="prefix", type=str,
                        help="supply output file prefix")

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    prefix = os.path.join(args.outdir, args.prefix)
    d_outf = open(prefix + "_domain.tsv", "w")
    k_outf = open(prefix + "_kingdom.tsv", "w")
    p_outf = open(prefix + "_phylum.tsv", "w")
    c_outf = open(prefix + "_class.tsv", "w")
    o_outf = open(prefix + "_order.tsv", "w")
    f_outf = open(prefix + "_family.tsv", "w")
    g_outf = open(prefix + "_genus.tsv", "w")
    s_outf = open(prefix + "_species.tsv", "w")
    
    for line in args.stdin.readlines():
        data = line[:-1].split("\t")
        taxon = data[0]
        counts = data[1:]
        taxonomy = taxon.split("|")
        if "d__" in taxonomy[-1]:
            d_outf.write("\t".join([taxon] + counts) + "\n")
        elif "k__" in taxonomy[-1]:
            k_outf.write("\t".join([taxon] + counts) + "\n")
        elif "p__" in taxonomy[-1]:
            p_outf.write("\t".join([taxon] + counts) + "\n")
        elif "c__" in taxonomy[-1]:
            c_outf.write("\t".join([taxon] + counts) + "\n")
        elif "o__" in taxonomy[-1]:
            o_outf.write("\t".join([taxon] + counts) + "\n")
        elif "f__" in taxonomy[-1]:
            f_outf.write("\t".join([taxon] + counts) + "\n")
        elif "g__" in taxonomy[-1]:
            g_outf.write("\t".join([taxon] + counts) + "\n")
        elif "s__" in taxonomy[-1]:
            s_outf.write("\t".join([taxon] + counts) + "\n")

    # write footer and output benchmark information.
    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
