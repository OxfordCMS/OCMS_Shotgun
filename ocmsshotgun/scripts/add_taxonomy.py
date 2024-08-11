'''
add_taxonomy.py
==================

:Author: Nick Ilott
:Tags: Python

Purpose
-------
Add taxonomy information as the feature ids to the output from bracken merged tables. The
idea is that this makes it easier to track clades etc without multiple files.

Usage
-----
.. Example use case

Example::
    ocms_shotgun add_taxonomy -c counts_table -t taxonomy --log=add_taxonomy.log

Type::
    ocms_shotgun add_taxonomy --help
    
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
    """

    if argv is None:
        argv = sys.argv

    # set up command line parser
    parser = E.ArgumentParser(description = __doc__)

    parser.add_argument("-c", "--counts", dest="counts", type=str,
                        help="merged bracken counts table"),
    parser.add_argument("-t", "--taxonomy", dest="taxonomy", type=str,
                        help="taxonomy information from bracken2mpataxonomy.py")

    # add common options(-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    # read in taxonomy: name + "-" + taxid as index
    taxonomy = {}
    tax = open(args.taxonomy)
    for line in tax.readlines():
        data = line.strip("\n").split("\t")
        name = data[0].split(" ")
        taxid = name[-1]
        name = " ".join([x for x in name if x not in ["taxid", taxid]])
        ind = name + "-" + taxid
        taxonomic_anno = data[1].split(";")
        taxonomic_anno = ";".join([x for x in taxonomic_anno if "unassigned" not in x])
        taxonomy[ind] = taxonomic_anno

    # read in counts: 1st column as index
    # and replace the 1st column
    counts = open(args.counts)
    header = counts.readline().strip("\n").split("\t")
    header = ["taxon"] + header[1:]
    sys.stdout.write("\t".join(header) + "\n")
    for line in counts.readlines():
        data = line.strip("\n").split("\t")
        new_name = taxonomy[data[0]]
        sys.stdout.write(new_name + "\t" + "\t".join(data[1:]) + "\n")

if __name__ == "__main__":
    sys.exit(main(sys.argv))
