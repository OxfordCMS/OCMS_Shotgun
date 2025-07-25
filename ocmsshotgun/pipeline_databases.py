"""===========================
pipeline_databases.py
===========================

Overview
========

This pipeline creates relevant shotgun metagenomics databases for use with the various pipelines. It includes downloading and indexing.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use cgat pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.yml` file.

Default configuration files can be generated by executing:

   ocms_shotgun databases config

Input files
-----------

The only required input file is the pipeline.yml

Dependencies
------------

Dependencies will depend a bit on which tasks you want to run.


Pipeline output
===============


Glossary
========

.. glossary::


Code
====

"""
import sys
import os
import re
import glob
import yaml
from pathlib import Path
from ruffus import *
import ocmsshotgun.modules.Databases as DB
from cgatcore import pipeline as P
from cgatcore import iotools as IOTools

PARAMS = P.get_parameters(["pipeline.yml"])

########################################################
########################################################
# get the general information on versions etc
########################################################
########################################################


try:
    IOTools.open_file("pipeline.yml")
except FileNotFoundError as e:
    # allow config to make yml
    pass

gcc_version = PARAMS["gcc"]
python_version = PARAMS["python"]
human_build = PARAMS["genomes_human"]
mouse_build = PARAMS["genomes_mouse"]
srprism_version = PARAMS["srprism_version"]
bmtool_version = PARAMS["bmtool_version"]
sortmerna_version = PARAMS["sortmerna_version"]
kraken2_version = PARAMS["kraken2_version"]
metaphlan_version = PARAMS["metaphlan_version"]
hisat2_version = PARAMS["hisat2_version"]
minimap2_version = PARAMS["minimap2_version"]
trimmomatic_version = PARAMS["trimmomatic_version"]

########################################################
########################################################
########################################################
# Mammalian datasets: plain fasta files 
########################################################
########################################################
########################################################

@follows(mkdir(f"genomes/human/{human_build}"))
@follows(mkdir(f"genomes/mouse/{mouse_build}"))
@split(None,
       [f"genomes/mouse/{mouse_build}/{mouse_build}.fa.gz",
        f"genomes/human/{human_build}/{human_build}.fa.gz"])
def getMammalianGenomes(infile, outfiles):
    '''
    wget mammalian genomes: at the moment we deal
    with just human and mouse
    '''
    # has to be run locally on node that
    # has internet access
    without_cluster=True
    statement = DB.get_mammalian_genomes() % globals()
    P.run(statement)

########################################################
########################################################
########################################################
# Mammalian datasets: plain fasta files 
########################################################
########################################################
########################################################

@follows(mkdir(f"HISAT2/human/{human_build}/GCC-{gcc_version}/hisat2-{hisat2_version}"))
@follows(mkdir(f"HISAT2/mouse/{mouse_build}/GCC-{gcc_version}/hisat2-{hisat2_version}"))
@transform(getMammalianGenomes,
           regex("genomes/(\S+)/(\S+)/(\S+).fa.gz"),
           r"HISAT2/\1/\2/GCC-%(gcc_version)s/hisat2-%(hisat2_version)s/\2.ht2" % globals())
def buildHisat2(infile, outfile):
    '''
    build hisat2 indexes
    '''
    statement = DB.HISAT2(PARAMS, "hisat2").build_statement(infile, outfile) + "; touch %(outfile)s"
    P.run(statement)

    
########################################################
########################################################
########################################################
# SRPRISM
########################################################
########################################################
########################################################

@follows(mkdir(f"SRPRISM/human/{human_build}/GCC-{gcc_version}/srprism-{srprism_version}"))
@follows(mkdir(f"SRPRISM/mouse/{mouse_build}/GCC-{gcc_version}/srprism-{srprism_version}"))
@transform(getMammalianGenomes,
           regex("genomes/(\S+)/(\S+)/(\S+).fa.gz"),
           r"SRPRISM/\1/\2/GCC-%(gcc_version)s/srprism-%(srprism_version)s/\3_index.sentinel" % globals())
def buildSrprismIndex(infile, outfile):
    '''
    build srprism indexes
    '''
    statement = DB.SRPRISMdb(PARAMS, "srprism").build_statement(infile, outfile) + "; touch %(outfile)s"
    P.run(statement)

########################################################
########################################################
########################################################
# bmtagger
########################################################
########################################################
########################################################

@follows(mkdir(f"bmtool/human/{human_build}/GCC-{gcc_version}/bmtool-{bmtool_version}"))
@follows(mkdir(f"bmtool/mouse/{mouse_build}/GCC-{gcc_version}/bmtool-{bmtool_version}"))
@transform(getMammalianGenomes,
           regex("genomes/(\S+)/(\S+)/(\S+).fa.gz"),
           r"bmtool/\1/\2/GCC-%(gcc_version)s/bmtool-%(bmtool_version)s/\3.bitmask" % globals())
def buildBitmask(infile, outfile):
    '''
    build bmtool indexes
    '''
    statement = DB.BMTOOLdb(PARAMS, "bmtool").build_statement(infile, outfile)
    P.run(statement)

########################################################
########################################################
########################################################
# sortmerrna
########################################################
########################################################
########################################################

@follows(mkdir(f"sortmerna//GCC-{gcc_version}/sortmerna-{sortmerna_version}/rrna"))
@split(None, f"sortmerna/GCC-{gcc_version}/sortmerna-{sortmerna_version}/rrna/*.fasta")
def getSortmerna(infile, outfiles):
    '''
    get pre-made sortmerna fasta datasets
    '''
    # dummy outfiles
    outfiles = ["sortmerna/GCC-%(gcc_version)s/sortmerna-%(sortmerna_version)s/rrna/dummy.fasta" % globals()] 
    statement = DB.SORTMERNAdb(PARAMS, "sortmerna").get_plain_fasta(outfiles)

    without_cluster = True
    P.run(statement)

########################################################
########################################################
########################################################

@follows(mkdir(f"sortmerna/GCC-{gcc_version}/sortmerna-{sortmerna_version}/index"))
@transform(getSortmerna,
           regex("sortmerna/(\S+)/(\S+)/(\S+)/(\S+).fasta"),
           r"sortmerna/GCC-%(gcc_version)s/sortmerna-%(sortmerna_version)s/index/\4" % globals())
def buildSortmerna(infile, outfile):
    '''
    build the sortmerna indexes
    '''
    statement = DB.SORTMERNAdb(PARAMS, "sortmerna").build_statement(infile, outfile)
    P.run(statement)

########################################################
########################################################
########################################################
# Preprocess database files 
########################################################
########################################################
########################################################

@follows(buildHisat2,
         buildSrprismIndex,
         buildBitmask,
         buildSortmerna)
def buildPreprocessDatabases():
    pass


########################################################
########################################################
########################################################

@follows(mkdir(f"kraken2/kraken2-{kraken2_version}"))
@split(None, "kraken2/kraken2-%(kraken2_version)s/*/*.input" % globals())
def getDummyInputs(infile, outfiles):
    '''
    get some dummy inputs for downloading the kraken2 databases
    '''
    statement = DB.KRAKEN2db(PARAMS, "kraken2").make_dummy_files()
    P.run(statement)

########################################################
########################################################

@transform(getDummyInputs, regex("(\S+).input"), r"\1_index.sentinel")
def getKraken2Index(infile, outfile):
    '''
    get the kraken2 indexes
    '''
    without_cluster = True
    statement = DB.KRAKEN2db(PARAMS, "kraken2").build_statement(infile, outfile)
    P.run(statement)


########################################################
########################################################
########################################################
# Kraken2 databases
########################################################
########################################################
########################################################

@follows(getKraken2Index)
def buildKraken2Databases():
    pass



########################################################
########################################################
########################################################

metaphlan_version = metaphlan_version.split(".")[0]
@follows(mkdir(f"metaphlan/python-{python_version}/metaphlan-{metaphlan_version}"))
@split(None, "metaphlan/python-%(python_version)s/metaphlan-%(metaphlan_version)s/*.bt2*" % globals())
def getMetaphlanIndex(infile, outfiles):
    '''
    get some dummy inputs for downloading the kraken2 databases
    '''
    without_cluster = True
    statement = DB.METAPHLANdb(PARAMS, "metaphlan").build_statement(infile, outfiles)
    P.run(statement)

########################################################
########################################################
########################################################
# metaphlan databases
########################################################
########################################################
########################################################

@follows(getMetaphlanIndex)
def buildMetaphlanDatabases():
    pass

########################################################
########################################################
########################################################
# minimap2 databases
########################################################
########################################################
########################################################

@follows(mkdir(f"minimap2/human/{human_build}/GCC-{gcc_version}/minimap2-{minimap2_version}"))
@follows(mkdir(f"minimap2/mouse/{mouse_build}/GCC-{gcc_version}/minimap2-{minimap2_version}"))
@transform(getMammalianGenomes,
           regex("genomes/(\S+)/(\S+)/(\S+).fa.gz"),
           fr"minimap2/\1/\2/GCC-{gcc_version}/minimap2-{minimap2_version}/\3.mmi")
def getMinimap2Index(infile, outfile):
    '''
    get minimap2 indexes
    '''
    tool = DB.MINIMAP2db(PARAMS, 'minimap2')
    statement = tool.build_statement(infile, outfile)

    P.run(statement)

@follows(getMinimap2Index)
def buildMinimap2Databases():
    pass

########################################################
########################################################
# Trimmomatic adpater fasta files
########################################################
########################################################
########################################################

@follows(mkdir(f"Trimmomatic/trimmomatic-{trimmomatic_version}"))
@split(None, "Trimmomatic/trimmomatic-{trimmomatic_version}/*.fasta")
def getTrimmomaticSequences(infile, outfiles):
    '''
    just downloads adapter fasta files for
    generic use
    '''
    without_cluster = True
    tool = DB.Trimmomatic(PARAMS, "trimmomatic")
    statement = tool.build_statement(infile, outfiles)

    P.run(statement)


# ---------------------------------------------------
# Generic pipeline tasks
@follows(buildPreprocessDatabases,
         buildKraken2Databases, 
         buildMinimap2Databases,
         getTrimmomaticSequences)
def full():
    pass

def main(argv=None): 
    if argv is None: 
        argv = sys.argv 
    P.main(argv)


if __name__ == "__main__": 
    sys.exit(P.main(sys.argv))
