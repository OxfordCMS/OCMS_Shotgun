"""
=============================
Metagenomic assembly pipeline
=============================

:Author: Jethro Johnson
:Release: $Id$
:Date: |today|
:Tags: Python

This pipeline imports cleaned, trimmed reads from one or more NGS 
experiments and creates a metagenomic assembly. Pooling of samples
is optional. 


Overview
========


Principal targets
-----------------

Optional targets
----------------



Usage
=====

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.


Input
-----

Reads
+++++


Optional inputs
+++++++++++++++

Requirements
-------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+---------+------------+------------------------------------------------+
|*Program*|*Version*   |*Purpose*                                       |
+---------+------------+------------------------------------------------+
|         |            |                                                |
+---------+------------+------------------------------------------------+
|         |            |                                                |
+---------+------------+------------------------------------------------+

Pipeline output
===============

The major output is in the database file :file:`csvdb`.

Code
====

removeDuplicateReadPairs  : bio/cdhit
assembleWithMetaSpades    : bio/spades
assembleWithMegaHit       : bio/MEGAHIT
runQuast                  : conda
mergeAssemblies           : bio/megamerge

...




"""

###############################################################################
# load modules
from ruffus import *
from ruffus.combinatorics import *

import sys
import os
import re
import glob
import shutil
import pickle
import sqlite3
import collections
import numpy as np
import pandas as pd

import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.IOTools as IOTools
import CGAT.BamTools as BamTools
import CGAT.FastaIterator as FastaIterator

import CGATPipelines.Pipeline as P
import CGATPipelines.PipelineMapping as PipelineMapping

import PipelineMetaAssembly as PMA

###############################################################################
# Pipeline configuration
P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"],
    defaults={
        'paired_end': False})

PARAMS = P.PARAMS

# determine the location of the input files (reads).
try:
    PARAMS["input"]
except NameError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]
        
###############################################################################
# Utility functions
def connect():
    '''connect to database.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])

    return dbh

def getNLines(infile):
    '''Predict the number of lines per read based on file suffix'''
    suffix = PMA.FetchData().getFormat(infile)

    if re.match('fastq', suffix):
        nb = '4'
    elif re.match('fasta', suffix):
        nb = '2'
    else:
        raise IOError('Unrecognized input file %s' % infile)

    return nb

###############################################################################
# Main Pipeline
###############################################################################
#### PREPROCESSING AND SAMPLE QC
###############################################################################
# The assumption is that the samples will have gone through minimal
# preprocessing (e.g. quality trimming) prior to being run through the
# pipeline. Howevever, the option to remove duplicates (cd-hit-dup) is included
# here. Sample read depths are also calculated. 
# In some cases it may be desirable to perform assembly on a pooled set of
# samples (e.g. for genome binning approaches such as MetaBat). The final
# preprocessing step allows for pooling some or all of the input samples based
# on file name regular expressions supplied in pipeline.ini. 
###############################################################################
## Option to remove duplicates 
###############################################################################
@follows(mkdir('input_deduped.dir'))
@follows(mkdir('read_stats.dir'))
@transform(os.path.abspath(os.path.join(DATADIR, '*.fast*.1.gz')),
           regex('.+/(.+).fast(q|a).1.gz'),
           r'input_deduped.dir/\1.fast\2.1.gz')
def removeDuplicateReadPairs(infile, outfile):
    '''Remove duplicate reads or read pairs'''

    if PARAMS['preprocess_remove_duplicates']:
        suffix = PMA.FetchData().getFormat(infile)
        
        # First count the number of reads before deduplication
        nreads = P.snip(outfile, suffix, strip_path=True) + 'prededup.nreads'
        nreads = os.path.join('read_stats.dir', nreads)
        nb = getNLines(infile)
        statement = ("zcat %(infile)s |"
                     " awk '{n+=1;} END {printf(n/%(nb)s);}'"
                     " > %(nreads)s")
        P.run()

        # Next run deduplication
        infiles = PMA.FetchData().getPairedFiles(infile)

        if len(infiles) == 1 or len(infiles) == 3:
            # Process single reads
            inf = infiles.pop()
            tmpf = P.snip(inf, '.gz')
            outf = os.path.join('input_deduped.dir', os.path.basename(tmpf))
            out_clstr = outf + '*clstr'
            
            E.info('Deduplicating file %s' % inf)
            statement = ("zcat %(inf)s > %(tmpf)s &&"
                         " cd-hit-dup"
                         "  -i %(tmpf)s"
                         "  -o %(outf)s"
                         "  %(cdhit_dup_options)s"
                         "  &> %(outfile)s.log &&"
                         " gzip %(outf)s &&"
                         " rm -f %(tmpf)s %(out_clstr)s")
            cluster_options = PARAMS['cdhit_dup_cluster_options']
            P.run()

        if len(infiles) == 2:
            # Process paired reads
            inf2 = infiles.pop()
            inf1 = infiles.pop()
            assert inf1.endswith('.1.gz')
            assert inf2.endswith('.2.gz')
            tmpf1 = P.snip(inf1, '.gz')
            tmpf2 = P.snip(inf2, '.gz')    
            outf1 = os.path.join('input_deduped.dir', os.path.basename(tmpf1))
            outf_clstr = outf1 + '*clstr'
            outf2 = os.path.join('input_deduped.dir', os.path.basename(tmpf2))

            E.info('Deduplicating files %s and %s' % (inf1, inf2))
            statement = ("zcat %(inf1)s > %(tmpf1)s &&"
                         " zcat %(inf2)s > %(tmpf2)s &&"
                         " cd-hit-dup"
                         "  -i %(tmpf1)s"
                         "  -i2 %(tmpf2)s"
                         "  -o %(outf1)s"
                         "  -o2 %(outf2)s"
                         "  %(cdhit_dup_options)s"
                         "  &>> %(outfile)s.log &&"
                         " gzip %(outf1)s &&"
                         " gzip %(outf2)s &&"
                         " rm -f %(tmpf1)s %(tmpf2)s %(outf_clstr)s")
            cluster_options = PARAMS["cdhit_dup_cluster_options"]
            P.run()
                                 
        if infiles:
            raise ValueError('Unexpected infile: %s' ' '.join(infiles))

    else:
        infiles = PMA.FetchData().getPairedFiles(infile)
        for inf in infiles:
            out_dir = os.path.dirname(outfile)
            outf = os.path.join(out_dir, os.path.basename(inf))
            p, inf = os.path.split(inf) 
            p = os.path.join('..', os.path.basename(p))
            inf = os.path.join(p, inf)
            os.symlink(inf, outf)
    
###############################################################################
## Summary statistics on input samples
###############################################################################
@follows(mkdir('read_stats.dir'))
@transform(removeDuplicateReadPairs,
           regex('.+/(.+).fast(q|a).1.gz'),
           r'read_stats.dir/\1.count')
def countReads(infile, outfile):
    '''Count the number of reads in the input files'''

    nb = getNLines(infile)
    
    statement = ("zcat %(infile)s |"
                 " awk '{n+=1;} END {printf(n/%(nb)s);}'"
                 " > %(outfile)s")
    cluster_options = "-l walltime=00:30:00,mem=2GB,nodes=1:ppn=1"
    P.run()


@active_if(PARAMS['preprocess_remove_duplicates'])
@follows(removeDuplicateReadPairs)
@merge('read_stats.dir/*.prededup.nreads',
       'read_stats.dir/prededup_sample_sizes.tsv.gz')
def mergePreDedupReadCounts(infiles, outfile):
    '''Collate read counts across samples before deduplication'''
    
    infiles = ' '.join(infiles)
    
    statement = ("cgat combine_tables"
                 " --cat SampleID"
                 " --regex-filename='.+/(.+)\.count'"
                 " --no-titles"
                 " --log=%(outfile)s.log"
                 " %(infiles)s |"
                 " gzip > %(outfile)s")
    P.run()

    
@merge(countReads, 'read_stats.dir/sample_sizes.tsv.gz')
def mergeReadCounts(infiles, outfile):
    '''Collate read counts across samples'''
    
    infiles = ' '.join(infiles)
    
    statement = ("cgat combine_tables"
                 " --cat SampleID"
                 " --regex-filename='.+/(.+)\.count'"
                 " --no-titles"
                 " --log=%(outfile)s.log"
                 " %(infiles)s |"
                 " gzip > %(outfile)s")
    P.run()

@follows(mergePreDedupReadCounts, mergeReadCounts)
def summarizeReadDepths():
    pass

###############################################################################
## Option to pool samples before assembly
###############################################################################
@follows(mkdir('input_pooled.dir'))
@collate(removeDuplicateReadPairs,
         regex('.+/' + PARAMS['preprocess_pool_input_regex']),
         os.path.join('input_pooled.dir',
                      PARAMS['preprocess_pool_output_regex']))
def poolSamples(infiles, out_fastq1):
    '''Collapse samples based on supplied regular expression'''

    out_fastq2 = P.snip(out_fastq1, '.1.gz') + '.2.gz'
    out_fastq3 = P.snip(out_fastq1, '.1.gz') + '.3.gz'

    in_fastqs1 = infiles
    in_fastqs2 = [P.snip(x, '.1.gz') + '.2.gz' for x in in_fastqs1]
    for i in in_fastqs2:
        assert os.path.exists(i)
    in_fastqs3 = [P.snip(x, '.1.gz') + '.3.gz' for x in in_fastqs1]
    for i in in_fastqs3:
        assert os.path.exists(i)

    if len(in_fastqs1) == 1:
        os.symlink(os.path.join('..', in_fastqs1[0]), out_fastq1)
        os.symlink(os.path.join('..', in_fastqs2[0]), out_fastq2)
        os.symlink(os.path.join('..', in_fastqs3[0]), out_fastq3)

    else:
        in_fastqs1 = ' '.join(in_fastqs1)
        in_fastqs2 = ' '.join(in_fastqs2)
        in_fastqs3 = ' '.join(in_fastqs3)

        statement = (" cat %(in_fastqs1)s > %(out_fastq1)s &&"
                     " cat %(in_fastqs2)s > %(out_fastq2)s &&"
                     " cat %(in_fastqs3)s > %(out_fastq3)s")
        to_cluster = True
        cluster_options = '-l walltime=12:00:00'
        P.run()

    
###############################################################################
#### METAGENOME ASSEMBLY
###############################################################################
# This section of the pipeline runs different metagenome assemblers. For legacy
# reasons the option to run read correction (BayesHammer) is included. But this
# isn't recommended when input files are likely to contain sequences from
# multiple closely-related genomes.
# Following assembly, QC stats (n50 etc) are calculated manually and via QUAST.
# Finally, if multiple assemblies are performed separately, they are merged
# so that genome binning tools such as MetaBat (which requires a single assembly
# to work with) can be run downstream. 
# A single genome assembly could be achieved by pooling all input samples prior
# to assembling (see task poolSamples above). However, in cases where either the
# canopy approach is being used on inidivual assemblies, or the pooled assembly
# exceeds available RAM, it may be necessary to initially assemble individual
# samples, or sample subsets.
#
# NB. An alternative to MeGAMerge should be found, it's slow and it relies on
# Newbler which has not been available since the decline of 454. 
###############################################################################
@active_if(PARAMS['preprocess_read_correction'])
@follows(mkdir('spades_read_correction.dir'))
#@transform(os.path.abspath(os.path.join(DATADIR, '*.fast*.1.gz')),
@transform(poolSamples,
           regex('.+/(.+).fast(q|a).1.gz'),
           r'spades_read_correction.dir/\1.fast\2.1.gz')
def runReadCorrection(infile, outfile):
    '''Run BayesHammer read correction prior to assembly'''

    cluster_options = PARAMS['spades_error_correction_run_options']
    assembler = PMA.SpadesReadCorrection()
    statement = assembler.build(infile, outfile, **PARAMS)
    # print statement
    P.run()
    
    # files output by spades are in corrected.yaml
    assembler = PMA.fetchSpadesProcessedReads()
    statement = assembler(infile, outfile)
    # print statement
    P.run()


if PARAMS['preprocess_read_correction']:
    @follows(mkdir('spades_assembly.dir'))
    @transform(runReadCorrection,
               regex('.+/(.+).fast(q|a).1.gz'),
               r'spades_assembly.dir/\1.spades.contigs.fasta')
    def assembleWithMetaSpades(infile, outfile):
        '''
        Run SPAdes --meta following read error correction
        '''
        cluster_options = PARAMS['spades_run_options']
        assembler = PMA.runMetaSpades()
        statement = assembler.build(infile, outfile, **PARAMS)

        P.run()

        
    @follows(mkdir('megahit_assembly.dir'))
    @transform(runReadCorrection,
               regex('.+/(.+).fast(q|a).1.gz'),
               r'megahit_assembly.dir/\1.megahit.contigs.fasta') 
    def assembleWithMegaHit(infile, outfile):
        '''
        Run MEGAHIT following read error correction
        '''
        cluster_options = PARAMS['megahit_run_options']
        assembler = PMA.runMegaHit()
        statement = assembler.build(infile, outfile, **PARAMS)
        
        P.run()


    @follows(mkdir('minia_assembly.dir'))
    @originate('minia_assembly.dir/stool-all-R000.minia.contigs.fasta.gz')
    def assembleWithMinia(outfile):
        '''Hack, insert minia assembly into the pipeline'''
        pass
    
else:
    @follows(mkdir('spades_assembly.dir'))
    # @transform(os.path.abspath(os.path.join(DATADIR, '*.fast*.1.gz')),
    @transform(poolSamples,
               regex('.+/(.+).fast(q|a).1.gz'),
               r'spades_assembly.dir/\1.spades.contigs.fasta')
    def assembleWithMetaSpades(infile, outfile):
        '''
        Run SPAdes --meta without read error correction
        '''
        cluster_options = PARAMS['spades_run_options']
        assembler = PMA.runMetaSpades()
        statement = assembler.build(infile, outfile, **PARAMS)

        P.run()

    
    @follows(mkdir('megahit_assembly.dir'))
    # @transform(os.path.abspath(os.path.join(DATADIR, '*.fast*.1.gz')),
    @transform(poolSamples,
               regex('.+/(.+).fast(q|a).1.gz'),
               r'megahit_assembly.dir/\1.megahit.contigs.fasta') 
    def assembleWithMegaHit(infile, outfile):
        '''
        Run MEGAHIT without read error correction
        '''
        cluster_options = PARAMS['megahit_run_options']
        assembler = PMA.runMegaHit()
        statement = assembler.build(infile, outfile, **PARAMS)

        P.run()


    @follows(mkdir('minia_assembly.dir'))
    @originate('minia_assembly.dir/stool-all-R000.minia.contigs.fasta')
    def assembleWithMinia(outfile):
        '''Hack, insert minia assembly into the pipeline'''
        pass

###############################################################################
# Run Assembly
ASSEMBLY_TARGETS = []
assembly_targets = {'metaspades': (assembleWithMetaSpades,),
                    'megahit': (assembleWithMegaHit,),
                    'minia': (assembleWithMinia,)
                    }

for tool in PARAMS['assemblers'].split(','):
    ASSEMBLY_TARGETS.extend(assembly_targets[tool])

@follows(*ASSEMBLY_TARGETS)
def assembleMetaGenome():
    '''Perform assembly using tools specified in config file'''
    pass

###############################################################################
# Calculate Assembly Statistics
@transform(ASSEMBLY_TARGETS,
           regex('(.+)/(.+).contigs.fasta'),
           r'\1/\2.contigs.tsv.gz')
def getContigStats(contig_file, outfile):
    '''Run fasta2table to get contig/scaffold length/ngaps'''

    scaffold_file = P.snip(contig_file, '.contigs.fasta') + '.scaffolds.fasta'
    scaffold_out = P.snip(outfile, 'contigs.tsv.gz') + 'scaffolds.tsv.gz'
    
    statement = ("cat %(contig_file)s |"
                 " cgat fasta2table"
                 "  --section=length,gaps"
                 "  --log=%(outfile)s.log |"
                 " gzip > %(outfile)s")

    statement2 = ("cat %(scaffold_file)s |"
                 " cgat fasta2table"
                 "  --section=length,gaps"
                 "  --log=%(outfile)s.log |"
                 " gzip > %(scaffold_out)s")
    if os.path.exists(scaffold_file):
        statement = statement + ' && ' + statement2
        
    P.run()


@collate(getContigStats,
         regex('(.+)/(.+)\.(.+)\.contigs.tsv.gz'),
         r'\1/\3_contig_stats.tsv.gz')
def collateContigStats(infiles, outfile):
    '''Collate contig stats for each assembler'''

    scaffold_files = ' '.join([P.snip(x, '.contigs.tsv.gz') + \
                              '.scaffolds.tsv.gz' for x in infiles])
    scaffold_out = P.snip(outfile, 'contig_stats.tsv.gz') + 'scaffold_stats.tsv.gz'
    contig_files = ' '.join(infiles)
    
    statement = ("cgat combine_tables"
                 " --cat SampleID"
                 " --regex-filename='.+/(.+)\..+\.contigs.tsv.gz'"
                 " --log=%(outfile)s.log"
                 " %(contig_files)s |"
                 " gzip > %(outfile)s")
    
    statement2 = ("cgat combine_tables"
                  "  --cat SampleID"
                  "  --regex-filename='.+/(.+)\..+\.scaffolds.tsv.gz'"
                  "  --log=%(outfile)s.log"
                  "  %(scaffold_files)s |"
                  "  gzip > %(scaffold_out)s")
    if os.path.exists(scaffold_files.split()[0]):
        statement = statement + ' && ' + statement2

    P.run()


@merge(collateContigStats, 'assembly_contig_stats.tsv.gz')
def collateContigStatsAcrossAssemblers(infiles, outfile):
    '''Combine contig stats for different assemblers'''

    # scaffolds are not reported for all assemblers
    scaffold_files = [P.snip(x, '_contig_stats.tsv.gz') + \
                      '_scaffold_stats.tsv.gz' for x in infiles]
    scaffold_files = [x for x in scaffold_files if os.path.exists(x)] 
    scaffold_out = P.snip(outfile, 'contig_stats.tsv.gz') + 'scaffold_stats.tsv.gz'

    contig_files = ' '.join(infiles)
    scaffold_files = ' '.join(scaffold_files)
    
    # HACK: Assumes combine_tables is in the pipeline directory
    call_dir = os.path.dirname(os.path.realpath(__file__))

    statement = ("cgat combine_tables"
                 " --cat assembler"
                 " --regex-filename='.+/(.+)_contig_stats.tsv.gz'"
                 " --log=%(outfile)s.log"
                 " %(contig_files)s |"
                 " gzip > %(outfile)s")
 
    statement2 = ("cgat combine_tables"
                  "  --cat assembler"
                  "  --regex-filename='.+/(.+)_scaffold_stats.tsv.gz'"
                  "  --log=%(outfile)s.log"
                  "  %(scaffold_files)s |"
                  "  gzip > %(scaffold_out)s")
    if scaffold_files:
        statement = statement + ' && ' + statement2

    P.run()


# @transform(collateContigStatsAcrossAssemblers, suffix('.tsv.gz'), '.load')
# def loadContigStats(contig_file, outfile):

#     assert os.path.exists(contig_file)
#     contig_table = P.snip(contig_file, '.tsv.gz', strip_path=True)
#     df_contig = pd.read_table(contig_file, index_col=0)
#     df_contig.to_sql(name=contig_table, con=connect(), if_exists='replace')

#     open(outfile, 'w').close()

#     scaffold_file = re.sub('_contig_', '_scaffold_', contig_file)
#     if os.path.exists(scaffold_file):
#         scaffold_table = P.snip(scaffold_file, '.tsv.gz', strip_path=True)
#         df_scaff = pd.read_table(scaffold_file, index_col=0)
#         df_scaff.to_sql(name=scaffold_table, con=connect(), if_exists='replace')
        
#         open(re.sub('_contig_', '_scaffold_', outfile), 'w').close()


#@split(loadContigStats, 'assembly_*_summary_stats.tsv')
@split(collateContigStatsAcrossAssemblers, 'assembly_*_summary_stats.tsv') 
def calcAssemblyStats(infile, outfiles):
    '''Calculate summary statistics for each genome assembly, median 
    contig size, n50 etc.
    NEEDS TO BE EXTRACTED.... SLOW!
    '''

#    contig_table = P.snip(infile, '.load', strip_path=True)
    contig_table = infile
    contig_out = 'assembly_contig_summary_stats.tsv'
    tables = [[contig_table, contig_out],]
    
    if os.path.exists(re.sub('_contig_', '_scaffold_', infile)):
        scaff_table = re.sub('_contig_', '_scaffold_', contig_table)
        scaff_out = 'assembly_scaffold_summary_stats.tsv'
        tables.append([scaff_table, scaff_out])
        
    for inf, outf in tables:
        E.info('LOADING  TABLE %s INTO PANDAS' % inf)
#        statement = 'SELECT * FROM {}'.format(inf)
#        df = pd.read_sql(statement,
#                         con=connect(),
#                         index_col=['SampleID', 'assembler', 'id'])
        df = pd.read_table(inf,
                           compression='gzip',
                           index_col=['SampleID', 'assembler', 'id'])
        E.info('LOADED TABLE %s INTO PANDAS' % inf)
        
        # calculate assembly statistics
        sample_assembler = df.groupby(level=['SampleID', 'assembler'])
        medians = sample_assembler.median()[['length', 'ngaps']]
        medians.columns = ['median_length', 'median_ngaps']
        counts = sample_assembler.count()[['length']]
        counts.columns = ['n_contigs']
        maxs = sample_assembler.max()[['length', 'ngaps']]
        maxs.columns = ['max_length', 'max_ngaps']
        totals = sample_assembler.sum()[['length']]
        totals.columns = ['total_length']
        
        lengths = df[['length']].groupby(level=['SampleID', 'assembler'])
        n50 = lambda x: (np.median([item for sublist in [[n,]* n for n in x]\
                                    for item in sublist]))
        n90 = lambda x: (np.percentile([item for sublist in [[n,]* n for n in x]\
                                        for item in sublist], 10))
        n50s = lengths.aggregate(n50)
        n50s.columns = ['n50']
        n90s = lengths.aggregate(n90)
        n90s.columns = ['n90']
        
        df_out = medians.join([counts, maxs, totals, n50s, n90s])
        df_out.to_csv(outf, sep='\t')
        

@transform(calcAssemblyStats, suffix('.tsv'), '.load')
def loadAssemblyStats(infile, outfile):

    table_name = P.snip(infile, '.tsv', strip_path=True)
    df = pd.read_table(infile, index_col=0)
    df.to_sql(name=table_name, con=connect(), if_exists='replace')

    open(outfile, 'w').close()


###############################################################################
# Calculate QUAST assembly statistics
@transform(ASSEMBLY_TARGETS,
            regex('(.+)/(.+).fasta'),
            r'\1/\2_quast.dir/report.tsv')
def runQuast(infile, outfile):
    '''Run QUAST with a set of reference genomes to compare against'''

    out_dir = os.path.dirname(outfile)
    out_log = P.snip(out_dir, '.dir') + '.log'

    statement = ("metaquast.py %(infile)s"
                 " --output-dir %(out_dir)s"
                 " -r %(quast_reference)s"
                 " %(quast_options)s"
                 " &> %(out_log)s")
    cluster_options = PARAMS['quast_run_options']
    P.run()
   
###############################################################################
# MergeAssemblies
@collate(ASSEMBLY_TARGETS,
         regex('(.+)_assembly.dir/(.+)\.(.+)\.contigs.fasta'),
         r'\1_assembly_merged.dir/merged_\1_contigs.fasta')
def mergeAssemblies(infiles, outfile):
    '''If necessary, pool multiple assemblies into a single multi-fasta
    using MeGAMerge.
    '''

    outd = os.path.dirname(outfile)
    if not os.path.exists(outd):
        os.mkdir(outd)
    
    if len(infiles) == 1:
        os.symlink(os.path.relpath(infiles[0], os.path.dirname(outfile)),
                   outfile)

    else:
        infiles = ' '.join(infiles)
        out_dir, outf = os.path.split(outfile)
        statement = ("MeGAMerge-1.1.pl"
                     " -single_genome=0" # This pipeline is for metagenomes
                     " -force"
                     " -d"
                     " %(megamerge_options)s"
                     " -o=%(outf)s"
                     " %(out_dir)s"
                     " %(infiles)s"
                     " &> %(outfile)s.log")
        cluster_options = PARAMS['megamerge_run_options']
        P.run()
    
    
###############################################################################
# Pipeline Targets
###############################################################################
@follows(assembleMetaGenome)
def Assemble():
    pass


@follows(loadAssemblyStats)
def full():
    pass


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))



