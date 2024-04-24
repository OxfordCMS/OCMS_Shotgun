# Module to run Kraken2 and Bracken on paired end reads

import sys
import subprocess
import os
import re
import glob
from pathlib import Path
from ruffus import *
from cgatcore import pipeline as P

import ocmsshotgun.modules.Utility as utility

class kraken2(utility.matchReference):
    def __init__(self, fastq1, outfile, **PARAMS):
        super().__init__(fastq1, outfile, **PARAMS)
    
    def builtStatement(self):
        # Note that at the moment I only deal with paired-end
        # reads
        p1 = infile
        
        prefix = P.snip(self.outfile, ".k2.report.tsv")
    
        db = self.PARAMS.get("kraken2_db")
        options = self.PARAMS.get("kraken2_options")
    
        kraken_statement = ('kraken2',
                            '--db %(db)s',
                            '--output %(prefix)s.classified.tsv',
                            '--report %(prefix)s.k2.report.tsv',
                            '--use-names',
                            '--threads %(job_threads)s')

        # paired end reads
        if self.fastq2:
            p2 = p1.replace(".fastq.1.gz", ".fastq.2.gz")
            statement_entry = ["--paired",
                               "--gzip-compressed %s %s" % (p1, p2)]
        # single end reads
        else:
            statement_entry = ["--gzip-compressed %s" % p1]

        # build kraken statement
        kraken_statement = ' '.join(list(kraken_statement) + statement_entry)
        
        if options != '':
            statement_entry = (*kraken_statement, '%(options)s')

        # add additional commands
        statement_entry = ['gzip %(prefix)s.classified.tsv']
        
        statement = ';'.join([kraken_statement] +  statement_entry)
        
        return statement

def generate_parallel_params():
    '''
    use output of runKraken2 to generate list  of bracken inputs and 
    outputs for every taxononmic level    
    '''
    kraken_outfiles = glob.glob("kraken2.dir/*.k2.report.tsv")
    levels = ['species','genus','family','class','order','phylum','domain']

    for kraken_outfile in kraken_outfiles:
        # generate output files for each taxonomic level
        sample_name = os.path.basename(kraken_outfile)
        bracken_outfiles = [sample_name.replace(".k2.report.tsv",".abundance.%s.tsv" % x) for x in levels]
        bracken_outfiles = [os.path.join("bracken.dir", x) for x in bracken_outfiles]

        for bracken_outfile in bracken_outfiles:
            yield kraken_outfile, bracken_outfile

def files_to_check():
    '''
    check all bracken at all taxonomic levels has been run
    '''
    expected_files = [x for x in list(generate_parallel_params())]
    expected_files = [x[1] for x in expected_files]
    return expected_files

def check_bracken_levels(expected_files, outfile):

    missing_files = [f for f in expected_files if not os.path.exists(f)]

    if missing_files:
        raise Exception("Abundance file not found:\n%s" % '\n'.join(missing_files))
    else:
        open(outfile, 'a').close()

class bracken(utility.matchReference):
    def __init__(self, fastq1, outfile, **PARAMS):
        super().__init__(fastq1, outfile, **PARAMS)

    def buildStatement(self):
        '''
        convert read classifications into abundance with Bracken
        '''
        # get sample name
        pattern = r"(.*)\.abundance\.(species|genus|family|class|order|phylum|domain)\.tsv$"
        match = re.search(pattern, os.path.basename(self.outfile))
        try:
            prefix = match.group(1)
        except:
            raise Exception("Sample name not found with regex of outfile")
    
        # taxonomic level
        level = match.group(2)
        level_param = level[0].capitalize()
        
        # check if sentinel file exists
        #if os.path.exists(os.path.join("bracken.dir", sentinel)):
        #    return
    
        # bracken parameters
        db = self.PARAMS.get("bracken_db")
        read_len = self.PARAMS.get("bracken_read_len")
        options = self.PARAMS.get("bracken_options")

    
        # run bracken at every taxonomic level
        statement = '''bracken 
                    -d %(db)s
                    -i %(infile)s
                    -o bracken.dir/%(prefix)s.abundance.%(level)s.tsv
                    -w bracken.dir/%(prefix)s.k2b.report.%(level)s.tsv
                    -l %(level_param)s
                    %(options)s
                    '''
        return statement
