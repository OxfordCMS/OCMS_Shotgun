# Module to run Kraken2 and Bracken on paired end reads

import sys
import subprocess
import os
import re
import glob
from pathlib import Path
from ruffus import *
from cgatcore import pipeline as P

class matchReference(object):
    """
    Base class for generating run statements to match mWGS reads to 
    reference sequences. Intended to work with single, paired, or
    paired + singleton fastn files. 

    Some options are  assumed to be passed via kwargs, as this and 
    inherited classes are writtento work with a PARAMS dict 
    generated from a pipeline.yml config file.

    ** Options:
    fn_suffix - option to pass something other than .fastq.1.gz

    


    """

    def __init__(self, fastn1, outfile, **PARAMS):
        self.outdir = os.path.dirname(outfile)

        self.fastn1 = fastn1

        # Assume that files are fastq and end .fastq.1.gz
        if PARAMS.get('fn_suffix', None):
            self.fn_suffix = PARAMS.get('fn_suffix')
        else:
            assert self.fastn1.endswith('.fastq.1.gz'), \
                'Please supply fastn suffix using keyword fn_suffix='
            self.fn_suffix = '.fastq.1.gz'

        # Find mate pair file
        n = self.fn_suffix.rfind('1')
        fn2_suffix =  self.fn_suffix[:n] + '2' + self.fn_suffix[n+1:]
        fastn2 = P.snip(self.fastn1, self.fn_suffix) + fn2_suffix

        if os.path.exists(fastn2):
            self.fastn2 = fastn2
        else:
            self.fastn2 = None

        # Find singleton file
        fn3_suffix =  self.fn_suffix[:n] + '3' + self.fn_suffix[n+1:]
        fastn3 = P.snip(self.fastn1, self.fn_suffix) + fn3_suffix

        if os.path.exists(fastn3):
            assert self.fastn2, "Can't have singletons without mate pairs"
            self.fastn3 = fastn3
        else:
            self.fastn3 = None

class kraken2(matchReference):
    def __init__(self, fastn1, outfile, **PARAMS):
        super().__init__(fastn1, outfile, **PARAMS)

    def run(self, infile, outfile, *args, **PARAMS):
        '''classify reads with kraken2
        '''
        
        # Note that at the moment I only deal with paired-end
        # reads
        p1 = infile
        
        prefix = P.snip(outfile, ".k2.report.tsv")
    
        db = PARAMS.get("kraken2_db")
        job_threads = PARAMS.get("kraken2_job_threads")
        job_memory = PARAMS.get("kraken2_job_mem")
        options = PARAMS.get("kraken2_options")
    
        kraken_statement = ('kraken2',
                            '--db %(db)s',
                            '--output %(prefix)s.classified.tsv',
                            '--report %(prefix)s.k2.report.tsv',
                            '--use-names',
                            '--threads %(job_threads)s')

        # paired end reads
        if self.fastn2:
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
        
        P.run(statement)


class utility():

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
        check all backen at all taxonomic levels has been run
        '''
        expected_files = [x for x in list(utility.generate_parallel_params())]
        expected_files = [x[1] for x in expected_files]
        return expected_files

    def check_bracken_levels(expected_files, outfile):
    
        missing_files = [f for f in expected_files if not os.path.exists(f)]
    
        if missing_files:
            raise Exception("Abundance file not found:\n%s" % '\n'.join(missing_files))
        else:
            open(outfile, 'a').close()

class bracken():

    def run(infile, outfile, **PARAMS):
        '''
        convert read classifications into abundance with Bracken
        '''
        # get sample name
        pattern = r"(.*)\.abundance\.(species|genus|family|class|order|phylum|domain)\.tsv$"
        match = re.search(pattern, os.path.basename(outfile))
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
        db = PARAMS.get("bracken_db")
        read_len = PARAMS.get("bracken_read_len")
        job_threads = PARAMS.get("bracken_job_threads")
        job_memory = PARAMS.get("bracken_job_mem")
        options = PARAMS.get("bracken_options")

    
        # run bracken at every taxonomic level
        statement = '''bracken 
                    -d %(db)s
                    -i %(infile)s
                    -o bracken.dir/%(prefix)s.abundance.%(level)s.tsv
                    -w bracken.dir/%(prefix)s.k2b.report.%(level)s.tsv
                    -l %(level_param)s
                    %(options)s
                    '''
        P.run(statement)
