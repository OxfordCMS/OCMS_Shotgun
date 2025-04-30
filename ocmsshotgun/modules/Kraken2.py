# Module to run Kraken2 and Bracken on paired end reads

import sys
import subprocess
import os
import re
import glob
from pathlib import Path
from ruffus import *
from cgatcore import pipeline as P

import ocmstoolkit.modules.Utility as Utility

class Kraken2(Utility.BaseTool):

    def buildStatement(self, fastn_obj):
        # Note that at the moment I only deal with paired-end
        # reads
        
        prefix = P.snip(self.outfile, ".k2.report.tsv")
    
        db = self.PARAMS["kraken2_db"]
        job_threads = self.PARAMS["kraken2_job_threads"]
        options = self.PARAMS["kraken2_options"]
        
        assert "--threads" not in options, (
            "Kraken2 multi-threading is set with job_memory and job_threads"
        )
        
        kraken_statement = ('kraken2'
                            f' --db {db}'
                            f' --output {prefix}.classified.tsv'
                            f' --report {prefix}.k2.report.tsv'
                            ' --use-names'
                            f' --threads {job_threads}')

        # paired end reads
        if fastn_obj.fastn2:
            statement_entry = (" --paired"
                               f" --gzip-compressed {fastn_obj.fastn1} {fastn_obj.fastn2}")
        # single end reads
        else:
            statement_entry = f" --gzip-compressed {fastn_obj.fastn1}"

        # build kraken statement
        kraken_statement = kraken_statement + statement_entry
        
        # add in options
        kraken_statement = kraken_statement + ' ' + options
            
        # add additional commands
        statement_entry = f'gzip {prefix}.classified.tsv'
        
        statement = kraken_statement + '; ' +  statement_entry
        
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

class Bracken(Utility.BaseTool):

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
        db = self.PARAMS["bracken_db"]
        options = self.PARAMS["bracken_options"]
        
        # run bracken at every taxonomic level
        statement = ('bracken' 
                    f' -d {db}'
                    f' -i {self.infile}'
                    f' -o bracken.dir/{prefix}.abundance.{level}.tsv'
                    f' -w bracken.dir/{prefix}.k2b.report.{level}.tsv'
                    f' -l {level_param}'
                    f' {options}')
        return statement
