# Module to run MetaPHlAn4 on paired end reads

import sys
import subprocess
import os
import re
import glob
from pathlib import Path
from ruffus import *
from cgatcore import pipeline as P

import ocmsshotgun.modules.Utility as utility

class metaphlan(utility.matchReference):

    def buildStatement(self):
        p1 = self.fastq1
        prefix = P.snip(self.outfile, ".metaphlan_profile.txt")
    
        db_path = os.path.join(
            self.PARAMS["metaphlan_db"], 
            self.PARAMS["metaphlan_database_version"])
            
        job_threads = self.PARAMS["metaphlan_threads"]
        options = self.PARAMS["metaphlan_options"]
        
        metaphlan_statement = ('metaphlan'
                             ' --bowtie2db %(db_path)s'
                             ' --nproc %(job_threads)s'
                             ' --input_type fastq'
                             ' --index %(db_version)s'
                             ' --bowtie2out %(prefix)s.bowtie2.bz2'
                             ' --tax_lev %(tax_lev)s'
                             ' -o %(outfile)s' % {
                                 'db_path': db_path,
                                 'job_threads': job_threads,
                                 'db_version': self.PARAMS["metaphlan_database_version"],
                                 'prefix': prefix,
                                 'tax_lev': self.PARAMS["metaphlan_tax_lev"],
                                 'outfile': self.outfile
                             })

        if self.fastq2:
            p2 = p1.replace(".fastq.1.gz", ".fastq.2.gz")
            statement_entry = " %s,%s" % (p1, p2)
        else:
            statement_entry = " %s" % p1

        statement = metaphlan_statement + statement_entry

        if options:
            statement += " " + options
        
        return statement

def extract_taxonomy_level(infile, outfile, level):
    statement = '''grep -E "%(level)s|clade" %(infile)s | 
                  grep "%(level)s" | 
                  sed 's/^.*%(level)s//g' > %(outfile)s
    ''' % {'level': level, 'infile': infile, 'outfile': outfile}
    
    return statement

def merge_tables(infiles, outfile):
    infiles_str = " ".join(infiles)
    statement = '''merge_metaphlan_tables.py 
                  %(infiles)s
                  > %(outfile)s
    ''' % {'infiles': infiles_str, 'outfile': outfile}
    
    return statement