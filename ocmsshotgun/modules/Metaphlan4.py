# Module to run MetaPHlAn4 on paired end reads

import sys
import subprocess
import os
import re
import globs
from pathlib import Path
from ruffus import *
from cgatcore import pipeline as P

import ocmsshotgun.modules.Utility as utility

class metaphlan(utility.matchReference):

    def buildStatement(self):
        # Handle paired-end reads
        p1 = self.fastq1
        prefix = P.snip(self.outfile, ".metaphlan_profile.txt")
    
        db_path = os.path.join(
            self.PARAMS["metaphlan_db"], 
            self.PARAMS["metaphlan_database_version"])
            
        job_threads = self.PARAMS["metaphlan_threads"]
        
        metaphlan_statement = ('metaphlan'
                             ' --bowtie2db %(db_path)s'
                             ' --nproc %(job_threads)s'
                             ' --input_type fastq'
                             ' --index %(db_version)s'
                             ' --bowtie2out %(prefix)s.bowtie2.bz2'
                             ' --min_cu_len %(min_cu_len)s'
                             ' --min_alignment_len %(min_align_len)s'
                             ' --tax_lev %(tax_lev)s'
                             ' --add_viruses %(add_vir)s'
                             ' --stat_q %(stat_q)s'
                             ' --perc_nonzero %(perc_nonzero)s'
                             ' --avoid_disqm %(avoid_disqm)s'
                             ' --ignore_viruses %(ignore_vir)s'
                             ' -o %(outfile)s' % {
                                 'db_path': db_path,
                                 'job_threads': job_threads,
                                 'db_version': self.PARAMS["metaphlan_database_version"],
                                 'prefix': prefix,
                                 'min_cu_len': self.PARAMS["metaphlan_min_cu_len"],
                                 'min_align_len': self.PARAMS["metaphlan_min_alignment_len"],
                                 'tax_lev': self.PARAMS["metaphlan_tax_lev"],
                                 'add_vir': self.PARAMS["metaphlan_add_viruses"],
                                 'stat_q': self.PARAMS["metaphlan_stat_q"],
                                 'perc_nonzero': self.PARAMS["metaphlan_perc_nonzero"],
                                 'avoid_disqm': self.PARAMS["metaphlan_avoid_disqm"],
                                 'ignore_vir': self.PARAMS["metaphlan_ignore_viruses"],
                                 'outfile': self.outfile
                             })

        # paired end reads
        if self.fastq2:
            p2 = p1.replace(".fastq.1.gz", ".fastq.2.gz")
            statement_entry = " %s,%s" % (p1, p2)
        # single end reads
        else:
            statement_entry = " %s" % p1

        # build final statement
        statement = metaphlan_statement + statement_entry
        
        return statement

def extract_taxonomy_level(infile, outfile, level):
    """
    Extract specific taxonomy level from MetaPHlAn output
    level should be one of: k__, p__, c__, o__, f__, g__, s__
    """
    statement = '''grep -E "%(level)s|clade" %(infile)s | 
                  grep "%(level)s" | 
                  sed 's/^.*%(level)s//g' > %(outfile)s
    ''' % {'level': level, 'infile': infile, 'outfile': outfile}
    
    return statement

def merge_tables(infiles, outfile):
    """
    Merge multiple MetaPHlAn profile files
    """
    infiles_str = " ".join(infiles)
    statement = '''merge_metaphlan_tables.py 
                  %(infiles)s
                  > %(outfile)s
    ''' % {'infiles': infiles_str, 'outfile': outfile}
    
    return statement