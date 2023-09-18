# Module to run Kraken2 and Bracken on paired end reads

import sys
import subprocess
import os
import re
import glob
from pathlib import Path
from ruffus import *
from cgatcore import pipeline as P

class kraken2():
    
    def run(infile, outfile, **PARAMS):
        '''classify reads with kraken2
        '''
        # Note that at the moment I only deal with paired-end
        # reads
        p1 = infile
        p2 = p1.replace(".fastq.1.gz", ".fastq.2.gz")

        prefix = P.snip(outfile, ".k2.report.tsv")
    
        db = PARAMS.get("kraken2_db")
        job_threads = PARAMS.get("kraken2_job_threads")
        job_memory = PARAMS.get("kraken2_job_mem")
        options = PARAMS.get("kraken2_options")
    
        statement = '''kraken2
                   --db %(db)s
                   --output %(prefix)s.classified.tsv
                   --paired
                   --report %(prefix)s.k2.report.tsv
                   --use-names
                   --threads %(job_threads)s
                   --gzip-compressed %(p1)s %(p2)s
                   %(options)s; 
                   gzip %(prefix)s.classified.tsv
                '''

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
