# Module to concatenate paired end fastq files to one file

import os
from cgatcore import pipeline as P

class concatFastq():
        
    def run(infiles, outfile):
        """Expects paired end fastq files to be in format fastq.1.gz, fastq.2.gz
        """

        infiles = ' '.join(infiles)
        tempfile = os.path.splitext(outfile)[0]
        
        # concatenate files together
        statement = '''zcat %(infiles)s >> %(tempfile)s && gzip %(tempfile)s'''
        P.run(statement)
