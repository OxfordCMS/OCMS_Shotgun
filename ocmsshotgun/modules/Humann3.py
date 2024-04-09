# Module to run Humann3

import os
from cgatcore import pipeline as P

class humann3(object):
    def __init__(self, infile, outfile, taxonomic_profile=None, **PARAMS):
        self.infile = infile
        self.prefix = P.snip(outfile, '_pathcoverage.tsv.gz', strip_path=True)
        self.outpath = os.path.dirname(os.path.abspath(outfile))
        self.PARAMS = PARAMS
        self.taxonomic_profile = taxonomic_profile
        
    def buildStatement(self):
        infile = self.infile
        prefix = self.prefix
        outpath = self.outpath
        db_metaphlan_path = self.PARAMS["humann3_db_metaphlan_path"]
        db_metaphlan_id = self.PARAMS["humann3_db_metaphlan_id"]
        db_nucleotide = self.PARAMS["humann3_db_nucleotide"]
        db_protein = self.PARAMS['humann3_db_protein']
        search_mode = self.PARAMS['humann3_search_mode']
        options = self.PARAMS["humann3_options"]

        # the option to add a taxonomic profile for restricting mapping
        if self.taxonomic_profile:
            options = '--taxonomic-profile %s' % self.taxonomic_profile \
                + ' ' + options
        
        statement = ("humann"
                     " --input %(infile)s"
                     " --output %(outpath)s"
                     " --nucleotide-database %(db_nucleotide)s"
                     " --protein-database %(db_protein)s"
                     " --search-mode %(search_mode)s"
                     " --metaphlan-options  \"--index %(db_metaphlan_id)s --bowtie2db=%(db_metaphlan_path)s\""
                     " %(options)s 2> %(outpath)s/%(prefix)s.log" % locals())
        
        return statement

    def postProcess(self):
        prefix = self.prefix
        outpath = self.outpath

        if self.taxonomic_profile:
            options = ""
        else:
            options = (" gzip %(outpath)s/%(prefix)s_humann_temp/%(prefix)s_metaphlan_bugs_list.tsv &&"
                       " mv %(outpath)s/%(prefix)s_humann_temp/%(prefix)s_metaphlan_bugs_list.tsv.gz %(outpath)s &&" % locals())
        
        statement =  ("mv %(outpath)s/%(prefix)s_humann_temp/%(prefix)s.log %(outpath)s &&"
                      " gzip %(outpath)s/%(prefix)s_pathcoverage.tsv &&"
                      " gzip %(outpath)s/%(prefix)s_pathabundance.tsv &&" 
                      " gzip %(outpath)s/%(prefix)s_genefamilies.tsv &&"
                      " %(options)s"
                      " tar -zcvf %(outpath)s/%(prefix)s_humann_temp.tar.gz %(outpath)s/%(prefix)s_humann_temp &&"
                      " rm -rf %(outpath)s/%(prefix)s_humann_temp" % locals())

        return statement
    
    def run(self):
        '''functional profile with humann3
        '''
        
        statement = self.buildStatement() + ' && ' + self.postProcess()

        return statement
