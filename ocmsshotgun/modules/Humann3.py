# Module to run Humann3

import os
from cgatcore import pipeline as P

class humann3():
    def statement(infile, outfile, **PARAMS):
        prefix = P.snip(os.path.basename(outfile), "_pathcoverage.tsv.gz")
        outpath = os.path.dirname(os.path.abspath(outfile))
                
        db_metaphlan_path = PARAMS.get("humann3_db_metaphlan_path")
        db_metaphlan_id = PARAMS.get("humann3_db_metaphlan_id")
        db_nucleotide = PARAMS.get("humann3_db_nucleotide")
        db_protein = PARAMS.get('humann3_db_protein')
        search_mode = PARAMS.get('humann3_search_mode')
        options = PARAMS.get("humann3_options")
        statement = '''humann
                    --input %(infile)s
                    --output humann3.dir/%(prefix)s
                    --nucleotide-database %(db_nucleotide)s
                    --protein-database %(db_protein)s
                    --metaphlan-options "--index %(db_metaphlan_id)s --bowtie2db=%(db_metaphlan_path)s"
                    %(options)s;
                    gzip %(outpath)s/%(prefix)s_pathcoverage.tsv;
                    gzip %(outpath)s/%(prefix)s_pathabundance.tsv; 
                    gzip %(outpath)s/%(prefix)s_genefamilies.tsv;
                    gzip %(outpath)s/%(prefix)s_humann_temp/%(prefix)s_metaphlan_bugs_list.tsv;
                    mv %(outpath)s/%(prefix)s_humann_temp/%(prefix)s_metaphlan_bugs_list.tsv.gz %(outpath)s;
                    tar -zcvf %(outpath)s/%(prefix)s_humann_temp.tar.gz %(outpath)s/%(prefix)s_humann_temp;
                    rm -r %(outpath)s/%(prefix)s_humann_temp
        '''
        return statement

    def run(infile, outfile, **PARAMS):
        '''functional profile with humann3
        '''
        
        statement = self.statement(infile, outfile, **PARAMS)
        P.run(statement,
              job_memory=PARAMS.get("humann3_job_memory"),
              job_threads=PARAMS.get("humann3_job_threads"),
              job_options=PARAMS.get('humann3_job_options','')
