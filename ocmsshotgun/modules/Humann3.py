# Module to run Humann3

import os
from cgatcore import pipeline as P
import ocmstoolkit.modules.Utility as utility

class humann3(utility.metaFastn):
    def __init__(self, infile, outfile, taxonomic_profile=None, **PARAMS):
        super().__init__(infile, outfile, **PARAMS)
        self.prefix = P.snip(outfile, '_pathcoverage.tsv.gz', strip_path=True)
        self.outpath = os.path.dirname(os.path.abspath(outfile))
        self.taxonomic_profile = taxonomic_profile
        
    def buildStatement(self):
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
                     f" --input {self.fastq1}"
                     f" --output {self.outpath}"
                     f" --nucleotide-database {db_nucleotide}"
                     f" --protein-database {db_protein}"
                     f" --search-mode {search_mode}"
                     f" --metaphlan-options  \"--index {db_metaphlan_id}"
                     f" --bowtie2db={db_metaphlan_path}\""
                     f" {options} 2> {self.outpath}/{self.prefix}.log")
        
        return statement

    def postProcess(self):
        prefix = self.prefix
        outpath = self.outpath

        if self.taxonomic_profile:
            options = ""
        else:
            metaphlan_bugs_list = (f"{self.outpath}/{self.prefix}_humann_temp/"
                                   f"{self.prefix}_metaphlan_bugs_list.tsv")
            options = (f" gzip {metaphlan_bugs_list} &&"
                       f" mv {metaphlan_bugs_list}.gz {self.outpath} &&")
        
        humann_log = (f"{self.outpath}/{self.prefix}_humann_temp/"
                      f"{self.prefix}.log")
        humann_tmp = f"{self.outpath}/{self.prefix}_humann_temp"
        statement =  (
            f"mv {humann_log} {self.outpath} &&"
            f" gzip {self.outpath}/{self.prefix}_pathcoverage.tsv &&"
            f" gzip {self.outpath}/{self.prefix}_pathabundance.tsv &&" 
            f" gzip {self.outpath}/{self.prefix}_genefamilies.tsv &&"
            f" {options}"
            f" tar -zcvf {humann_tmp}.tar.gz {humann_tmp} &&"
            f" rm -rf {humann_tmp}"
        )
        
        remove_humann_temp = self.PARAMS.get('humann_remove_humann_temp',True)
        if remove_humann_temp:
            statement = statement + f' && rm -f {humann_tmp}.tar.gz'
        return statement
    
