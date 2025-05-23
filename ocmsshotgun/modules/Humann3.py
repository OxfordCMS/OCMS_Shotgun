# Module to run Humann3

from cgatcore import pipeline as P
import os
import ocmstoolkit.modules.Utility as Utility
from importlib.metadata import version, PackageNotFoundError
from packaging.version import Version
from packaging.specifiers import SpecifierSet

class Humann3(Utility.BaseTool):
    '''
    Humann3 tool class, contains methods to building commandline statement
    to run humann3 and post processing of humann3 outputs
    '''
    def __init__(self, infile, outfile, taxonomic_profile=None, **PARAMS):
        super().__init__(infile, outfile, **PARAMS)

        self.taxonomic_profile = taxonomic_profile

        # Detect installed versions
        try:
            # Wrap version string in Version() to ensure proper semantic comparison,
            # including handling of versions like "3.10.1" or pre-releases like "4.0b1"
            self.metaphlan_version = Version(version("metaphlan"))
            self.humann_version = Version(version("humann"))
        except PackageNotFoundError as e:
            raise ImportError(f"Required tool not found: {e.name}")
        
        # Define version requirements
        metaphlan_required = SpecifierSet(">=4.0")
        humann_required = SpecifierSet(">=3.9")

        # Versions sanity check
        if self.metaphlan_version in metaphlan_required and self.humann_version not in humann_required:
            raise ValueError(
                f"MetaPhlAn v{self.metaphlan_version} requires HUMAnN >=3.9, "
                f"but found HUMAnN v{self.humann_version}"
            )

    def concat_fastqs(self, fastn_obj):
        '''
        Method to concatenate fastq files - symlink if there is
        no paired one
        '''
        if not fastn_obj.fastn2:
            statement = f"ln -s {self.infile} {self.outfile}"
        else:
            fastqs = [fastn_obj.fastn1, fastn_obj.fastn2, fastn_obj.fastn3]
            fastqs = [x for x in fastqs if os.path.exists(x)]
            fastqs = ' '.join(fastqs)
            statement = f"cat {fastqs} > {self.outfile}"
        return(statement)
            
    def build_statement(self, fastn_obj):
        '''
        build statement for executing humann3
        '''
        db_metaphlan_path = self.PARAMS["humann3_db_metaphlan_path"]
        db_metaphlan_id = self.PARAMS["humann3_db_metaphlan_id"]
        db_nucleotide = self.PARAMS["humann3_db_nucleotide"]
        db_protein = self.PARAMS['humann3_db_protein']
        search_mode = self.PARAMS['humann3_search_mode']
        options = self.PARAMS["humann3_options"]
        threads = self.PARAMS["humann3_job_threads"]

        # make sure system requreiments not set outside of 
        # job_options and job_memory
        assert "--threads" not in options, (
            "Do not set --threads: Humann3 multi-threading is already set with job_memory and job_threads"
        )
        
        # the option to add a taxonomic profile for restricting mapping
        if self.taxonomic_profile:
            options = '--taxonomic-profile %s' % self.taxonomic_profile \
                + ' ' + options

        #prefix = fastn_obj.prefix
        statement = ("humann"
                     f" --input {fastn_obj.fastn1}"
                     f" --output {self.outdir}"
                     f" --nucleotide-database {db_nucleotide}"
                     f" --protein-database {db_protein}"
                     f" --search-mode {search_mode}"
                     f" --threads {threads}"
                     f" --metaphlan-options \"--index {db_metaphlan_id}"
                     f" --bowtie2db={db_metaphlan_path}\""
                     f" {options} 2> {self.outdir}/{fastn_obj.prefix}.log")
        
        return statement

    def post_process(self):

        prefix = P.snip(os.path.basename(self.outfile), '_pathcoverage.tsv.gz')
        if self.taxonomic_profile:
            options = ""
        else:
            metaphlan_bugs_list = (f"{self.outdir}/{prefix}_humann_temp/"
                                   f"{prefix}_metaphlan_bugs_list.tsv")
            options = (f" gzip {metaphlan_bugs_list} &&"
                       f" mv {metaphlan_bugs_list}.gz {self.outdir} &&")

        humann_log = (f"{self.outdir}/{prefix}_humann_temp/"
                      f"{prefix}.log")
        humann_tmp = f"{self.outdir}/{prefix}_humann_temp"
        statement =  (
            f"mv {humann_log} {self.outdir} &&"
            f" gzip {self.outdir}/{prefix}_pathcoverage.tsv &&"
            f" gzip {self.outdir}/{prefix}_pathabundance.tsv &&" 
            f" gzip {self.outdir}/{prefix}_genefamilies.tsv &&"
            f" {options}"
            f" tar -zcvf {humann_tmp}.tar.gz {humann_tmp} &&"
            f" rm -rf {humann_tmp}"
        )
        
        remove_humann_temp = self.PARAMS.get('humann_remove_humann_temp',True)
        if remove_humann_temp:
            statement = statement + f' && rm -f {humann_tmp}.tar.gz'
        return statement

