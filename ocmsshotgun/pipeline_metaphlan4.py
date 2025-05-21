import sys
import os
import re
import glob
from pathlib import Path
from ruffus import *
from cgatcore import pipeline as P
import ocmsshotgun.modules.Utility as utility

PARAMS = P.get_parameters(["pipeline.yml"])

FASTQ1 = utility.check_input()

SEQUENCEFILES = ("*.fastq.1.gz")

SEQUENCEFILES_REGEX = regex(
    r"(\S+).(fastq.1.gz)")

@follows(mkdir("metaphlan.dir"))
@transform(SEQUENCEFILES,
          SEQUENCEFILES_REGEX,
          r"metaphlan.dir/\1.metaphlan_profile.txt")
def runMetaphlan(infile, outfile):
    read1 = infile
    read2 = infile.replace(".1.gz", ".2.gz")
    
    sample = os.path.basename(P.snip(infile, ".fastq.1.gz"))
    db_path = PARAMS["metaphlan_db"]  # Using the full path directly
    
    statement = '''metaphlan {read1},{read2} 
                  --input_type fastq 
                  --bowtie2db {db_path}
                  --nproc {threads}
                  --index {db_ver}
                  --bowtie2out metaphlan.dir/{sample}.bowtie2.bz2
                  --tax_lev {tax_lev}
                  -o {outfile}
    '''.format(read1=read1,
               read2=read2,
               db_path=db_path,
               db_ver=PARAMS["metaphlan_database_version"],
               sample=sample,
               threads=PARAMS["metaphlan_threads"],
               tax_lev=PARAMS["metaphlan_tax_lev"],
               outfile=outfile)
    
    if PARAMS["metaphlan_options"]:
        statement += ' ' + PARAMS["metaphlan_options"]
               
    P.run(statement,
          job_threads=PARAMS["metaphlan_threads"],
          job_memory=PARAMS["metaphlan_memory"],
          job_options=PARAMS["metaphlan_cluster_options"])

@follows(runMetaphlan)
@merge(runMetaphlan, "metaphlan.dir/merged_abundance_table.txt")
def mergeMetaphlanTables(infiles, outfile):
    infiles_str = " ".join(infiles)
    
    statement = '''merge_metaphlan_tables.py 
                  {}
                  > {}
    '''.format(infiles_str, outfile)
    
    P.run(statement)

@follows(mkdir("species.dir"))
@transform(mergeMetaphlanTables,
          regex(r"metaphlan.dir/merged_abundance_table.txt"),
          r"species.dir/species_abundance.txt")
def extractSpeciesAbundance(infile, outfile):
    statement = '''grep -E "s__|clade" {} | 
                  grep "s__" | 
                  sed 's/^.*s__//g' > {}
    '''.format(infile, outfile)
    
    P.run(statement)

@follows(mkdir("taxonomy.dir"))
@transform(mergeMetaphlanTables,
          regex(r"metaphlan.dir/merged_abundance_table.txt"),
          [r"taxonomy.dir/kingdom.txt",
           r"taxonomy.dir/phylum.txt",
           r"taxonomy.dir/class.txt",
           r"taxonomy.dir/order.txt",
           r"taxonomy.dir/family.txt",
           r"taxonomy.dir/genus.txt",
           r"taxonomy.dir/species.txt",
           r"taxonomy.dir/strain.txt"])
def extractTaxonomyLevels(infile, outfiles):
    levels = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__','t__']
    for level, outfile in zip(levels, outfiles):
        statement = '''grep -E "{level}|clade" {} | 
                      grep "{level}" | 
                      sed 's/^.*{level}//g' > {}
        '''.format(level=level, infile=infile, outfile=outfile)
        
        P.run(statement)

@follows(extractSpeciesAbundance,
         extractTaxonomyLevels)
def full():
    pass

def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
