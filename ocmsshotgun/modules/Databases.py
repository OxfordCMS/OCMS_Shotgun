##############################################
##############################################
##############################################
# Classes and functions for use with
# pipeline_databases.py
##############################################
##############################################
##############################################
import os

def get_mammalian_genomes():
    '''
    get the plain fasta files from the UCSC
    '''
    to_cluster=False
    statement = '''
    wget https://hgdownload.soe.ucsc.edu/goldenPath/%(human_build)s/bigZips/%(human_build)s.fa.gz
    -P genomes/human/%(human_build)s;
    wget https://hgdownload.soe.ucsc.edu/goldenPath/%(mouse_build)s/bigZips/%(mouse_build)s.fa.gz
    -P genomes/mouse/%(mouse_build)s
    '''
    return(statement)

##############################################
##############################################
##############################################

class DB:
    def __init__(self, PARAMS, tool):
        self.PARAMS = PARAMS
        self.gcc_version = PARAMS["gcc"]
        self.python_version = PARAMS["python"]
        self.tool = tool
        self.version = PARAMS[self.tool+"_version"]
        self.outdir = os.path.join(self.tool,
                                   "GCC-"+self.gcc_version,
                                   self.tool+"-"+self.version)

class SRPRISMdb(DB):

    def build_statement(self, infile, outfile):
        outdir = os.path.dirname(outfile)
        build = os.path.basename(infile).replace(".fa.gz", "")
        outprefix = os.path.join(outdir, build + ".srprism")
        statement = '''srprism mkindex -i %(infile)s -o %(outprefix)s 
                         ''' % locals()
        return(statement)

class BMTOOLdb(DB):

    def build_statement(self, infile, outfile):
        self.tool = "bmtool"
        statement = '''bmtool -d <(zcat %(infile)s) -o %(outfile)s 
                         ''' % locals()
        return(statement)

class SORTMERNAdb(DB):

    def get_plain_fasta(self, outfiles):
        # at the moment sortmerna is installed as version 4.3.<number> and so this should
        # work for any version that is installed on BMRC. If sortmerna gets updated then
        # need to change this harcoding.
        outdir = os.path.dirname(outfiles[0])
        statement = '''wget https://github.com/sortmerna/sortmerna/releases/download/v4.3.3/database.tar.gz
                   -O %(outdir)s/database.tar.gz;
                   tar -xzvf %(outdir)s/database.tar.gz -C %(outdir)s;
                   rm -rf %(outdir)s/database.tar.gz
                ''' % locals()
        return(statement)

    def build_statement(self, infile, outfile):

        prefix = os.path.basename(infile).replace(".fasta", "")
        indir = os.path.dirname(infile)
        outdir = os.path.dirname(outfile)
        idx_dir = os.path.join(outdir, prefix)    
        statement = '''mkdir %(idx_dir)s;
                       sortmerna --ref %(infile)s
                                 --index 1
                                 --idx-dir %(idx_dir)s
                    ''' % locals()
        return(statement)

class KRAKEN2db(DB):

    def make_dummy_files(self):
        filenames = self.PARAMS["kraken2_filenames"].split(",")
        outdir = self.outdir

        # I don't think we need the GCC version
        # for kraken2 so would just create 
        # duplicate databases which not desirable
        outdir = os.path.join(os.path.dirname(os.path.dirname(outdir)),
                              os.path.basename(outdir))
        dummies = [
            os.path.join(outdir,
                         x.replace(".tar.gz", ""),
                     ) for x in filenames
        ]
        for d in dummies:
            try:
                os.mkdir(d)
            except FileExistsError:
                continue

        dummy_files = " ".join([x + "/" + os.path.basename(x) + ".input" for x in dummies]        )
        statement = '''touch %(dummy_files)s''' % locals()
        return(statement)

    def build_statement(self, infile, outfile):
        url = "https://genome-idx.s3.amazonaws.com/kraken/"
        outf = infile.replace(".input", ".tar.gz")
        infile = os.path.join(url, os.path.basename(infile).replace(".input", ".tar.gz"))
        outdir = os.path.dirname(outf)
        statement = '''
                   wget %(infile)s -P %(outdir)s;
                   tar -xzvf %(outf)s -C %(outdir)s;
                   touch %(outfile)s 
                   ''' % locals()
        return(statement)

class METAPHLANdb(DB):
    def build_statement(self, infile, outfile):
        metaphlan_version = self.version.split(".")[0]
        outdir = self.outdir.split("/")
        outdir = os.path.join(outdir[0],
                              "python-"+self.python_version,
                              "metaphlan-"+metaphlan_version)
        print(outdir)        

        if metaphlan_version == "3":
            statement = '''
                        metaphlan --install
                                  --index
                                  mpa_v31_CHOCOPhlAn_201901
                                  --bowtie2db 
                                  %(outdir)s
                        ''' % locals()
        elif metaphlan_version == "4":
            statement = '''
                        metaphlan --install
                                  --index
                                  mpa_vOct22_CHOCOPhlAnSGB_202403
                                  --bowtie2db 
                                  %(outdir)s
                        ''' % locals()
        return(statement)
