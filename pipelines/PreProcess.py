# Module for generic shotgun preprocessing steps

import os,re
from cgatcore import pipeline as P
from cgatcore import iotools as IOTools
from cgatcore import experiment as E


def removeContaminants(in_fastn, out_fastn, method, **PARAMS):
    """
    Remove sequences of non-microbiome origin
    """

    if method == "sortmerna":
        method = sortMeRNA()
    else:
        raise ValueError("Method {} not implemented".format(method))

    method.run(in_fastn, out_fastn, method, **PARAMS)



class matchReference(object):
    """
    Base class for generating run statements to match mWGS reads to 
    reference sequences. Intended to work with single, paired, or
    paired + singleton fastn files. 

    Some options are  assumed to be passed via kwargs, as this and 
    inherited classes are writtento work with a PARAMS dict 
    generated from a pipeline.yml config file.

    ** Options:
    fn_suffix - option to pass something other than .fastq.1.gz

    


    """

    def __init__(self, fastn1, outfile, **PARAMS):
        self.outdir = os.path.dirname(outfile)

        self.fastn1 = fastn1

        # Assume that files are fastq and end .fastq.1.gz
        if PARAMS.get('fn_suffix', None):
            self.fn_suffix = PARAMS.get('fn_suffix')
        else:
            assert self.fastn1.endswith('.fastq.1.gz'), \
                'Please supply fastn suffix using keyword fn_suffix='
            self.fn_suffix = '.fastq.1.gz'

        # Find mate pair file
        n = self.fn_suffix.rfind('1')
        fn2_suffix =  self.fn_suffix[:n] + '2' + self.fn_suffix[n+1:]
        fastn2 = P.snip(self.fastn1, self.fn_suffix) + fn2_suffix

        if os.path.exists(fastn2):
            self.fastn2 = fastn2
        else:
            self.fastn2 = None

        # Find singleton file
        fn3_suffix =  self.fn_suffix[:n] + '3' + self.fn_suffix[n+1:]
        fastn3 = P.snip(self.fastn1, self.fn_suffix) + fn3_suffix

        if os.path.exists(fastn3):
            assert self.fastn2, "Can't have singletons without mate pairs"
            self.fastn3 = fastn3
        else:
            self.fastn3 = None

    def buildStatement(self, *args, **kwargs):
        return ""

    def postProcess(self, *args, **kwargs):
        return ""

    def run(self, *args, **PARAMS):
        
        # Custom command to run reference matching tool.
        statement, run_options = self.buildStatement(**PARAMS)

        # Logging
        runfiles = '\t'.join([os.path.basename(x) for x in (self.fastn1, \
                                                            self.fastn2, \
                                                            self.fastn3) if x])
        E.info("Running sortMeRNA for files: {}".format(runfiles))

        P.run(statement, job_options=run_options)

        # Post process results into generic output for downstream tasks.
        statement = self.postProcess(**PARAMS)
        if statement:
            P.run(statement, run_options)



class runSortMeRNA(matchReference):
    """
    Run sortMeRNA. 
    Assumes that reference indexes have been created in advance in a 
    specified location.
    """


    def buildStatement(self, *args, **PARAMS):
        """
        Generate run statement for processing single, paired, or paired
        + singleton samples. 

        Required arguments: 
        index
        reference
        
        """

        run_options = PARAMS["sortmerna_run_options"]
        threads = PARAMS["sortmerna_threads"]
        
        # A comma separated list of references
        references = PARAMS["sortmerna_reference"]
        references = ' --ref '.join(references.split(','))
        # All listed references must be pre-indexed in this location
        index_dir = PARAMS["sortmerna_index"] # Check this isn't automatically passed. 
        
        tmpf = P.get_temp_dir('.')
        tmpf_kvdb = os.path.join(tmpf, 'kvdb')
        tmpf_readb = os.path.join(tmpf, 'readb')

        if not self.fastn2:
            # Run sortMeRNA for single reads
            in_fastn1 = self.fastn1
            in_prefix = P.snip(in_fastn1, self.fn_suffix, strip_path=True)
            out_prefix = os.path.join(self.outdir, in_prefix)

            # Run sortMeRNA for single reads
            statement = ("sortmerna"
                         " --index 0" # skip indexing, assume in idx-dir
                         " --fastx"
                         " --reads %(in_fastn1)s"
                         " --ref %(references)s"
                         " --idx-dir %(index_dir)s" # location of reference indexes
                         " --aligned %(out_prefix)s_aligned" # output location of aligned seq
                         " --other %(out_prefix)s_unaligned" # output location of unalinged seq
                         " --readb %(tmpf_readb)s" # location of tmp file for reads
                         " --kvdb %(tmpf_kvdb)s" # location of tmp file for kv pairs
                         " --threads %(threads)s"
                         " --zip-out" % locals())

        else:
            # Run sortMeRNA for paired reads
            in_fastn1 = self.fastn1
            in_fastn2 = self.fastn2
            in_prefix = P.snip(in_fastn1, self.fn_suffix, strip_path=True)
            out_prefix = os.path.join(self.outdir, in_prefix)
            # Run sortMeRNA for single reads
            statement = ("sortmerna"
                         " --index 0" # skip indexing, assume in idx-dir
                         " --fastx"
                         " --reads %(in_fastn1)s" # First read file
                         " --reads %(in_fastn2)s" # Second read file
                         " --ref %(references)s"
                         " --idx-dir %(index_dir)s" # location of reference indexes
                         " --aligned %(out_prefix)s_aligned" # output location of aligned seq
                         " --other %(out_prefix)s_unaligned" # output location of unalinged seq
                         " --readb %(tmpf_readb)s" # location of tmp file for reads
                         " --kvdb %(tmpf_kvdb)s" # location of tmp file for kv pairs
                         " --paired_in" # If one read is aligned, both are output to aligned file
                         " --out2" # Output paired reads to separate files
                         " --threads %(threads)s"
                         " --zip-out" % locals())

        if self.fastn3 and not PARAMS.get('sortmerna_skip_singletons', False):
            in_fastn3 = self.fastn3
            statement_2 = ("sortmerna"
                           " --index 0" # skip indexing, assume in idx-dir
                           " --fastx"
                           " --reads %(in_fastn3)s"
                           " --idx-dir %(index_dir)s" # location of reference indexes
                           " --ref %(references)s"
                           " --aligned %(out_prefix)s_aligned_singleton" # output location of aligned seq
                           " --other  %(out_prefix)s_unaligned_singleton" # output location of unalinged seq
                           " --readb %(tmpf_readb)s" # location of tmp file for reads
                           " --kvdb %(tmpf_kvdb)s" # location of tmp file for kv pairs
                           " --threads %(threads)s"
                           " --zip-out" % locals()) 
            
            statement = " && ".join([statement, 
                                     "rm -rf %(tmpf)s/*" % locals(), # location of tmp_readb & kvdb
                                     statement_2,
                                     "rm -rf %(tmpf)s" % locals()])

        return statement, run_options

    def postProcess(self, *args, **PARAMS):
        ''' Rename files output by sortmeRNA to appropriate suffix
        At some point this would be good to become more flexible wrt FQ_SUFFIX'''

        outf_prefix = os.path.join(self.outdir, 
                                   P.snip(self.fastn1, self.fn_suffix, strip_path=True))

        # rename fastn1 files
        os.rename(outf_prefix + '_aligned_fwd.fq.gz', outf_prefix + '_rRNA.fastq.1.gz')
        os.rename(outf_prefix + '_unaligned_fwd.fq.gz', outf_prefix + '_rRNAremoved.fastq.1.gz')

        # rename fastn2 files
        if self.fastn2:
            os.rename(outf_prefix + '_aligned_rev.fq.gz', outf_prefix + '_rRNA.fastq.2.gz')
            os.rename(outf_prefix + '_unaligned_rev.fq.gz', outf_prefix + '_rRNAremoved.fastq.2.gz')

        # rename fastn3 files
        if self.fastn3:
            os.rename(outf_prefix + '_aligned_singleton.fq.gz', 
                      outf_prefix + '_rRNA.fastq.3.gz')
            os.rename(outf_prefix + '_unaligned_singleton.fq.gz',
                      outf_prefix + '_rRNAremoved.fastq.3.gz')
        
        return None
