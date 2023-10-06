# Module for generic shotgun preprocessing steps

import os,sys,re
import shutil
import itertools
import sqlite3
import distutils
import pandas as pd
from cgatcore import pipeline as P
from cgatcore import iotools as IOTools
from cgatcore import experiment as E


class utility():
    def params_setup():

        PARAMS = P.get_parameters(
            ["pipeline.yml"])
        
        # determine the location of input files (reads)
        try:
            PARAMS["location_fastq"]
        except KeyError:
            DATADIR = "."
        else:
            if PARAMS["location_fastq"] == 0:
                DATADIR = "."
            elif PARAMS["location_fastq"] == 1:
                DATADIR = "input.dir"
            else:
                DATADIR = PARAMS["location_fastq"]  # not recommended practice.

        # Check that the input files correspond
        # (this is set up so that input file suffixes can be altered if necessary)
        FASTQ1_SUFFIX = 'fastq.1.gz'
        FASTQ2_SUFFIX = 'fastq.2.gz'

        fq1_regex = re.compile('(\S+).(%s)' % FASTQ1_SUFFIX)
        mask1 = list(map(lambda x: bool(fq1_regex.match(x)), os.listdir(DATADIR)))
        FASTQ1S = [os.path.join(DATADIR, i) \
                   for i in itertools.compress(os.listdir(DATADIR), mask1)]

        if sum(mask1):
            fq2_regex = re.compile('(\S+).(%s)' % FASTQ2_SUFFIX)
            mask2 = list(map(lambda x: bool(fq2_regex.match(x)), os.listdir(DATADIR)))
            fastq2s = [os.path.join(DATADIR, i) \
                       for i in itertools.compress(os.listdir(DATADIR), mask2)]
            if sum(mask2):
                assert sum(mask1) == sum(mask2), 'Not all input files have pairs'
                fq1_stubs = [fq1_regex.match(x).group(1) for x in FASTQ1S]
                fq2_stubs = [fq2_regex.match(x).group(1) for x in fastq2s]
                assert sorted(fq1_stubs) == sorted(fq2_stubs), \
                    "First and second read pair files do not correspond"        
        else:
            raise ValueError("No input files detected... check the file suffixes"
                             " or specify input directory location in config file")

        return PARAMS, FASTQ1S, FASTQ1_SUFFIX, FASTQ2_SUFFIX

    def symlnk(inf, outf):
        try:
            os.symlink(os.path.abspath(inf), outf)
        except OSError as e:
            if e.errno == errno.EEXIST:
                os.remove(outf)
                os.symlink(inf, outf)

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
        self.fastn2 = None
        self.fastn3 = None
        self.outfile = outfile

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
        

        # Find singleton file
        fn3_suffix =  self.fn_suffix[:n] + '3' + self.fn_suffix[n+1:]
        fastn3 = P.snip(self.fastn1, self.fn_suffix) + fn3_suffix

        if os.path.exists(fastn3):
            assert self.fastn2, "Can't have singletons without mate pairs"
            self.fastn3 = fastn3
        
class cdhit(matchReference):

    def __init__(self, fastn1, outfile, fastq1_suffix, fastq2_suffix, **PARAMS):
        # initialize inherited attributes
        super().__init__(fastn1, outfile, **PARAMS)
        self.fastq1_suffix = fastq1_suffix
        self.fastq2_suffix = fastq2_suffix

    def run(self, fastq1, outfile, *args, **PARAMS):
        '''Filter exact duplicates, if specified in config file'''
        FASTQ1_SUFFIX = self.fastq1_suffix
        FASTQ2_SUFFIX = self.fastq2_suffix

        if self.fastn2:
            fastq2 = P.snip(fastq1, FASTQ1_SUFFIX) + FASTQ2_SUFFIX
            outfile1 = P.snip(outfile, '.gz') 
            outfile2 = P.snip(outfile1, '.fastq.1') + '.fastq.2'
            logfile = P.snip(outfile1, '.fastq.1') + '.log'
            cluster_file = P.snip(outfile1, '1') + '*.clstr'
        
            to_filter = PARAMS['cdhit_dedup']
            if to_filter:
                tmpf1 = P.get_temp_filename('.')
                tmpf2 = P.get_temp_filename('.')
                statement = ("zcat %(fastq1)s > %(tmpf1)s &&"
                             " zcat %(fastq2)s > %(tmpf2)s &&"
                             " cd-hit-dup"
                             "  -i %(tmpf1)s"
                             "  -i2 %(tmpf2)s"
                             "  -o %(outfile1)s"
                             "  -o2 %(outfile2)s"
                             "  %(cdhit_options)s"
                             " &> %(logfile)s &&"
                             " gzip %(outfile1)s &&"
                             " gzip %(outfile2)s &&"
                             " gzip %(logfile)s &&"
                             " rm -f %(tmpf1)s &&"
                             " rm -f %(tmpf2)s &&"
                             " rm -f %(cluster_file)s")
            else:
                E.warn('Deduplication step is being skipped for: %s' % fastq1)
                utility.symlnk(fastq1, outfile)
                utility.symlnk(fastq2, outfile2 + '.gz')

        else:
            outfile1 = P.snip(outfile, '.gz') 
            logfile = P.snip(outfile1, '.fastq.1') + '.log'
            cluster_file = P.snip(outfile1, '1') + '*.clstr'
        
            to_filter = PARAMS['preprocess_dedup']
            if to_filter:
                tmpf1 = P.get_temp_filename('.')
                statement = ("zcat %(fastq1)s > %(tmpf1)s"
                             " cd-hit-dup"
                             "  -i %(tmpf1)s"
                             "  -o %(outfile1)s"
                             "  %(cdhit_options)s"
                             " &> %(logfile)s &&"
                             " gzip %(outfile1)s &&"
                             " gzip %(logfile)s &&"
                             " rm -f %(tmpf1)s &&"
                             " rm -f %(cluster_file)s")

            else:
                E.warn('Deduplication step is being skipped for: %s' % fastq1)
                utiltity.symlnk(fastq1, outfile)     

        run_options = PARAMS.get('cdhit_job_options')
        P.run(statement,
              job_threads=PARAMS['cdhit_job_threads'],
              job_memory=PARAMS['cdhit_job_memory'])

class trimmomatic(matchReference):

    def __init__(self, fastn1, outfile, fastq1_suffix, fastq2_suffix, **PARAMS):
        # initialize inherited attributes
        super().__init__(fastn1, outfile,  **PARAMS)
        self.fastq1_suffix = fastq1_suffix
        self.fastq2_suffix = fastq2_suffix

    def run(self, fastq1, outfile1, *args, **PARAMS):
        '''Remove adapters using Trimmomatic'''
        FASTQ1_SUFFIX = self.fastq1_suffix
        FASTQ2_SUFFIX = self.fastq2_suffix

        if self.fastn2:
            fastq2 = P.snip(fastq1, FASTQ1_SUFFIX) + FASTQ2_SUFFIX
            outfile2 = P.snip(outfile1, '.fastq.1.gz') + '.fastq.2.gz'
            outf1_singletons = P.snip(outfile1, '.fastq.1.gz') + '.fastq.1s.gz'
            outf2_singletons = P.snip(outfile1, '.fastq.1.gz') + '.fastq.2s.gz'
            outf_singletons = P.snip(outfile1, '.fastq.1.gz') + '.fastq.3.gz'
            logfile = P.snip(outfile1, '.fastq.1.gz') + '.trim.log'
            logfile2 = P.snip(outfile1, '.fastq.1.gz') + '.log'
                        
            statement = ("java -Xmx5g -jar %(trimmomatic_jar_path)s PE"
                         " -threads %(trimmomatic_n_threads)s"
                         " -phred%(phred_format)s"
                         " -trimlog %(logfile)s"
                         " %(fastq1)s" # input read 1
                         " %(fastq2)s" # input read 2
                         " %(outfile1)s" # output read 1
                         " %(outf1_singletons)s" # output unpaired read 1
                         " %(outfile2)s" # output read 2
                         " %(outf2_singletons)s" # output unpaired read 2
                         " ILLUMINACLIP:"
                         "%(trimmomatic_adapters)s:"
                         "%(trimmomatic_seed_mismatches)s:"
                         "%(trimmomatic_score_palendromic)s:"
                         "%(trimmomatic_score_simple)s:"
                         "%(trimmomatic_min_adapter_len)s:"
                         "%(trimmomatic_keep_both_reads)s"
                         " LEADING:%(trimmomatic_quality_leading)s"
                         " TRAILING:%(trimmomatic_quality_trailing)s"
                         " MINLEN:%(trimmomatic_minlen)s"
                         " &> %(logfile2)s &&"
                         " gzip -f %(logfile)s &&"
                         " cat %(outf1_singletons)s %(outf2_singletons)s "
                         "  > %(outf_singletons)s &&"
                         " rm -f %(outf1_singletons)s && rm -f %(outf2_singletons)s")

        else:
            logfile = P.snip(outfile1, '.fastq.1.gz') + '.trim.log'
            logfile2 = P.snip(outfile1, '.fastq.1.gz') + '.log'
        
            statement = ("java -Xmx5g -jar %(trimmomatic_jar_path)s PE"
                         " -threads %(trimmomatic_n_threads)s"
                         " -phred%(phred_format)s"
                         " -trimlog %(logfile)s"
                         " %(fastq1)s" # input read 1
                         " %(outfile1)s" # output read 1
                         " ILLUMINACLIP:"
                         "%(trimmomatic_adapters)s:"
                         "%(trimmomatic_seed_mismatches)s:"
                         "%(trimmomatic_score_palendromic)s:"
                         "%(trimmomatic_score_simple)s"
                         "%(trimmomatic_min_adapter_len)s:"
                         "%(trimmomatic_keep_both_reads)s"
                         " LEADING:%(trimmomatic_quality_leading)s"
                         " TRAILING:%(trimmomatic_quality_trailing)s"
                         " MINLEN:%(trimmomatic_minlen)s"
                         " &> %(logfile2)s &&"
                         " gzip -f %(logfile)s")

        run_options = PARAMS.get('trimmomatic_job_options')
        P.run(statement, 
              job_threads=PARAMS['trimmomatic_job_threads'],
              job_memory=PARAMS['trimmomatic_job_memory'])

class runSortMeRNA(matchReference):
    """
    Run sortMeRNA. 
    Assumes that reference indexes have been created in advance in a 
    specified location.
    """

    def __init__(self, fastn1, outfile, **PARAMS):
        super().__init__(fastn1, outfile, **PARAMS)
        self.sortmerna_options = PARAMS['sortmerna_options']
        self.sortmerna_index = PARAMS['sortmerna_index']
        self.sortmerna_reference = PARAMS['sortmerna_reference']
        self.sortmerna_skip_singletons = PARAMS.get('sortmerna_skip_singletons', False)
        
    def buildStatement(self, *args, **PARAMS):
        """
        Generate run statement for processing single, paired, or paired
        + singleton samples. 

        Required arguments: 
        index
        reference
        
        """
        job_options = PARAMS.get("sortmerna_job_options")
        job_threads = PARAMS.get("sortmerna_job_threads")
        job_memory = PARAMS.get("sortmerna_job_memory")
        sortmerna_options = self.sortmerna_options

        # A comma separated list of references
        references = self.sortmerna_reference
        references = ':'.join(references.split(','))
        # All listed references must be pre-indexed in this location
        index_dir = self.sortmerna_index # Check this isn't automatically passed. 
        
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
                         #" --index 0" # skip indexing, assume in idx-dir
                         " --fastx"
                         " --reads %(in_fastn1)s"
                         " --ref %(references)s"
                         " %(index_dir)s" # location of reference indexes
                         " --aligned %(out_prefix)s_aligned" # output location of aligned seq
                         " --other %(out_prefix)s_unaligned" # output location of unalinged seq
                         " --readb %(tmpf_readb)s" # location of tmp file for reads
                         " --kvdb %(tmpf_kvdb)s" # location of tmp file for kv pairs
                         " --threads %(job_threads)s"
                         " --zip-out"
                         " %(sortmerna_options)" % locals())

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
                         " --threads %(job_threads)s"
                         " --zip-out"
                         " %(sortmerna_options)s" % locals())

        if self.fastn3 and not self.sortmerna_skip_singletons:
            in_fastn3 = self.fastn3
            statement_2 = ("sortmerna"
                           # " --index 0" # skip indexing, assume in idx-dir
                           " --fastx"
                           " --reads %(in_fastn3)s"
                           " --idx-dir %(index_dir)s" # location of reference indexes
                           " --ref %(references)s"
                           " --aligned %(out_prefix)s_aligned_singleton" 
                           " --other  %(out_prefix)s_unaligned_singleton"
                           " --readb %(tmpf_readb)s" # location of tmp file for reads
                           " --kvdb %(tmpf_kvdb)s" # location of tmp file for kv pairs
                           " --threads %(job_threads)s"
                           " --zip-out"
                           " %(sortmerna_options)s" % locals()) 
            
            statement = " && ".join([statement, 
                                     "rm -rf %(tmpf)s/*" % locals(), # location of tmp_readb & kvdb
                                     statement_2,
                                     "rm -rf %(tmpf)s" % locals()])
        else:
            statement = " && ".join([statement, 
                                     "rm -rf %(tmpf)s" % locals()])

        return statement, job_threads, job_memory
    
    def run(self, *args, **PARAMS):
        # Custom command to run reference matching tool.
        (statement, job_threads, job_memory) = self.buildStatement(**PARAMS)

        # Logging
        runfiles = '\t'.join([os.path.basename(x) for x in (self.fastn1, \
                                                            self.fastn2, \
                                                            self.fastn3) if x])
        E.info("Running sortMeRNA for files: {}".format(runfiles))

        P.run(statement, 
              job_threads=job_threads,
              job_memory=job_memory)

        # Post process results into generic output for downstream tasks.
        statement = self.postProcess(**PARAMS)
        if statement:
            print(statement)
#            P.run(statement, threads, job_memory)

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
        if self.fastn3 and not PARAMS.get('sortmerna_skip_singletons', False):
            os.rename(outf_prefix + '_aligned_singleton.fq.gz', 
                      outf_prefix + '_rRNA.fastq.3.gz')
            os.rename(outf_prefix + '_unaligned_singleton.fq.gz',
                      outf_prefix + '_rRNAremoved.fastq.3.gz')
        
        return None


class createSortMeRNAOTUs(runSortMeRNA):
    """
    Tweak SortMeRNA code to create closed reference taxonomic classification

    NOTE: This assumes that --otu_map True and -blast '1 cigar qcov' are
    passed to sortmerna_options in pipeline config.

    NOTE: This also ignores singletons by default
    """

    def __init__(self, fastn1, outfile, **PARAMS):
        super().__init__(fastn1, outfile, **PARAMS)
        self.sortmerna_options = PARAMS['sortmerna_otu_options']
        self.sortmerna_index = PARAMS['sortmerna_otu_index']
        self.sortmerna_reference = PARAMS['sortmerna_otu_reference']
        self.sortmerna_skip_singletons = True

        assert re.search('--otu_map', self.sortmerna_options), \
            'Need to set --otu_map in sortmerna'
        
        # output of --otu_map is hard coded as otu_map.txt
        # therefore writing everything to a temporary outfile
        tmpf = P.get_temp_dir(self.outdir)
        self.outdir = tmpf
        self.outfile = outfile

    def postProcess(self, *args, **PARAMS):
        '''Rename files output by sortmerna, including otu_map table.''' 

        tmpf_prefix = os.path.join(self.outdir, 
                                   P.snip(self.fastn1, self.fn_suffix, strip_path=True))
        outf_prefix = os.path.join(os.path.dirname(self.outdir),
                                   P.snip(self.fastn1, self.fn_suffix, strip_path=True))

        # rename fastn1 files
        os.rename(tmpf_prefix + '_aligned_fwd.fq.gz', outf_prefix + '_rRNA.fastq.1.gz')
        os.rename(tmpf_prefix + '_unaligned_fwd.fq.gz', outf_prefix + '_rRNAremoved.fastq.1.gz')

        # rename fastn2 files
        if self.fastn2:
            os.rename(tmpf_prefix + '_aligned_rev.fq.gz', outf_prefix + '_rRNA.fastq.2.gz')
            os.rename(tmpf_prefix + '_unaligned_rev.fq.gz', outf_prefix + '_rRNAremoved.fastq.2.gz')

        # rename 'denovo' otus
        if re.search('de_novo_otu', self.sortmerna_options):
            os.rename(tmpf_prefix + '_aligned_denovo_fwd.fq.gz', outf_prefix + '_denovo.fastq.1.gz')
            os.rename(tmpf_prefix + '_aligned_denovo_rev.fq.gz', outf_prefix + '_denovo.fastq.2.gz')

        # rename log file 
        os.rename(tmpf_prefix + '_aligned.log', outf_prefix + '.log')

        # rename otu_map.txt
        otu_file = os.path.join(self.outdir, 'otu_map.txt')
        os.rename(otu_file, outf_prefix + '_otu_map.txt')

        shutil.rmtree(self.outdir)
        
        return None

class bmtagger(matchReference):

    def __init__(self, fastn1, outfile, **PARAMS):
        # initialize inherited attributes
        super().__init__(fastn1, outfile, **PARAMS)

    def run(fastq1, outfile, **PARAMS):
        '''Remove host contamination using bmtagger'''
        job_options = PARAMS.get("bmtagger_job_options")
        outf_host = P.snip(outfile, '_dehost.fastq.1.gz') + '_host.txt'
        outf_host_stub = P.snip(outf_host, '.txt') + '_toremove'

        # Currently disabled. Has no effect. See drop_fastq.py
        # # Whether to keep pair if a read is identified as host.
        # if PARAMS['bmtagger_keep_pairs']:
        #     keep_pairs = True
        #     E.info("BMTagger: reads with a pair identified as host will be"
        #            " discarded")
        # else:
        #     keep_pairs = False
        #     E.info("BMTagger: reads with a pair identified as host will be"
        #            " kept as singletons (assuming they are not also identified"
        #            " as host)")

        job_options = PARAMS.get('bmtagger_job_options')

        if self.fastn2:
            fastq2 = P.snip(fastq1, '.1.gz') + '.2.gz'
            fastq3 = P.snip(fastq1, '.1.gz') + '.3.gz'
            
            to_remove_paired =  P.get_temp_filename('.')            
            to_remove_singletons = P.get_temp_filename('.')
        
            # In some cases, it may be desirable to screen against multiple hosts.
            indexes = zip(PARAMS['bmtagger_bitmask'].split(','),
                          PARAMS['bmtagger_srprism'].split(','))
            for n, indexes in enumerate(indexes, 1):
                n = str(n)
                bitmask, srprism = indexes

                # Screen the paired reads, then singletons
                tmpdir1 = P.get_temp_dir('.')
                tmpdir2 = P.get_temp_dir('.')
            
                tmpf1 = P.get_temp_filename('.')
                tmpf2 = P.get_temp_filename('.')
                tmpf3 = P.get_temp_filename('.')
    
                # bmtagger truncates fasta headers...  sed 's/[[:space:]]\+/__/g'
                # It won't accept... sed 's|[[:space:]].*$|/1|'
                # It also fails if fastq1 header differs from fastq2
                statement1 = ("zcat %(fastq1)s > %(tmpf1)s &&"
                              " zcat %(fastq2)s > %(tmpf2)s &&"
                              " bmtagger.sh"
                              "  -b %(bitmask)s"
                              "  -x %(srprism)s"
                              "  -T %(tmpdir1)s"
                              "  -q1" # Input is fastq
                              "  -1 %(tmpf1)s"
                              "  -2 %(tmpf2)s"
                              "  -o %(outf_host_stub)s_paired%(n)s"
                              "  &> %(outfile)s.log &&"
                              " cat %(outf_host_stub)s_paired%(n)s"
                              "  >> %(to_remove_paired)s &&"
                              " rm -rf %(tmpdir1)s %(tmpf1)s %(tmpf2)s"
                              "  %(outf_host_stub)s_paired%(n)s")

                # Screen the singletons
                if IOTools.open_file(fastq3).read(1):
                    statement2 = ("zcat %(fastq3)s > %(tmpf3)s &&"
                                  " bmtagger.sh"
                                  "  -b %(bitmask)s"
                                  "  -x %(srprism)s"
                                  "  -T %(tmpdir2)s"
                                  "  -q1" # Input is fastq
                                  "  -1 %(tmpf3)s"
                                  "  -o %(outf_host_stub)s_singletons%(n)s"
                                  " &>> %(outfile)s.log &&"
                                  " cat %(outf_host_stub)s_singletons%(n)s"
                                  "  >> %(to_remove_singletons)s &&"
                                  " rm -rf %(tmpdir2)s %(tmpf3)s"
                                  "  %(outf_host_stub)s_singletons%(n)s")
                else:
                    statement2 = ("touch  %(to_remove_singletons)s &&"
                                  " rm -rf %(tmpdir2)s %(tmpf3)s")

                statement = " && ".join([statement1, statement2])

                P.run(statement, 
                      job_threads=PARAMS['bmtagger_job_threads'],
                      job_memory=PARAMS['bmtagger_job_memory'])
            
            # Drop host contaminated reads
            # A hack due to the fact that BMTagger truncates fastq identifiers
            # TO DO: Look at bmtagger/.../bin/extract_fullseq
            drop_script = os.path.join(os.path.splitext(__file__)[0],
                                           'drop_fastqs.py')
        
            fastq1_out = outfile
            fastq2_out = P.snip(outfile, '.1.gz') + '.2.gz'
            fastq3_out = P.snip(outfile, '.1.gz') + '.3.gz'

            fastq1_host = P.snip(outfile, '_dehost.fastq.1.gz') + '_host.fastq.1.gz'
            fastq2_host = P.snip(outfile, '_dehost.fastq.1.gz') + '_host.fastq.2.gz'
            fastq3_host = P.snip(outfile, '_dehost.fastq.1.gz') + '_host.fastq.3.gz'
    
            statement = ("python %(drop_script)s"
                         " --fastq1 %(fastq1)s"
                         " --fastq2 %(fastq2)s"
                         " --fastq3 %(fastq3)s"
                         " --to-drop-paired %(to_remove_paired)s"
                         " --to-drop-single %(to_remove_singletons)s"
                         " --fastq-out1 %(fastq1_out)s"
                         " --fastq-out2 %(fastq2_out)s"
                         " --fastq-out3 %(fastq3_out)s"
                         " --fastq-drop1 %(fastq1_host)s"
                         " --fastq-drop2 %(fastq2_host)s"
                         " --fastq-drop3 %(fastq3_host)s"
                         " &>> %(outfile)s.log")

            P.run(statement)
        
            os.unlink(to_remove_paired)
            os.unlink(to_remove_singletons)

        
        else:
            indexes = zip(PARAMS['bmtagger_bitmask'].split(','),
                          PARAMS['bmtagger_srprism'].split(','))
            to_remove = P.get_temp_filename('.')
            
            for n, indexes in enumerate(indexes, 1):
                n = str(n)
                bitmask, srprism = indexes
                # Screen the singletons
                tmpdir1 = P.get_temp_dir('.')
                tmpf = P.get_temp_filename('.')
                
                statement = ("zcat %(fastq1)s > %(tmpf)s &&"
                             " bmtagger.sh"
                             "  -b %(bitmask)s"
                             "  -x %(srprism)s"
                             "  -T %(tmpdir1)s"
                             "  -q1" # Input is fastq
                             "  -1 %(tmpf)s"
                             "  -o %(outf_host_stub)s_%(n)s"
                             "  &>> %(outfile)s.log &&"
                             " cat %(outf_host_stub)s_%(n)s >> %(to_remove)s"
                             " rm -rf %(tmpdir1)s %(tmpf)s %(outf_host_stub)s_%(n)s")
                
                
                P.run(statement, 
                      job_threads=PARAMS['bmtagger_job_threads'],
                      job_memory=PARAMS['bmtagger_job_memory'])
            

            # Drop host contaminated reads
            drop_script = ps.path.join(os.path.splitext(__file__)[0],
                                       'drop_single_fastqs.py')

            fastq_host = P.snip(outfile, '_dehost.fastq.1.gz') + '_host.fastq.1.gz'
            
            statement = ("python %(drop_script)s"
                         " --fastq1 %(fastq1)s"
                         " --to-drop-single %(to_remove)s"
                         " --fastq-out1 %(outfile)s"
                         " --fastq-drop1 %(fastq_host)s"
                         " &>> %(outfile)s.log")
            
            P.run(statement)
            
            os.unlink(to_remove)

class bbtools(matchReference):
    def __init__(self, fastn1, outfile, **PARAMS):
        # initialize inherited attributes
        super().__init__(fastn1, outfile, **PARAMS)

    def run(fastq1, outfile, **PARAMS):
        '''Either softmask low complexity regions, or remove reads with a large
        proportion of low complexity. 

        Uses BBTools scripts bbduk.sh (removal), or bbmask.sh. 

        Entropy is calculated as shannon entropy for of kmers with a specified size
        within a sliding window. Ranges from 0: mask nothing, 0.0001: mask
        homopolymers, 1: mask everything.
        '''

        bb_options= ' '.join(PARAMS['dust_options'].split(','))
        job_options = PARAMS.get('dust_job_options')

        # bbmap assumes the file format based on the output being *fq.gz
        # I can't find any instructions as to how to override this.   
        if self.fastn2:
            fastq2 = P.snip(fastq1, '.1.gz') + '.2.gz'
            fastq3 = P.snip(fastq1, '.1.gz') + '.3.gz'
            
            outfile1 = P.snip(outfile, '.1.gz') + '.1.fq.gz'
            outfile2 = P.snip(outfile, '.1.gz') + '.2.fq.gz'
            outfile3 = P.snip(outfile, '.1.gz') + '.3.fq.gz'

            out_disc1 = P.snip(outfile, '_masked.fastq.1.gz') + '_discarded.fastq.1.fq.gz'
            out_disc2 = P.snip(outfile, '_masked.fastq.1.gz') + '_discarded.fastq.2.fq.gz'
            out_disc3 = P.snip(outfile, '_masked.fastq.1.gz') + '_discarded.fastq.3.fq.gz'
        
            if PARAMS['dust_discard_low_complexity']:
                statement1 = ("bbduk.sh"
                              "  in=%(fastq1)s"
                              "  in2=%(fastq2)s"
                              "  out=%(outfile1)s"
                              "  out2=%(outfile2)s"
                              "  outm=%(out_disc1)s"
                              "  outm2=%(out_disc2)s"
                              "  entropy=%(dust_entropy)s"
                              "  threads=%(dust_job_thread)s"
                              "  %(bb_options)s"
                              "  &> %(outfile)s.log")
                if IOTools.open_file(fastq3).read(1):
                    statement2 = (" bbduk.sh"
                                  "  in=%(fastq3)s"
                                  "  out=%(outfile3)s"
                                  "  outm=%(out_disc3)s"
                                  "  entropy=%(dust_entropy)s"
                                  "  threads=%(dust_job_threads)s"
                                  "  %(bb_options)s"
                                  "  &>> %(outfile)s.log")
                else:
                    statement2 = (" touch %(outfile3)s  %(out_disc3)s")

                statement = " && ".join([statement1, statement2])
            
                P.run(statement, 
                      job_threads=PARAMS['dust_job_threads'],
                      job_memory=PARAMS['dust_job_memory'])

            else:
                statement1 = ("bbmask.sh"
                              "  in=%(fastq1)s"
                              "  out=%(outfile1)s"
                              "  entropy=%(dust_entropy)s"
                              "  threads=%(dust_job_threads)s"
                              "  overwrite=t"
                              "  lowercase=t"
                              "  %(bb_options)s"
                              "  &> %(outfile)s.log &&"
                              " bbmask.sh"
                              "  in=%(fastq2)s"
                              "  out=%(outfile2)s"
                              "  entropy=%(dust_entropy)s"
                              "  threads=%(dust_job_threads)s"
                              "  overwrite=t"
                              "  lowercase=t"
                              "  %(bb_options)s"
                              "  &>> %(outfile)s.log")
                if IOTools.open_file(fastq3).read(1):           
                    statement2 = (" bbmask.sh"
                                  "  in=%(fastq3)s"
                                  "  out=%(outfile3)s"
                                  "  entropy=%(dust_entropy)s"
                                  "  threads=%(dust_job_threads)s"
                                  "  overwrite=t"
                                  "  lowercase=t"
                                  "  %(bb_options)s"
                                  "  &>> %(outfile)s.log")
                else:
                    statement2 = (" touch %(outfile3)s")

                statement = " && ".join([statement1, statement2])

                P.run(statement,
                      job_threads=PARAMS['dust_job_threads'],
                      job_memory=PARAMS['dust_job_memory'])


            # Renaming files because of bbmap idiosyncracies
            of1 = P.snip(outfile1, '.fq.gz') + '.gz'
            of2 = P.snip(outfile2, '.fq.gz') + '.gz'
            of3 = P.snip(outfile3, '.fq.gz') + '.gz'
            os.rename(outfile1, of1)
            os.rename(outfile2, of2)
            os.rename(outfile3, of3)

            if PARAMS['dust_discard_low_complexity']:
                od1 = P.snip(out_disc1, '.fq.gz') + '.gz'
                od2 = P.snip(out_disc2, '.fq.gz') + '.gz'
                od3 = P.snip(out_disc3, '.fq.gz') + '.gz'
                os.rename(out_disc1, od1)
                os.rename(out_disc2, od2)
                os.rename(out_disc3, od3)

        else:
            outfile1 = P.snip(outfile, '.gz') + '.fq.gz'
            out_disc = P.snip(outfile, '_masked.fastq.1.gz') + '_discarded.fastq.1.fq.gz'
        
            if PARAMS['dust_discard_low_complexity']:
                statement = ("bbduk.sh"
                             " in=%(fastq1)s"
                             " out=%(outfile1)s"
                             " outm=%(out_disc)s"
                             " entropy=%(dust_entropy)s"
                             " threads=%(dust_threads)s"
                             " lowercase=t"
                             " %(bb_options)s"
                             " &> %(outfile)s.log")

                P.run(statement, 
                      job_threads=PARAMS['dust_job_threads'],
                      job_memory=PARAMS['dust_job_memory'])

            else:
                statement = ("bbmask.sh"
                             " in=%(fastq1)s"
                             " out=%(outfile1)s"
                             " entropy=%(dust_entropy)s"
                             " threads=%(dust_job_threads)s"
                             " lowercase=t"
                             " %(bb_options)s"
                             " &> %(outfile.log")

                P.run(statement, 
                      job_threads=PARAMS['dust_job_threads'],
                      job_memory=PARAMS['dust_job_memory'])

            os.rename(outfile1, outfile)
            if PARAMS['dust_discard_low_complexity']:
                od1 = P.snip(out_disc, '.fq.gz') + '.gz'
            os.rename(out_disc, od1)
