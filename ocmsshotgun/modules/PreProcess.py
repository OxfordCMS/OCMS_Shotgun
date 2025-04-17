# Module for generic shotgun preprocessing steps

import os,sys,re
import shutil
import itertools
import sqlite3
import distutils
import pandas as pd
import errno
import json
from glob import glob
from cgatcore import pipeline as P
from cgatcore import iotools as IOTools
from cgatcore import experiment as E

import ocmsshotgun.modules.Utility as utility            
class cdhit(utility.matchReference):
        
    def buildStatement(self):
        '''Filter exact duplicates, if specified in config file'''
        
        fastq1 = self.fastq1
        fastq2 = self.fastq2
        sample_out = P.snip(self.outfile, self.fq1_suffix)
        
        outfile1 = P.snip(self.outfile, '.gz')
        logfile = sample_out + '.log'
        cluster_file = sample_out + '*.clstr'
        
        cdhit_options = self.PARAMS['cdhit_options']
        to_filter = self.PARAMS['cdhit_dedup']
        
        if self.fastq2:
            outfile2 = re.sub(self.fq1_suffix, self.fq2_suffix, self.outfile)
            outfile2 = P.snip(outfile2, '.gz')
            
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
                             " rm -f %(cluster_file)s" % locals())
            else:
                E.warn('Deduplication step is being skipped for: %s' % fastq1)
                utility.symlnk(fastq1, self.outfile)
                utility.symlnk(fastq2, outfile2 + '.gz')

        else:
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
                             " rm -f %(cluster_file)s" % locals())

            else:
                E.warn('Deduplication step is being skipped for: %s' % fastq1)
                utility.symlnk(fastq1, outfile)     
        
        return statement
        
class trimmomatic(utility.matchReference):

    def buildStatement(self):
        '''Remove adapters using Trimmomatic'''
        fastq1 = self.fastq1
        fastq2 = self.fastq2
        outfile1 = self.outfile
        sample_out = P.snip(self.outfile, self.fq1_suffix)
        logfile = sample_out + '.trim.log'
        logfile2 = sample_out + '.log'
        
        trimmomatic_jar_path = self.PARAMS["trimmomatic_jar_path"]
        trimmomatic_n_threads = self.PARAMS["trimmomatic_n_threads"]
        # >0.32 determines phred format automatically, here for legacy
        phred_format = '-phred' + str(self.PARAMS.get('phred_format', 33))
        
        trimmomatic_adapters = self.PARAMS["trimmomatic_adapters"]
        trimmomatic_seed_mismatches = self.PARAMS["trimmomatic_seed_mismatches"]
        trimmomatic_score_palendromic = self.PARAMS["trimmomatic_score_palendromic"]
        trimmomatic_score_simple = self.PARAMS["trimmomatic_score_simple"]
        trimmomatic_min_adapter_len = self.PARAMS["trimmomatic_min_adapter_len"]
        trimmomatic_keep_both_reads = self.PARAMS["trimmomatic_keep_both_reads"]
        trimmomatic_quality_leading = self.PARAMS["trimmomatic_quality_leading"]
        trimmomatic_quality_trailing = self.PARAMS["trimmomatic_quality_trailing"]
        trimmomatic_minlen = self.PARAMS["trimmomatic_minlen"]

        if self.fastq2:
            outfile2 = re.sub(self.fq1_suffix, self.fq2_suffix, self.outfile)
            outf1_singletons = sample_out + re.sub("1", "1s", self.fq1_suffix)
            outf2_singletons = sample_out + re.sub("2", "2s", self.fq2_suffix)
            outf_singletons = sample_out + self.fq3_suffix
            
            statement = ("java -Xmx5g -jar %(trimmomatic_jar_path)s PE"
                         " -threads %(trimmomatic_n_threads)s"
                         " %(phred_format)s"
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
                         " rm -f %(outf1_singletons)s && rm -f %(outf2_singletons)s" % locals())

        else:
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
                         " gzip -f %(logfile)s" % locals())
                
        return statement

#def removeContaminants(in_fastq, out_fastq, method, **PARAMS):
#    """
#    Remove sequences of non-microbiome origin
#    """
#
#    if method == "sortmerna":
#        tool = runSortMeRNA(in_fastq, out_fastq, **PARAMS)
#    else:
#        raise ValueError("Method {} not implemented".format(method))
#
#    tool.run()
#    # Post process results into generic output for downstream tasks.
#    tool.postProcess()

    
class runSortMeRNA(utility.matchReference):
    """
    Run sortMeRNA. 
    Assumes that reference indexes have been created in advance in a 
    specified location.
    """

    def __init__(self, fastq1, outfile, **PARAMS):
        super().__init__(fastq1, outfile, **PARAMS)
        self.sortmerna_skip_singletons = self.PARAMS.get('sortmerna_skip_singletons', False)
        
    def buildStatement(self):
        """
        Generate run statement for processing single, paired, or paired
        + singleton samples. 

        Required arguments: 
        index
        reference
        
        """
        sortmerna_options = self.PARAMS.get("sortmerna_options")

        # A comma separated list of references
        references = self.PARAMS.get("sortmerna_reference")
        references = ' --ref '.join(references.split(','))
        # All listed references must be pre-indexed in this location
        index_dir = self.PARAMS.get("sortmerna_index") # Check this isn't automatically passed. 
        
        tmpf = P.get_temp_dir('.')
        tmpf_kvdb = os.path.join(tmpf, 'kvdb')
        tmpf_readb = os.path.join(tmpf, 'readb')

        job_threads = self.PARAMS.get("sortmerna_job_threads")
        if not self.fastq2:
            # Run sortMeRNA for single reads
            in_fastq1 = self.fastq1
            in_prefix = P.snip(in_fastq1, '_deadapt'+self.fq1_suffix, strip_path=True)
            out_prefix = os.path.join(self.outdir, in_prefix)

            # Run sortMeRNA for single reads
            statement = ("sortmerna"
                         #" --index 0" # skip indexing, assume in idx-dir
                         " --fastx"
                         " --reads %(in_fastq1)s"
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
            in_fastq1 = self.fastq1
            in_fastq2 = self.fastq2
            in_prefix = P.snip(in_fastq1, '_deadapt'+self.fq1_suffix, strip_path=True)
            out_prefix = os.path.join(self.outdir, in_prefix)
            # Run sortMeRNA for single reads
            statement = ("sortmerna"
                         " --index 0" # skip indexing, assume in idx-dir
                         " --fastx"
                         " --reads %(in_fastq1)s" # First read file
                         " --reads %(in_fastq2)s" # Second read file
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
        
        if not self.sortmerna_skip_singletons and IOTools.open_file(self.fastq3).read(1):
            in_fastq3 = self.fastq3
            statement_2 = ("sortmerna"
                           # " --index 0" # skip indexing, assume in idx-dir
                           " --fastx"
                           " --reads %(in_fastq3)s"
                           " --idx-dir %(index_dir)s" # location of reference indexes
                           " --ref %(references)s"
                           " --aligned %(out_prefix)s_aligned_singleton" 
                           " --other  %(out_prefix)s_unaligned_singleton"
                           " --readb %(tmpf_readb)s" # location of tmp file for reads
                           " --kvdb %(tmpf_kvdb)s" # location of tmp file for kv pairs
                           " --threads %(job_threads)s"
                           " --zip-out 1"
                           " %(sortmerna_options)s" % locals()) 
            
            statement = " && ".join([statement, 
                                     "rm -rf %(tmpf)s/*" % locals(), # location of tmp_readb & kvdb
                                     statement_2,
                                     "rm -rf %(tmpf)s" % locals()])
        else:
            statement = " && ".join([statement, 
                                     "rm -rf %(tmpf)s" % locals()])

        return statement
        
    def postProcess(self):
        ''' Rename files output by sortmeRNA to appropriate suffix
        At some point this would be good to become more flexible wrt FQ1_SUFFIX'''

        outf_prefix = os.path.join(self.outdir, 
                                   P.snip(self.fastq1, '_deadapt'+self.fq1_suffix, strip_path=True))

        # rename fastq1 files
        os.rename(outf_prefix + '_aligned_fwd.fq.gz', 
                  outf_prefix + '_rRNA'+self.fq1_suffix)
        os.rename(outf_prefix + '_unaligned_fwd.fq.gz', 
                  outf_prefix + '_rRNAremoved'+self.fq1_suffix)

        # rename fastq2 files
        if self.fastq2:
            os.rename(outf_prefix + '_aligned_rev.fq.gz', 
                      outf_prefix + '_rRNA'+self.fq2_suffix)
            os.rename(outf_prefix + '_unaligned_rev.fq.gz', 
                      outf_prefix + '_rRNAremoved'+self.fq2_suffix)

        # rename fastq3 files
        if not self.sortmerna_skip_singletons and IOTools.open_file(self.fastq3).read(1):            
            os.rename(f3_aligned,
                      outf_prefix + '_rRNA' + self.fq3_suffix)

            f3_unaligned = glob(os.path.join(self.outdir, "*unaligned_singleton*"))[0]            
            os.rename(f3_unaligned,
                      outf_prefix +  '_rRNAremoved' + self.fq3_suffix)
        else:
            utility.symlnk(self.fastq3, 
                           os.path.join(outf_prefix + "_rRNAremoved" + self.fq3_suffix))
        return None


class createSortMeRNAOTUs(runSortMeRNA):
    """
    Tweak SortMeRNA code to create closed reference taxonomic classification

    NOTE: This assumes that --otu_map True and -blast '1 cigar qcov' are
    passed to sortmerna_options in pipeline config.

    NOTE: This also ignores singletons by default
    """

    def __init__(self, fastq1, outfile, **PARAMS):
        super().__init__(fastq1, outfile, **PARAMS)
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
    
    def postProcess(self):
        '''Rename files output by sortmerna, including otu_map table.''' 

        tmpf_prefix = os.path.join(self.outdir, 
                                   P.snip(self.fastq1, '_deadapt'+self.fq1_suffix, strip_path=True))
        outf_prefix = os.path.join(os.path.dirname(self.outdir),
                                   P.snip(self.fastq1, '_deadapt'+self.fq1_suffix, strip_path=True))

        # rename fastq1 files
        os.rename(tmpf_prefix + '_aligned_fwd.fq.gz', outf_prefix + '_rRNA'+self.fq1_suffix)
        os.rename(tmpf_prefix + '_unaligned_fwd.fq.gz', outf_prefix + '_rRNAremoved'+self.fq1_suffix)

        # rename fastq2 files
        if self.fastq2:
            os.rename(tmpf_prefix + '_aligned_rev.fq.gz', outf_prefix + '_rRNA'+self.fq2_suffix)
            os.rename(tmpf_prefix + '_unaligned_rev.fq.gz', outf_prefix + '_rRNAremoved'+self.fq2_suffix)

        # rename 'denovo' otus
        if re.search('de_novo_otu', self.sortmerna_options):
            os.rename(tmpf_prefix + '_aligned_denovo_fwd.fq.gz', outf_prefix + '_denovo'+self.fq1_suffix)
            os.rename(tmpf_prefix + '_aligned_denovo_rev.fq.gz', outf_prefix + '_denovo'+self.fq2_suffix)

        # rename log file 
        os.rename(tmpf_prefix + '_aligned.log', outf_prefix + '.log')

        # rename otu_map.txt
        otu_file = os.path.join(self.outdir, 'otu_map.txt')
        os.rename(otu_file, outf_prefix + '_otu_map.txt')

        shutil.rmtree(self.outdir)
        
        return None

class bmtagger(utility.matchReference):

    def buildStatement(self):
        '''Remove host contamination using bmtagger'''
        outf_host = P.snip(self.outfile, '_dehost'+self.fq1_suffix) + '_host.txt'
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

        fastq1 = self.fastq1
        outfile = self.outfile
        bmtagger_exec = self.PARAMS['bmtagger_executable']
        assert bmtagger_exec in ["bmtagger.sh", "bmtagger_mod.sh"], "must specify bmtagger.sh or bmtagger_mod.sh"

        if self.fastq2:
            fastq2 = self.fastq2
            fastq3 = self.fastq3
            
            to_remove_paired =  P.get_temp_filename('.')            
            to_remove_singletons = P.get_temp_filename('.')
            to_remove_tmp = [to_remove_paired, to_remove_singletons]
 
            # In some cases, it may be desirable to screen against multiple hosts.
            indexes = zip(self.PARAMS['bmtagger_bitmask'].split(','),
                          self.PARAMS['bmtagger_srprism'].split(','))
           
            statements=[]
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
                #### UPDATED LOCATION OF bmtagger.sh ####
                statement1 = ("zcat %(fastq1)s > %(tmpf1)s &&"
                              " zcat %(fastq2)s > %(tmpf2)s &&"
                              " %(bmtagger_exec)s"
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
                              "  %(outf_host_stub)s_paired%(n)s" % locals())

                # Screen the singletons
                #### UPDATED LOCATION OF bmtagger.sh ####
                if os.path.exists(self.fastq3) and IOTools.open_file(self.fastq3).read(1):
                    statement2 = ("zcat %(fastq3)s > %(tmpf3)s &&"
                                  " %(bmtagger_exec)s"
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
                                  "  %(outf_host_stub)s_singletons%(n)s" % locals())
                else:
                    statement2 = ("touch  %(to_remove_singletons)s &&"
                                  " rm -rf %(tmpdir2)s %(tmpf3)s" % locals())

                statement = " && ".join([statement1, statement2])
                statements.append(statement)

        else:
            indexes = zip(self.PARAMS['bmtagger_bitmask'].split(','),
                          self.PARAMS['bmtagger_srprism'].split(','))
            to_remove = P.get_temp_filename('.')
            to_remove_tmp = [to_remove]
            
            statements = []
            for n, indexes in enumerate(indexes, 1):
                n = str(n)
                bitmask, srprism = indexes
                # Screen the singletons
                tmpdir1 = P.get_temp_dir('.')
                tmpf = P.get_temp_filename('.')
                #### UPDATED LOCATION OF bmtagger.sh ####
                statement = ("zcat %(fastq1)s > %(tmpf)s &&"
                             " %(bmtagger_exec)s"
                             "  -b %(bitmask)s"
                             "  -x %(srprism)s"
                             "  -T %(tmpdir1)s"
                             "  -q1" # Input is fastq
                             "  -1 %(tmpf)s"
                             "  -o %(outf_host_stub)s_%(n)s"
                             "  &>> %(outfile)s.log &&"
                             " cat %(outf_host_stub)s_%(n)s >> %(to_remove)s"
                             " rm -rf %(tmpdir1)s %(tmpf)s %(outf_host_stub)s_%(n)s" % locals())
                statements.append(statement)
                
        return statements, to_remove_tmp

    def postProcess(self, to_remove_tmp):

        if self.fastq2:
            # Drop host contaminated reads
            # A hack due to the fact that BMTagger truncates fastq identifiers
            # TO DO: Look at bmtagger/.../bin/extract_fullseq
            
            fastq1 = self.fastq1
            fastq2 = self.fastq2
            
            fastq1_out = self.outfile
            fastq2_out = P.snip(self.outfile, self.fq1_suffix) + self.fq2_suffix
            
            fastq1_host = P.snip(self.outfile, '_dehost'+self.fq1_suffix) + '_host'+self.fq1_suffix
            fastq2_host = P.snip(self.outfile, '_dehost'+self.fq1_suffix) + '_host'+self.fq2_suffix
            
            fastq3 = self.fastq3
            fastq3_out = P.snip(self.outfile, self.fq1_suffix) + self.fq3_suffix
            fastq3_host = P.snip(self.outfile, '_dehost'+self.fq1_suffix) + '_host'+self.fq3_suffix
            to_remove_paired = to_remove_tmp[0]
            to_remove_singletons = to_remove_tmp[1]

            statement = ("ocms_shotgun drop_fastqs"
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
                         " &>> %(fastq1_out)s.log" % locals())

            to_unlink = [to_remove_paired, to_remove_singletons]

        else:
            
            fastq1 = self.fastq1
            outfile = self.outfile
            to_remove = to_remove_tmp[0]
            fastq_host = P.snip(self.outfile, '_dehost'+self.fq1_suffix) + '_host'+self.fq1_suffix
            statement = ("ocms_shotgun drop_fastqs"
                         " --fastq1 %(fastq1)s"
                         " --to-drop-single %(to_remove)s"
                         " --fastq-out1 %(outfile)s"
                         " --fastq-drop1 %(fastq_host)s"
                         " &>> %(outfile)s.log" % locals())
        
            P.run(statement)
            
            os.unlink(to_remove)
            to_unlink = [to_remove]

        return statement, to_unlink

class bbtools(utility.matchReference):

    def buildStatement(self):
        '''Either softmask low complexity regions, or remove reads with a large
        proportion of low complexity. 

        Uses BBTools scripts bbduk.sh (removal), or bbmask.sh. 

        Entropy is calculated as shannon entropy for of kmers with a specified size
        within a sliding window. Ranges from 0: mask nothing, 0.0001: mask
        homopolymers, 1: mask everything.
        '''

        # bbmap assumes the file format based on the output being *fq.gz
        # I can't find any instructions as to how to override this.
        entropy = self.PARAMS['dust_entropy']
        bb_options= ' '.join(self.PARAMS['dust_options'].split(','))
        job_threads = self.PARAMS.get("dust_job_threads")

        sample_out = P.snip(self.outfile, self.fq1_suffix)
        
        fastq1 = self.fastq1
        fastq2 = self.fastq2
        fastq3 = self.fastq3

        outfile = self.outfile
        outfile1 = sample_out + '.1.fq.gz'
        outfile2 = sample_out + '.2.fq.gz'
        outfile3 = sample_out + '.3.fq.gz'
        out_disc1 = P.snip(self.outfile, '_masked' + self.fq1_suffix) \
            + '_discarded.fastq.1.fq.gz'
        out_disc2 = P.snip(self.outfile, '_masked' + self.fq1_suffix) \
            + '_discarded.fastq.2.fq.gz'
        out_disc3 = P.snip(self.outfile, '_masked' + self.fq1_suffix) \
            + '_discarded.fastq.3.fq.gz'
        
        if self.fastq2:
            if self.PARAMS['dust_discard_low_complexity']:
                statement1 = ("bbduk.sh"
                              "  in=%(fastq1)s"
                              "  in2=%(fastq2)s"
                              "  out=%(outfile1)s"
                              "  out2=%(outfile2)s"
                              "  outm=%(out_disc1)s"
                              "  outm2=%(out_disc2)s"
                              "  entropy=%(entropy)s"
                              "  threads=%(job_threads)s"
                              "  %(bb_options)s"
                              "  &> %(outfile)s.log" % locals())
                if IOTools.open_file(fastq3).read(1):
                    statement2 = (" bbduk.sh"
                                  "  in=%(fastq3)s"
                                  "  out=%(outfile3)s"
                                  "  outm=%(out_disc3)s"
                                  "  entropy=%(entropy)s"
                                  "  threads=%(job_threads)s"
                                  "  %(bb_options)s"
                                  "  &>> %(outfile)s.log" % locals())
                else:
                    statement2 = (" touch %(outfile3)s  %(out_disc3)s" % locals())

                statement = " && ".join([statement1, statement2])

            else:
                statement1 = ("bbmask.sh"
                              "  in=%(fastq1)s"
                              "  out=%(outfile1)s"
                              "  entropy=%(entropy)s"
                              "  threads=%(job_threads)s"
                              "  overwrite=t"
                              "  lowercase=t"
                              "  %(bb_options)s"
                              "  &> %(outfile)s.log &&"
                              " bbmask.sh"
                              "  in=%(fastq2)s"
                              "  out=%(outfile2)s"
                              "  entropy=%(entropy)s"
                              "  threads=%(job_threads)s"
                              "  overwrite=t"
                              "  lowercase=t"
                              "  %(bb_options)s"
                              "  &>> %(outfile)s.log" % locals())
                if IOTools.open_file(fastq3).read(1):           
                    statement2 = (" bbmask.sh"
                                  "  in=%(fastq3)s"
                                  "  out=%(outfile3)s"
                                  "  entropy=%(entropy)s"
                                  "  threads=%(job_threads)s"
                                  "  overwrite=t"
                                  "  lowercase=t"
                                  "  %(bb_options)s"
                                  "  &>> %(outfile)s.log" % locals())
                else:
                    statement2 = (" touch %(outfile3)s")

                statement = " && ".join([statement1, statement2])

        else:
            if self.PARAMS['dust_discard_low_complexity']:
                statement = ("bbduk.sh"
                             " in=%(fastq1)s"
                             " out=%(outfile1)s"
                             " outm=%(out_disc)s"
                             " entropy=%(entropy)s"
                             " threads=%(job_threads)s"
                             " lowercase=t"
                             " %(bb_options)s"
                             " &> %(outfile)s.log" % locals())

            else:
                statement = ("bbmask.sh"
                             " in=%(fastq1)s"
                             " out=%(outfile1)s"
                             " entropy=%(entropy)s"
                             " threads=%(job_threads)s"
                             " lowercase=t"
                             " %(bb_options)s"
                             " &> %(outfile)s.log" % locals())

        return statement
    
    def postProcess(self):
        sample_out = P.snip(self.outfile, self.fq1_suffix)
        fastq1 = self.fastq1
        fastq2 = self.fastq2
        fastq3 = sample_out + self.fq3_suffix
        outfile = self.outfile
        outfile1 = sample_out + '.1.fq.gz'
        outfile2 = sample_out + '.2.fq.gz'
        outfile3 = sample_out + '.3.fq.gz'
        out_disc1 = P.snip(self.outfile, '_masked' + self.fq1_suffix) \
            + '_discarded.fastq.1.fq.gz'
        out_disc2 = P.snip(self.outfile, '_masked' + self.fq1_suffix) \
            + '_discarded.fastq.2.fq.gz'
        out_disc3 = P.snip(self.outfile, '_masked' + self.fq1_suffix) \
            + '_discarded.fastq.3.fq.gz'
        if self.fastq2:
            # Renaming files because of bbmap idiosyncracies
            of1 = P.snip(outfile1, '.1.fq.gz') + self.fq1_suffix
            of2 = P.snip(outfile2, '.2.fq.gz') + self.fq2_suffix
            of3 = P.snip(outfile3, '.3.fq.gz') + self.fq3_suffix
            os.rename(outfile1, of1)
            os.rename(outfile2, of2)
            os.rename(outfile3, of3)

            if self.PARAMS['dust_discard_low_complexity']:
                od1 = P.snip(out_disc1, '.1.fq.gz') + self.fq1_suffix
                od2 = P.snip(out_disc2, '.2.fq.gz') + self.fq2_suffix
                od3 = P.snip(out_disc3, '.3.fq.gz') + self.fq3_suffix
                os.rename(out_disc1, od1)
                os.rename(out_disc2, od2)
                os.rename(out_disc3, od3)

        else:
            os.rename(outfile1, outfile)
            if self.PARAMS['dust_discard_low_complexity']:
                od1 = P.snip(out_disc1, '.1.fq.gz') + self.fq1_suffix
            os.rename(out_disc, od1)


def summariseReadCounts(infiles, outfile):
    with IOTools.open_file(outfile, 'w') as outf:
        outf.write("sample_id\tinput_reads\toutput_reads\tduplicates\t"
                   "adapter_contamination\trRNA\thost\tlow_complexity\t"
                   "duplicates_percent\tadapters_percent\trrna_percent\t"
                   "host_percent\tlow_complexity_perc\tremaining_percent\n")
        for infile in infiles:
            sample_id = P.snip(os.path.basename(infile),
                               '_read_count_summary.tsv')
            E.info('Processing sample %s' % sample_id)
            
            df = pd.read_table(infile, index_col=0, header=None)
            deadapt = df.loc['deadapt', 1]
            deduped = df.loc['deduped', 1]
            rrna = df.loc['rRNAremoved', 1]
            dehost = df.loc['dehost', 1]
            masked = df.loc['masked', 1]
            input_reads = df.loc['input', 1]
            
            lost_dup = int(input_reads) - int(deduped)
            lost_adapt = int(deduped) - int(deadapt)
            lost_rrna = int(deadapt) - int(rrna)
            lost_host = int(rrna) - int(dehost)
            lost_mask = int(dehost) - int(masked)

            lost_dup_perc = round(lost_dup/float(input_reads) * 100, 2)
            lost_adapt_perc = round(lost_adapt/float(input_reads) * 100, 2)
            lost_rrna_perc = round(lost_rrna/float(input_reads) * 100, 2)
            lost_host_perc = round(lost_host/float(input_reads) * 100, 2)
            lost_mask_perc = round(lost_mask/float(input_reads) * 100, 2)
            output_perc = round(float(masked)/float(input_reads) * 100, 2)

            outf.write('\t'.join(map(str, [sample_id, input_reads, masked, 
                                           lost_dup, lost_adapt, lost_rrna, 
                                           lost_host, lost_mask, lost_dup_perc, 
                                           lost_adapt_perc, lost_rrna_perc, 
                                           lost_host_perc, lost_mask_perc, 
                                           output_perc])) + '\n')

# Class for aligning reads with host genome with Hisat2
# returning mapped and unmapped reads
class hisat2(utility.matchReference):
    def __init__(self, fastq1, outfile, **PARAMS):
        super().__init__(fastq1, outfile, **PARAMS)
        # resetting prefix to include _unmapped or _dehost
        suffix = os.path.basename(self.outfile)
        self.prefixstrip = suffix[suffix.rfind('_'):]

    def hisatStatement(self):
        
        ref_genome = self.PARAMS["hisat2_ref_genome"]
    
        statement = [f"hisat2 -x {ref_genome}"]
        
        if self.fastq2 is None:
            unmapped = self.outfile.replace(self.prefixstrip, "_unmapped.fastq.gz")
            mapped = self.outfile.replace(self.prefixstrip, "_mapped.fastq.gz")
            entry = (f"-U {self.fastq1} --un-gz {unmapped}"
                     f" --al-gz {mapped}")
            statement.append(entry)
        elif self.fastq3 is None:
            unmapped = self.outfile.replace(self.prefixstrip, "_unmapped")
            mapped = self.outfile.replace(self.prefixstrip, "_mapped")
            unmapped_unpaired = self.outfile.replace(self.prefixstrip, 
                                                     "_unmapped.fastq.3.gz")
            mapped_unpaired = self.outfile.replace(self.prefixstrip,
                                                   "_mapped.fastq.3.gz")
            entry = (f"-1 {self.fastq1} -2 {self.fastq2}"
                     f" --un-conc-gz {unmapped} --un-gz {unmapped_unpaired}"
                     f" --al-conc-gz {mapped} --al-gz {mapped_unpaired}")
            statement.append(entry)
        elif self.fastq3 is not None:
            unmapped = self.outfile.replace(self.prefixstrip, "_unmapped")
            mapped = self.outfile.replace(self.prefixstrip, "_mapped")
            unmapped_unpaired = self.outfile.replace(self.prefixstrip, 
                                                     "_unmapped.fastq.3.gz")
            mapped_unpaired = self.outfile.replace(self.prefixstrip,
                                                   "_mapped.fastq.3.gz")
            entry = (f"-1 {self.fastq1} -2 {self.fastq2} -U {self.fastq3}"
                     f" --un-conc-gz {unmapped} --un-gz {unmapped_unpaired}"
                     f" --al-conc-gz {mapped} --al-gz {mapped_unpaired}")
            statement.append(entry)

        if self.fileformat == "fasta":
            statement.append("-f")
        elif self.fileformat == "fastq":
            statement.append("-q")

        samfile = self.outfile.replace(self.prefixstrip, ".sam")
        hisat2_summary = self.outfile.replace(self.prefixstrip, "_hisat2_summary.log")

        nthreads = self.PARAMS["hisat2_job_threads"]
        entry = (f"-S {samfile} --summary-file {hisat2_summary}"
                 " --new-summary"
                 " --met-stderr "
                 f" -p {nthreads} --reorder")
        statement.append(entry)

        statement.append(self.PARAMS["hisat2_options"])

        statement = " ".join(statement)
        return statement

    # method to sort and convert sam output to bam output
    def sam2bamStatement(self):
        samfile = self.outfile.replace(self.prefixstrip, ".sam")
        bamfile = re.sub("sam$", "bam", samfile)
        statement = f"samtools sort {samfile} > {bamfile}"
        
        return statement

    # clean up hisat2 outputs
    def postProcess(self):
        samfile = self.outfile.replace(self.prefixstrip, ".sam")
        
        # rename hisat output of paired end reads
        hisat_fq = {
            self.outfile.replace(self.prefixstrip, "_mapped.1"):
            self.outfile.replace(self.prefixstrip, "_mapped.fastq.1.gz"),
            self.outfile.replace(self.prefixstrip, "_mapped.2"):
            self.outfile.replace(self.prefixstrip, "_mapped.fastq.2.gz"),
            self.outfile.replace(self.prefixstrip, "_unmapped.1"):
            self.outfile.replace(self.prefixstrip, "_unmapped.fastq.1.gz"),
            self.outfile.replace(self.prefixstrip, "_unmapped.2"):
            self.outfile.replace(self.prefixstrip, "_unmapped.fastq.2.gz")
        }
        statements = []
        iter = [x for x in hisat_fq.keys() if os.path.exists(x)]
        for key in iter:
            statements.append(f"mv {key} {hisat_fq[key]}")

        # delete sam and files
        statements.append(f"rm {samfile}")
        
        statement = " && ".join(statements)
        
        return statement, hisat_fq
    
    # post-processing for pipeline_preprocess
    def postProcessPP(self):
        # rename hisat outputs to end in fastq.1.gz notation (if paired end)
        postprocess_statement, hisat_fq = self.postProcess()

        # add fastq3 to dictionary
        unmapped_fq3 = self.outfile.replace(self.prefixstrip,
                                            "_unmapped.fastq.3.gz")
        mapped_fq3 = self.outfile.replace(self.prefixstrip,
                                          "_mapped.fastq.3.gz")
        hisat_fq[unmapped_fq3] = unmapped_fq3
        hisat_fq[mapped_fq3] = mapped_fq3
        
        # check if files exist
        hisat_fq_found = {}
        for file in hisat_fq.keys():
            if os.path.exists(file):
                hisat_fq_found[file] = hisat_fq[file]
        
        # rename hisat output to pipline expected outfile
        new = [x.replace("unmapped","dehost") for x in hisat_fq_found.values()]
        new = [x.replace("mapped", "host") for x in new]
        rename = zip(hisat_fq.values(), new)
        rename_statements = [f"mv {x[0]} {x[1]}" for x in rename]
        
        statements = [postprocess_statement] + rename_statements
        statement = " && ".join(statements)

        return statement
    
    # method to read hisat summary file to dictionary
    def readHisatSummary(self, summary_file):
        summary_dict = {}    
        with open(summary_file, "r") as file:
            for line in file:
                line = line.strip()
                
                # Match key-value pairs with optional counts and percentages
                match = re.match(r"(.+?):\s+([\d.]+)(?:\s+\(([\d.]+)%\))?", line)
                if match:
                    key, count, percent = match.groups() # ignore percent
                    key = key.strip()
                    key = key.replace(" ", "_")
                    summary_dict[key] = count
        return summary_dict

    # method to merge sample hisat summaries into one table    
    def mergeHisatSummary(self, summary_files, merged_summary):
        summary = {}
        for summary_file in summary_files:
                entry = self.readHisatSummary(summary_file)
                summary[summary_file] = entry
        
        # merged summary file samples in rows summary in columns
        header = ["sample_name"] + list(summary[summary_file].keys())
        with open(merged_summary, "w") as f:
            f.write("\t".join(header) + "\n")
            for summary_file in summary_files:
                sample_name = os.path.basename(summary_file)
                sample_name = sample_name.strip("_hisat2_summary.log")
                entry = [sample_name] + list(summary[summary_file].values())
                entry = [str(x) for x in entry]
                f.write("\t".join(entry) + "\n")

    # wrapper method for running hisat in pipelines
    def hisat2bam(self):
        hisat_statement = self.hisatStatement()
        h2bam_statement = self.sam2bamStatement()
        statements = [hisat_statement, h2bam_statement]
        statement = " && ".join(statements)
        
        return(statement)
    
    # remove hisat fastq outputs
    def clean(self, outfile):
        # initialize log file
        open(outfile, "w").close()

        to_remove = glob(f"{self.outdir}/*.fastq*gz")
        statements = [f"rm -f {x} && echo 'clean up: deleted {x}' >> {outfile}" 
                      for x in to_remove]

        statement = " && ".join(statements)
        return statement
    
    # remove bam files when running hisat in pipeline_proeprocess
    def cleanPP(self, infiles, outfile):
        # initialize log file
        open(outfile, "w").close()

        # identify bam file
        to_remove = glob(f"{self.outdir}/*.bam")

        statements = []
        for x in to_remove:
            entry = (f"rm -f {x} "
                     f"&& echo 'clean up: deleted {x}' >> {outfile}")
            statements.append(entry)
        
        statement = " && ".join(statements)

        return statement
    