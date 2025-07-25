# Module for generic shotgun preprocessing steps

import os,re
import shutil
import pandas as pd
from glob import glob
from cgatcore import pipeline as P
from cgatcore import iotools as IOTools
from cgatcore import experiment as E

import ocmstoolkit.modules.Utility as Utility            
class Cdhit(Utility.BaseTool):
        
    def build_statement(self, fastn_obj):
        '''Filter exact duplicates, if specified in config file'''
        
        fastq1 = fastn_obj.fastn1
        fastq2 = fastn_obj.fastn2
        sample_out = P.snip(self.outfile, fastn_obj.fn1_suffix)
        
        outfile1 = P.snip(self.outfile, '.gz')
        logfile = sample_out + '.log'
        cluster_file = sample_out + '*.clstr'
        
        cdhit_options = self.PARAMS['cdhit_options']
        to_filter = self.PARAMS['cdhit_dedup']
        
        if fastn_obj.fastn2:
            outfile2 = re.sub(fastn_obj.fn1_suffix, fastn_obj.fn2_suffix, self.outfile)
            outfile2 = P.snip(outfile2, '.gz')
            
            if to_filter:
                tmpf1 = P.get_temp_filename('.')
                tmpf2 = P.get_temp_filename('.')
                statement = (f"zcat {fastn_obj.fastn1} > {tmpf1} &&"
                             f" zcat {fastq2} > {tmpf2} &&"
                             " cd-hit-dup"
                             f"  -i {tmpf1}"
                             f"  -i2 {tmpf2}"
                             f"  -o {outfile1}"
                             f"  -o2 {outfile2}"
                             f"  {cdhit_options}"
                             f" &> {logfile} &&"
                             f" gzip {outfile1} &&"
                             f" gzip {outfile2} &&"
                             f" gzip {logfile} &&"
                             f" rm -f {tmpf1} &&"
                             f" rm -f {tmpf2} &&"
                             f" rm -f {cluster_file}")
            else:
                E.warn('Deduplication step is being skipped for: %s' % fastq1)
                Utility.relink(fastq1, outfile1 + '.gz')
                Utility.relink(fastq2, outfile2 + '.gz')

        else:
            if to_filter:
                tmpf1 = P.get_temp_filename('.')
                statement = (f"zcat {fastq1} > {tmpf1}"
                             " cd-hit-dup"
                             f"  -i {tmpf1}"
                             f"  -o {outfile1}"
                             f"  {cdhit_options}"
                             f" &> {logfile} &&"
                             f" gzip {outfile1} &&"
                             f" gzip {logfile} &&"
                             f" rm -f {tmpf1} &&"
                             f" rm -f {cluster_file}")

            else:
                E.warn('Deduplication step is being skipped for: %s' % fastq1)
                Utility.relink(fastq1, self.outfile)     
        
        return statement
        
class Trimmomatic(Utility.BaseTool):

    def build_statement(self, fastn_obj):
        '''Remove adapters using Trimmomatic'''
        fastq1 = fastn_obj.fastn1
        fastq2 = fastn_obj.fastn2
        outfile1 = self.outfile
        sample_out = P.snip(self.outfile, fastn_obj.fn1_suffix)
        logfile = sample_out + '.trim.log'
        logfile2 = sample_out + '.log'
        
        trimmomatic_jar_path = self.PARAMS["trimmomatic_jar_path"]
        trimmomatic_n_threads = self.PARAMS["trimmomatic_job_threads"]
        
        # >0.32 determines phred format automatically
        # if phred score needs to be speficied, use trimmomatic_options
        
        trimmomatic_adapters = self.PARAMS["trimmomatic_adapters"]
        trimmomatic_seed_mismatches = self.PARAMS["trimmomatic_seed_mismatches"]
        trimmomatic_score_palendromic = self.PARAMS["trimmomatic_score_palendromic"]
        trimmomatic_score_simple = self.PARAMS["trimmomatic_score_simple"]
        trimmomatic_min_adapter_len = self.PARAMS["trimmomatic_min_adapter_len"]
        trimmomatic_keep_both_reads = self.PARAMS["trimmomatic_keep_both_reads"]
        trimmomatic_quality_leading = self.PARAMS["trimmomatic_quality_leading"]
        trimmomatic_quality_trailing = self.PARAMS["trimmomatic_quality_trailing"]
        trimmomatic_minlen = self.PARAMS["trimmomatic_minlen"]
        trimmomatic_options = self.PARAMS["trimmomatic_options"]

        assert "-threads" not in trimmomatic_options, (
            "Trimmomatic multi-threading is set with job_threads parameter."
            )

        if fastn_obj.fastn2:
            outfile2 = re.sub(fastn_obj.fn1_suffix, fastn_obj.fn2_suffix, self.outfile)
            outf1_singletons = sample_out + re.sub("1", "1s", fastn_obj.fn1_suffix)
            outf2_singletons = sample_out + re.sub("2", "2s", fastn_obj.fn2_suffix)
            outf_singletons = sample_out + fastn_obj.fn3_suffix
            
            statement = (f"java -Xmx5g -jar {trimmomatic_jar_path} PE"
                         f" -threads {trimmomatic_n_threads}"
                         f" -trimlog {logfile}"
                         f" {fastq1}" # input read 1
                         f" {fastq2}" # input read 2
                         f" {outfile1}" # output read 1
                         f" {outf1_singletons}" # output unpaired read 1
                         f" {outfile2}" # output read 2
                         f" {outf2_singletons}" # output unpaired read 2
                         " ILLUMINACLIP:"
                         f"{trimmomatic_adapters}:"
                         f"{trimmomatic_seed_mismatches}:"
                         f"{trimmomatic_score_palendromic}:"
                         f"{trimmomatic_score_simple}:"
                         f"{trimmomatic_min_adapter_len}:"
                         f"{trimmomatic_keep_both_reads}"
                         f" LEADING:{trimmomatic_quality_leading}"
                         f" TRAILING:{trimmomatic_quality_trailing}"
                         f" MINLEN:{trimmomatic_minlen}"
                         f" {trimmomatic_options}"
                         f" &> {logfile2} &&"
                         f" gzip -f {logfile} &&"
                         f" cat {outf1_singletons} {outf2_singletons} "
                         f"  > {outf_singletons} &&"
                         f" rm -f {outf1_singletons} && rm -f {outf2_singletons}")

        else:
            statement = (f"java -Xmx5g -jar {trimmomatic_jar_path} SE"
                         f" -threads {trimmomatic_n_threads}"
                         f" -trimlog {logfile}"
                         f" {fastq1}" # input read 1
                         f" {outfile1}" # output read 1
                         " ILLUMINACLIP:"
                         f"{trimmomatic_adapters}:"
                         f"{trimmomatic_seed_mismatches}:"
                         f"{trimmomatic_score_palendromic}:"
                         f"{trimmomatic_score_simple}"
                         f"{trimmomatic_min_adapter_len}:"
                         f"{trimmomatic_keep_both_reads}"
                         f" LEADING:{trimmomatic_quality_leading}"
                         f" TRAILING:{trimmomatic_quality_trailing}"
                         f" MINLEN:{trimmomatic_minlen}"
                         f" {trimmomatic_options}"
                         f" &> {logfile2} &&"
                         f" gzip -f {logfile}")
                
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

    
class SortMeRNA(Utility.BaseTool):
    """
    Run sortMeRNA. 
    Assumes that reference indexes have been created in advance in a 
    specified location.
    """

    def __init__(self, fastq1, outfile, **PARAMS):
        super().__init__(fastq1, outfile, **PARAMS)
        self.sortmerna_skip_singletons = self.PARAMS.get('sortmerna_skip_singletons', False)
        
    def build_statement(self, fastn_obj):
        """
        Generate run statement for processing single, paired, or paired
        + singleton samples. 

        Required arguments: 
        index
        reference
        
        """
        sortmerna_options = self.PARAMS.get("sortmerna_options")
        assert "--threads" not in sortmerna_options, (
            "Sortmerna multi-threading is set with job_threads parameter."
            )
        
        # A comma separated list of references
        references = self.PARAMS.get("sortmerna_reference")
        references = ' --ref '.join(references.split(','))
        # All listed references must be pre-indexed in this location
        index_dir = self.PARAMS.get("sortmerna_index") # Check this isn't automatically passed. 
        
        tmpf = P.get_temp_dir('.')
        tmpf_kvdb = os.path.join(tmpf, 'kvdb')
        tmpf_readb = os.path.join(tmpf, 'readb')

        job_threads = self.PARAMS.get("sortmerna_job_threads")
        if not fastn_obj.fastn2:
            # Run sortMeRNA for single reads
            in_fastq1 = fastn_obj.fastn1
            in_prefix = P.snip(in_fastq1, '_deadapt'+fastn_obj.fn1_suffix, strip_path=True)
            out_prefix = os.path.join(self.outdir, in_prefix)

            # Run sortMeRNA for single reads
            statement = ("sortmerna"
                         #" --index 0" # skip indexing, assume in idx-dir
                         " --fastx"
                         f" --reads {in_fastq1}"
                         f" --ref {references}"
                         f" {index_dir}" # location of reference indexes
                         f" --aligned {out_prefix}_aligned" # output location of aligned seq
                         f" --other {out_prefix}_unaligned" # output location of unalinged seq
                         f" --readb {tmpf_readb}" # location of tmp file for reads
                         f" --kvdb {tmpf_kvdb}" # location of tmp file for kv pairs
                         f" --threads {job_threads}"
                         " --zip-out"
                         f" {sortmerna_options}")

        else:
            # Run sortMeRNA for paired reads
            in_fastq1 = fastn_obj.fastn1
            in_fastq2 = fastn_obj.fastn2
            in_prefix = P.snip(in_fastq1, '_deadapt'+fastn_obj.fn1_suffix, strip_path=True)
            out_prefix = os.path.join(self.outdir, in_prefix)
            # Run sortMeRNA for single reads
            statement = ("sortmerna"
                         " --index 0" # skip indexing, assume in idx-dir
                         " --fastx"
                         f" --reads {in_fastq1}" # First read file
                         f" --reads {in_fastq2}" # Second read file
                         f" --ref {references}"
                         f" --idx-dir {index_dir}" # location of reference indexes
                         f" --aligned {out_prefix}_aligned" # output location of aligned seq
                         f" --other {out_prefix}_unaligned" # output location of unalinged seq
                         f" --readb {tmpf_readb}" # location of tmp file for reads
                         f" --kvdb {tmpf_kvdb}" # location of tmp file for kv pairs
                         " --paired_in" # If one read is aligned, both are output to aligned file
                         " --out2" # Output paired reads to separate files
                         f" --threads {job_threads}"
                         " --zip-out"
                         f" {sortmerna_options}")
        
        if not self.sortmerna_skip_singletons and IOTools.open_file(fastn_obj.fastn3).read(1):
            in_fastq3 = fastn_obj.fastn3
            statement_2 = ("sortmerna"
                           # " --index 0" # skip indexing, assume in idx-dir
                           " --fastx"
                           f" --reads {in_fastq3}"
                           f" --idx-dir {index_dir}" # location of reference indexes
                           f" --ref {references}"
                           f" --aligned {out_prefix}_aligned_singleton" 
                           f" --other  {out_prefix}_unaligned_singleton"
                           f" --readb {tmpf_readb}" # location of tmp file for reads
                           f" --kvdb {tmpf_kvdb}" # location of tmp file for kv pairs
                           f" --threads {job_threads}"
                           " --zip-out 1"
                           f" {sortmerna_options}") 
            
            statement = " && ".join([statement, 
                                     f"rm -rf {tmpf}/*", # location of tmp_readb & kvdb
                                     statement_2,
                                     f"rm -rf {tmpf}"])
        else:
            statement = " && ".join([statement, 
                                     f"rm -rf {tmpf}"])

        return statement
        
    def post_process(self, fastn_obj):
        ''' Rename files output by sortmeRNA to appropriate suffix
        At some point this would be good to become more flexible wrt FQ1_SUFFIX'''

        outf_prefix = os.path.join(self.outdir, 
                                   P.snip(fastn_obj.fastn1, 
                                          '_deadapt'+fastn_obj.fn1_suffix, 
                                          strip_path=True))

        # rename fastq1 files
        os.rename(outf_prefix + '_aligned_fwd.fq.gz', 
                  outf_prefix + '_rRNA'+fastn_obj.fn1_suffix)
        os.rename(outf_prefix + '_unaligned_fwd.fq.gz', 
                  outf_prefix + '_rRNAremoved'+fastn_obj.fn1_suffix)

        # rename fastq2 files
        if fastn_obj.fastn2:
            os.rename(outf_prefix + '_aligned_rev.fq.gz', 
                      outf_prefix + '_rRNA'+fastn_obj.fn2_suffix)
            os.rename(outf_prefix + '_unaligned_rev.fq.gz', 
                      outf_prefix + '_rRNAremoved'+fastn_obj.fn2_suffix)

        # rename fastq3 files
        if not self.sortmerna_skip_singletons and IOTools.open_file(fastn_obj.fastn3).read(1):            
            # for some reason, singletons are in fasta instead of fastq
            f3_aligned = glob(os.path.join(self.outdir, "*aligned_singleton*"))[0]

            os.rename(f3_aligned,
                      outf_prefix + '_rRNA' + fastn_obj.fn3_suffix)

            f3_unaligned = glob(os.path.join(self.outdir, "*unaligned_singleton*"))[0]            
            os.rename(f3_unaligned,
                      outf_prefix +  '_rRNAremoved' + fastn_obj.fn3_suffix)
        else:
            Utility.relink(fastn_obj.fastn3, 
                           outf_prefix + "_rRNAremoved" + fastn_obj.fn3_suffix)
        return None


class createSortMeRNAOTUs(SortMeRNA):
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
    
    def post_process(self, fastn_obj):
        '''Rename files output by sortmerna, including otu_map table.''' 

        tmpf_prefix = os.path.join(self.outdir, 
                                   P.snip(fastn_obj.fastn1, '_deadapt'+fastn_obj.fn1_suffix, strip_path=True))
        outf_prefix = os.path.join(os.path.dirname(self.outdir),
                                   P.snip(fastn_obj.fastn1, '_deadapt'+fastn_obj.fn1_suffix, strip_path=True))

        # rename fastq1 files
        os.rename(tmpf_prefix + '_aligned_fwd.fq.gz', outf_prefix + '_rRNA'+fastn_obj.fn1_suffix)
        os.rename(tmpf_prefix + '_unaligned_fwd.fq.gz', outf_prefix + '_rRNAremoved'+fastn_obj.fn1_suffix)

        # rename fastq2 files
        if fastn_obj.fastn2:
            os.rename(tmpf_prefix + '_aligned_rev.fq.gz', outf_prefix + '_rRNA'+fastn_obj.fn2_suffix)
            os.rename(tmpf_prefix + '_unaligned_rev.fq.gz', outf_prefix + '_rRNAremoved'+fastn_obj.fn2_suffix)

        # rename 'denovo' otus
        if re.search('de_novo_otu', self.sortmerna_options):
            os.rename(tmpf_prefix + '_aligned_denovo_fwd.fq.gz', outf_prefix + '_denovo'+fastn_obj.fn1_suffix)
            os.rename(tmpf_prefix + '_aligned_denovo_rev.fq.gz', outf_prefix + '_denovo'+fastn_obj.fn2_suffix)

        # rename log file 
        os.rename(tmpf_prefix + '_aligned.log', outf_prefix + '.log')

        # rename otu_map.txt
        otu_file = os.path.join(self.outdir, 'otu_map.txt')
        os.rename(otu_file, outf_prefix + '_otu_map.txt')

        shutil.rmtree(self.outdir)
        
        return None

class Bmtagger(Utility.BaseTool):

    def build_statement(self, fastn_obj):
        '''Remove host contamination using bmtagger'''
        outf_host = P.snip(self.outfile, '_dehost'+fastn_obj.fn1_suffix) + '_host.txt'
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

        fastq1 = fastn_obj.fastn1
        outfile = self.outfile
        bmtagger_exec = self.PARAMS['bmtagger_executable']
        assert bmtagger_exec in ["bmtagger.sh", "bmtagger_mod.sh"], (
            "must specify bmtagger.sh or bmtagger_mod.sh"
        )

        if fastn_obj.fastn2:
            fastq2 = fastn_obj.fastn2
            fastq3 = fastn_obj.fastn3
            
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
                statement1 = (f"zcat {fastq1} > {tmpf1} &&"
                              f" zcat {fastq2} > {tmpf2} &&"
                              f" {bmtagger_exec}"
                              f"  -b {bitmask}"
                              f"  -x {srprism}"
                              f"  -T {tmpdir1}"
                              "  -q1" # Input is fastq
                              f"  -1 {tmpf1}"
                              f"  -2 {tmpf2}"
                              f"  -o {outf_host_stub}_paired{n}"
                              f"  &> {outfile}.log &&"
                              f" cat {outf_host_stub}_paired{n}"
                              f"  >> {to_remove_paired} &&"
                              f" rm -rf {tmpdir1} {tmpf1} {tmpf2}"
                              f"  {outf_host_stub}_paired{n}")

                # Screen the singletons
                if os.path.exists(fastn_obj.fastn3) and IOTools.open_file(fastn_obj.fastn3).read(1):
                    statement2 = (f"zcat {fastq3} > {tmpf3} &&"
                                  f" {bmtagger_exec}"
                                  f"  -b {bitmask}"
                                  f"  -x {srprism}"
                                  f"  -T {tmpdir2}"
                                  "  -q1" # Input is fastq
                                  f"  -1 {tmpf3}"
                                  f"  -o {outf_host_stub}_singletons{n}"
                                  f" &>> {outfile}.log &&"
                                  f" cat {outf_host_stub}_singletons{n}"
                                  f"  >> {to_remove_singletons} &&"
                                  f" rm -rf {tmpdir2} {tmpf3}"
                                  f"  {outf_host_stub}_singletons{n}")
                else:
                    statement2 = (f"touch  {to_remove_singletons} &&"
                                  f" rm -rf {tmpdir2} {tmpf3}")

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
                
                statement = (f"zcat {fastq1} > {tmpf} &&"
                             f" {bmtagger_exec}"
                             f"  -b {bitmask}"
                             f"  -x {srprism}"
                             f"  -T {tmpdir1}"
                             "  -q1" # Input is fastq
                             f"  -1 {tmpf}"
                             f"  -o {outf_host_stub}_{n}"
                             f"  &>> {outfile}.log &&"
                             f" cat {outf_host_stub}_{n} >> {to_remove}"
                             f" rm -rf {tmpdir1} {tmpf} {outf_host_stub}_{n}")
                statements.append(statement)
                
        return statements, to_remove_tmp

    def post_process(self, fastn_obj, to_remove_tmp):

        if fastn_obj.fastn2:
            # Drop host contaminated reads
            # A hack due to the fact that BMTagger truncates fastq identifiers
            # TO DO: Look at bmtagger/.../bin/extract_fullseq
            
            fastq1 = fastn_obj.fastn1
            fastq2 = fastn_obj.fastn2
            
            fastq1_out = self.outfile
            fastq2_out = P.snip(self.outfile, fastn_obj.fn1_suffix) + fastn_obj.fn2_suffix
            
            fastq1_host = P.snip(self.outfile, '_dehost'+fastn_obj.fn1_suffix) + '_host'+fastn_obj.fn1_suffix
            fastq2_host = P.snip(self.outfile, '_dehost'+fastn_obj.fn1_suffix) + '_host'+fastn_obj.fn2_suffix
            
            fastq3 = fastn_obj.fastn3
            fastq3_out = P.snip(self.outfile, fastn_obj.fn1_suffix) + fastn_obj.fn3_suffix
            fastq3_host = P.snip(self.outfile, '_dehost'+fastn_obj.fn1_suffix) + '_host'+fastn_obj.fn3_suffix
            to_remove_paired = to_remove_tmp[0]
            to_remove_singletons = to_remove_tmp[1]

            statement = ("ocms_shotgun drop_fastqs"
                         f" --fastq1 {fastq1}"
                         f" --fastq2 {fastq2}"
                         f" --fastq3 {fastq3}"
                         f" --to-drop-paired {to_remove_paired}"
                         f" --to-drop-single {to_remove_singletons}"
                         f" --fastq-out1 {fastq1_out}"
                         f" --fastq-out2 {fastq2_out}"
                         f" --fastq-out3 {fastq3_out}"
                         f" --fastq-drop1 {fastq1_host}"
                         f" --fastq-drop2 {fastq2_host}"
                         f" --fastq-drop3 {fastq3_host}"
                         f" &>> {fastq1_out}.log")

            to_unlink = [to_remove_paired, to_remove_singletons]

        else:
            
            fastq1 = fastn_obj.fastn1
            outfile = self.outfile
            to_remove = to_remove_tmp[0]
            fastq_host = P.snip(self.outfile, '_dehost'+fastn_obj.fn1_suffix) + '_host'+fastn_obj.fn1_suffix
            statement = ("ocms_shotgun drop_fastqs"
                         f" --fastq1 {fastq1}"
                         f" --to-drop-single {to_remove}"
                         f" --fastq-out1 {outfile}"
                         f" --fastq-drop1 {fastq_host}"
                         f" &>> {outfile}.log")
        
            P.run(statement)
            
            os.unlink(to_remove)
            to_unlink = [to_remove]

        return statement, to_unlink


# Class for aligning reads with host genome with Hisat2
# returning mapped and unmapped reads
class Hisat2(Utility.BaseTool):
    def __init__(self, fastq1, outfile, **PARAMS):
        super().__init__(fastq1, outfile, **PARAMS)
        # resetting prefix to include _unmapped or _dehost
        suffix = os.path.basename(self.outfile)
        self.prefixstrip = suffix[suffix.rfind('_'):]
    def hisat_statement(self, fastn_obj):
        
        ref_genome = self.PARAMS["hisat2_ref_genome"]
    
        # job threads set with job_options so want to make sure threads is not
        # seperately defined
        banned_options = ["-p","--threads"]
        check_options = [x not in self.PARAMS["hisat2_options"] for x in banned_options]
        assert any(check_options), (
            "Hisat2 multi-threading is set with job_options"
        )
        statement = [f"hisat2 -x {ref_genome}"]
        
        # single end reads
        if fastn_obj.fastn2 is None:
            unmapped = self.outfile.replace(self.prefixstrip, "_unmapped.fastq.gz")
            mapped = self.outfile.replace(self.prefixstrip, "_mapped.fastq.gz")
            entry = (f"-U {fastn_obj.fastn1} --un-gz {unmapped}"
                     f" --al-gz {mapped}")
            statement.append(entry)
        # paired end reads, no singletons
        elif not fastn_obj.has_singleton:
            unmapped = self.outfile.replace(self.prefixstrip, "_unmapped")
            mapped = self.outfile.replace(self.prefixstrip, "_mapped")
            unmapped_unpaired = self.outfile.replace(self.prefixstrip, 
                                                     "_unmapped.fastq.3.gz")
            mapped_unpaired = self.outfile.replace(self.prefixstrip,
                                                   "_mapped.fastq.3.gz")
            entry = (f"-1 {fastn_obj.fastn1} -2 {fastn_obj.fastn2}"
                     f" --un-conc-gz {unmapped} --un-gz {unmapped_unpaired}"
                     f" --al-conc-gz {mapped} --al-gz {mapped_unpaired}")
            statement.append(entry)
        # paired end reads with singletons
        elif fastn_obj.has_singleton:
            unmapped = self.outfile.replace(self.prefixstrip, "_unmapped")
            mapped = self.outfile.replace(self.prefixstrip, "_mapped")
            unmapped_unpaired = self.outfile.replace(self.prefixstrip, 
                                                     "_unmapped.fastq.3.gz")
            mapped_unpaired = self.outfile.replace(self.prefixstrip,
                                                   "_mapped.fastq.3.gz")
            entry = (f"-1 {fastn_obj.fastn1} -2 {fastn_obj.fastn2} -U {fastn_obj.fastn3}"
                     f" --un-conc-gz {unmapped} --un-gz {unmapped_unpaired}"
                     f" --al-conc-gz {mapped} --al-gz {mapped_unpaired}")
            statement.append(entry)

        if fastn_obj.fileformat == "fasta":
            statement.append("-f")
        elif fastn_obj.fileformat == "fastq":
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
    def sam2bam_statement(self):
        samfile = self.outfile.replace(self.prefixstrip, ".sam")
        bamfile = re.sub("sam$", "bam", samfile)
        statement = f"samtools sort {samfile} > {bamfile}"
        
        return statement

    # clean up hisat2 outputs
    def post_process(self):
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
    def post_process_pp(self):
        # rename hisat outputs to end in fastq.1.gz notation (if paired end)
        postprocess_statement, hisat_fq = self.post_process()

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
    def read_hisat_summary(self, summary_file):
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
    def merge_hisat_summary(self, summary_files, merged_summary):
        summary = {}
        for summary_file in summary_files:
                entry = self.read_hisat_summary(summary_file)
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
    def hisat2bam(self, fastn_obj):
        hisat_statement = self.hisat_statement(fastn_obj)
        h2bam_statement = self.sam2bam_statement()
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
    def clean_pp(self, outfile):
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
    
class Bbtools(Utility.BaseTool):

    def build_statement(self, fastn_obj):
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
        bb_options= self.PARAMS['dust_options']
        job_threads = self.PARAMS.get("dust_job_threads")

        # make sure threads is not set outside of job_threads
        assert "threads" not in bb_options, (
            "BBtools multi-threading is set with jop_options"
        )

        sample_out = P.snip(self.outfile, fastn_obj.fn1_suffix)
        
        fastq1 = fastn_obj.fastn1
        fastq2 = fastn_obj.fastn2
        fastq3 = fastn_obj.fastn3

        outfile = self.outfile
        outfile1 = sample_out + '.1.fq.gz'
        outfile2 = sample_out + '.2.fq.gz'
        outfile3 = sample_out + '.3.fq.gz'
        out_disc1 = outfile1.replace('_masked', '_discarded')
        out_disc2 = out_disc1.replace('.1.fq.gz', '.2.fq.gz') 
        out_disc3 = out_disc1.replace('.1.fq.gz', '.3.fq.gz')
        
        if fastn_obj.fastn2:
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
    
    def post_process(self, fastn_obj):
        sample_out = P.snip(self.outfile, fastn_obj.fn1_suffix)
        fastq1 = fastn_obj.fastn1
        fastq2 = fastn_obj.fastn2
        fastq3 = fastn_obj.fn3_suffix
        outfile = self.outfile
        outfile1 = sample_out + '.1.fq.gz'
        outfile2 = sample_out + '.2.fq.gz'
        outfile3 = sample_out + '.3.fq.gz'
        out_disc1 = outfile1.replace('_masked', '_discarded')
        out_disc2 = out_disc1.replace('.1.fq.gz', '.2.fq.gz') 
        out_disc3 = out_disc1.replace('.1.fq.gz', '.3.fq.gz')

        if fastn_obj.fastn2:
            # Renaming files because of bbmap idiosyncracies
            of1 = P.snip(outfile1, '.1.fq.gz') + fastn_obj.fn1_suffix
            of2 = P.snip(outfile2, '.2.fq.gz') + fastn_obj.fn2_suffix
            of3 = P.snip(outfile3, '.3.fq.gz') + fastn_obj.fn3_suffix
            os.rename(outfile1, of1)
            os.rename(outfile2, of2)
            os.rename(outfile3, of3)

            if self.PARAMS['dust_discard_low_complexity']:
                od1 = P.snip(out_disc1, '.1.fq.gz') + fastn_obj.fn1_suffix
                od2 = P.snip(out_disc2, '.2.fq.gz') + fastn_obj.fn2_suffix
                od3 = P.snip(out_disc3, '.3.fq.gz') + fastn_obj.fn3_suffix
                os.rename(out_disc1, od1)
                os.rename(out_disc2, od2)
                os.rename(out_disc3, od3)

        else:
            os.rename(outfile1, outfile)
            if self.PARAMS['dust_discard_low_complexity']:
                od1 = P.snip(out_disc1, '.1.fq.gz') + fastn_obj.fn1_suffix
            os.rename(out_disc1, od1)


def summarise_read_counts(infiles, outfile):
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
