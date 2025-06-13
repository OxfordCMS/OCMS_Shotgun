##############################################################################
#

#   $Id$
#
#   Copyright (C) 2017 Jethro Johnson
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################

"""
==========================
PipelineMetagenomeAssembly
==========================

Classes and utility functions for metagenomics assembly

Metagenome assembly tools
-------------------------
SPAdes Meta


"""

import os, errno
import re
import yaml
import pickle
import glob
import collections
import random
import shutil
import numpy as np
import distutils.spawn as spawn

import cgatcore.iotools as IOTools
import cgatcore.experiment as E
import cgat.FastaIterator as FastaIterator
import cgat.Fastq as Fastq

import cgatcore.pipeline as P

###############################################################################
# Base classes

class FetchData(object):

    '''
    class for assessing and retrieving data files
    '''

    def __init__(self):

        self.format = None
        self.paired = False
        self.singletons = False

    def getTrack(self, infile):
        '''
        return the track for the file
        '''
        return P.snip(os.path.basename(infile),
                      ".%s" % self.getFormat(infile),
                      strip_path=True)

    def getFormat(self, infile):
        possible_formats = ["fasta", "fastq", "fasta.gz",
                            "fastq.gz", "fasta.1.gz", "fastq.1.gz"]

        for f in possible_formats:
            if infile.endswith(f):
                self.format = f
        assert self.format, ("File %s is not of correct format. Needs to be"
                        " one of: %s" % (infile, ' '.join(possible_formats)))
        
        return self.format

        
        self.format = self.fileFormat(infile)
        return self.format

    def getPairedFiles(self, infile, check_exists=True):
        '''
        check whether paired files (.2.) and singleton files (.3.)
        exist for the appropriate file extensions, return a list of
        input files.
        '''
        file_list = [infile,]
        
        format = self.getFormat(infile)
        
        if format in ['fastq.1', 'fasta.1',
                      'fastq.1.gz', 'fasta.1.gz']:
            read2 = P.snip(infile, format) + re.sub('1', '2', format)
            read3 = P.snip(infile, format) + re.sub('1', '3', format)
           
            if check_exists:
                if os.path.exists(os.path.abspath(read2)):
                    self.paired = True
                    file_list.append(read2)
                    if os.path.exists(os.path.abspath(read3)):
                        self.singletons = True
                        file_list.append(read3)
                elif os.path.exists(os.path.abspath(read3)):
                    raise IOError("Singleton file exists %s, but paired file (%s)"
                    "not found" % (read3, read2))
                else:
                    pass
            else:
                file_list.extend([read2, read3])
                
        return file_list

    
class MetaAssembler(FetchData):
    '''Generic assembler class. '''
    
    def preProcess(self, infiles, outfile, **PARAMS):
        '''
        detect files containing reads, determine whether this
        is a hybrid assembly
        '''
        return ""

    def assembler(self, infiles, outfile, **PARAMS):
        '''
        Receive a list of up to three infiles, representing read1, read2
        and singletons, generate a run command for a specific assembler.
        '''
        return ""
        
    def postProcess(self, infiles, outfile, **PARAMS):
        '''
        return assembled contigs in a generic format
        '''
        return ""

    def build(self, infile, outfile, **PARAMS):
        '''Create run statement for assembler'''
        
        infiles = self.getPairedFiles(infile)
        cmds = []
        cmd_preprocess = self.preProcess(infiles, outfile, **PARAMS)
        cmd_assembly = self.assembler(infiles, outfile, **PARAMS)
        cmd_postprocess = self.postProcess(infiles, outfile, **PARAMS)

        cmds = [x for x in (cmd_preprocess, cmd_assembly, cmd_postprocess) if x]
        statement = ' && '.join(cmds)

        return statement

###############################################################################
# Run SPAdes

class SpadesReadCorrection(MetaAssembler):
    '''Run BayesHammer via Spades, dynamically compute --threads and --memory
        from spades_ec_threads and spades_ec_memory in PARAMS.'''

    def fetch_run_statement(self, libraries, out_sub_dir, **PARAMS):
        '''Generate commandline statement for running BayesHammer'''
        
        ec_threads = PARAMS['spades_ec_threads']  # e.g., 4
        ec_memory_per_thread = PARAMS['spades_ec_memory']  # e.g., '20G'

        # Remove 'G' and calculate total memory in GB
        memory_value = int(ec_memory_per_thread.strip().upper().replace('G', ''))
        total_memory = ec_threads * memory_value

        # Build run options string
        run_options = f"--threads {ec_threads} --memory {total_memory}"
                     
        statement = ("spades.py"
                     " --only-error-correction"
                     " {libraries}"
                     " -o {out_sub_dir}"
                     " {run_options}".format(**locals()))

        return statement
        
    def assembler(self, infiles, outfile, **PARAMS):
        ''' '''
        # set up sub-directories for spades run output 
        out_dir = os.path.dirname(outfile)
        out_sub_dir = os.path.join(out_dir, self.getTrack(infiles[0]))
        if not os.path.exists(out_sub_dir):
            os.mkdir(out_sub_dir)

        # spades expects either fastq, fasta, fastq.gz, fasta.gz suffix
        sym_files = [os.path.basename(x) for x in infiles]
        sym_files = [re.sub('.fast(.).(.).gz', r'.\2.fast\1.gz', x) for x in sym_files]
        sym_files = [os.path.join(out_dir, x) for x in sym_files]
        for i, infile in enumerate(infiles):
            try:
                os.symlink(os.path.abspath(infile), sym_files[i])
            except OSError as e:
                if e.errno == errno.EEXIST:
                    os.remove(sym_files[i])
                    os.symlink(os.path.abspath(infile), sym_files[i])
                else:
                    raise e

        # specify single or p.e. libraries
        if len(infiles) == 1:
            assert self.paired == False
            libraries = '--s ' + infiles[0]
        elif len(infiles) == 3:
            assert self.paired == True and self.singletons == True
            libraries = zip(['--pe1-1', '--pe1-2', '--pe1-s'], sym_files)
            libraries = ' '.join([' '.join(x) for x in libraries])
        else:
            assert self.paired == True and self.singletons == False
            libraries = zip(['--pe1-1', '--pe1-2'], sym_files)
            libraries = ' '.join([' '.join(x) for x in libraries])

        # specify run command
        statement = self.fetch_run_statement(libraries, out_sub_dir, **PARAMS)

        return statement

    def postProcess(self, infiles, outfile, **PARAMS):
        '''Remove the renamed symlinks used as spades input'''
        out_dir = os.path.dirname(outfile)
        sym_files = [os.path.basename(x) for x in infiles]
        sym_files = [re.sub('.fast(.).(.).gz', r'.\2.fast\1.gz', x) for x in sym_files]
        sym_files = [os.path.join(out_dir, x) for x in sym_files]

        statement = 'rm -f %s' % ' '.join(sym_files)

        return statement

    
class fetchSpadesProcessedReads(SpadesReadCorrection):
    '''Using BayesHammer outputs a yaml file containing names of cleaned
    read files. Rather than predict what they're called, read yaml.
    '''
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, infile, outfile):
        infiles = self.getPairedFiles(infile)
        outfiles = self.getPairedFiles(outfile, check_exists=False)
        
        out_yaml = os.path.join(os.path.dirname(outfile),
                                self.getTrack(infiles[0]),
                                'corrected',
                                'corrected.yaml')
        cleaned_dict = yaml.load(open(out_yaml), Loader=yaml.Loader)[0]

        statement = ''
        if len(infiles) == 1:
            infile = ' '.join(cleaned_dict['single reads'])
            statement += 'cat %s > %s' % (infile, outfiles[0])
        elif len(infiles) == 2:
            infile1 = ' '.join(cleaned_dict['left reads'])
            statement += 'cat %s > %s' % (infile1, outfiles[0])
            infile2 = ' '.join(cleaned_dict['right reads'])
            statement += '; cat %s > %s' % (infile2, outfiles[1])
        elif len(infiles) == 3:
            infile1 = ' '.join(cleaned_dict['left reads'])
            statement += 'cat %s > %s' % (infile1, outfiles[0])
            infile2 = ' '.join(cleaned_dict['right reads'])
            statement += '; cat %s > %s' % (infile2, outfiles[1])
            infile3 = ' '.join(cleaned_dict['single reads'])
            statement += '; cat %s > %s' % (infile3, outfiles[2])
        else:
            raise IOError("Unexpected number of files specified when merging"
                          " spades output. Script currently only takes 3:"
                          " %s" % ' '.join(infiles))
        
        return statement

    
class runMetaSpades(SpadesReadCorrection):
    '''Run spades --meta'''

    def fetch_run_statement(self, libraries, out_sub_dir, **PARAMS):
        '''Generate commandline statement for running metaSPAdes without
        BayesHammer'''
        
        meta_threads = PARAMS['spades_meta_threads']  
        meta_memory_per_thread = PARAMS['spades_meta_memory']  # e.g., '30G'

        # Remove 'G' and calculate total memory in GB
        memory_value = int(meta_memory_per_thread.strip().upper().replace('G', ''))
        total_memory = meta_threads * memory_value

        # Construct run options dynamically
        run_options = f"--threads {meta_threads} --memory {total_memory}"
        
        statement1 = ("spades.py"
                     " --meta"
                     " --only-assembler"
                     " {libraries}"
                     " -o {out_sub_dir}"
                     " {run_options}".format(**locals()))

        # fetch scaffolds and contigs
        out_dir = os.path.dirname(out_sub_dir)
        out_prefix = os.path.basename(out_sub_dir)

        contigs = os.path.join('..', out_sub_dir, 'contigs.fasta')
        scaffolds = os.path.join('..', out_sub_dir, 'scaffolds.fasta')

        cont_out = os.path.join(out_dir,
                                out_prefix + '.spades.contigs.fasta')
        scaf_out = os.path.join(out_dir,
                                out_prefix + '.spades.scaffolds.fasta')

        statement2 = (" ln -s {contigs} {cont_out};"
                      " ln -s {scaffolds} {scaf_out}".format(**locals()))

        statements = ' && '.join([statement1, statement2])
        
        return statements

    

###############################################################################
# Run MEGAHIT
class runMegaHit(MetaAssembler):
    '''Run MEGAHIT assembler'''
    
    
    def assembler(self, infiles, outfile, **PARAMS):
        ''' '''
        # set up sub-directories for run output
        out_dir = os.path.dirname(outfile)
        out_sub_dir = os.path.join(out_dir, self.getTrack(infiles[0]))
        # if not os.path.exists(out_sub_dir):
        #     os.mkdir(out_sub_dir)

        if len(infiles) == 1:
            # check whether input files are interlaced
            if self.checkPairs(infiles[0]):
                libraries =  '--12 ' + infiles[0]
            else:
                libraries = '-r ' + infiles[0]
        else:
            libraries = zip(['-1', '-2', '-r'], infiles)
            libraries = ' '.join([' '.join(x) for x in libraries])


        #mega_threads = PARAMS['megahit_meta_threads']
        #mega_memory_per_thread = PARAMS['megahit_meta_memory']

        # Parse and compute total memory in GB
        #memory_value = int(mega_memory_per_thread.strip().upper().replace('G', ''))
        #total_memory = mega_threads * memory_value

        # Dynamically compute run options
        #run_options = f"-t {mega_threads} -m {total_memory}"
                
        # Megahit is randomly failing when large numbers of jobs submitted
        # /usr/lib64/libgomp.so.1: version `GOMP_4.0' not found
        delay = random.randint(1, 1000)
        statement1 = ("sleep {delay} &&"
                      " megahit"
                     # "  {run_options}"
                      "  {libraries}"
                      "  -o {out_sub_dir}"
                      "  2> {outfile}.log".format(**locals()))


        # Fetch contigs
        contigs = os.path.join('..', out_sub_dir, 'final.contigs.fa')
        out_prefix = os.path.basename(out_sub_dir)

        cont_out = os.path.join(out_dir,
                                out_prefix + '.megahit.contigs.fasta')

        statement2 = ("ln -s {contigs} {cont_out}".format(**locals()))

        statements = ' && '.join([statement1, statement2])
        
        return statements

    
###############################################################################
# Miscellaneous

# ASSEMBLERS = {'bayesHammer': SpadesReadCorrection(),
#               'spades': Spades(),}
# def mapper(infile, outfile, ref_db, **PARAMS):
#     '''Wrapper for using different metagenome assemblers'''
#     mapping_tool = PARAMS.get('assembler', 'spades')

#     if mapping_tool == 'rtg':
#         mapper = MapWithRTG()
#     else:
#         raise Exception( "The mapping tool %s is not configured for the"
#                          " virome pipeline" % mapping_tool)

#     statement = mapper.build(infile, outfile, ref_db, **PARAMS)

#     return statement


###############################################################################
# Annotation tools

def runAnnotationTool(infile, out_fasta, out_gff, tool, **PARAMS):
    '''Run annotation tool to predict genes in assembly.
    Current options are:
                        MetaGeneMark
    '''

    if tool == 'mgm':
        annotation_tool = AnnotateWithMGM()
    elif tool == 'prodigal':
        annotation_tool = AnnotateWithProdigal()
    else:
        raise Exception("The annotation tool %s is not configured to"
                        " work with this pipeline" % tool)

    statement = annotation_tool.build(infile, out_fasta, out_gff, **PARAMS)

    return statement
    

class AnnotateMetaAssembly(object):
    '''Base class for simple handling of different annotation tools'''

    def annotate(self, infile, out_fasta, out_gff, **PARAMS):
        '''Create annotation'''
        pass
    
    def postProcess(self, infile, out_fasta, out_gff, **PARAMS):
        '''Post-processing step to extract annotated sequence'''
        pass

    def build(self, infile, out_fasta, out_gff, **PARAMS):
        cmd_annotate, cluster_options = self.annotate(infile,
                                                      out_fasta,
                                                      out_gff,
                                                      **PARAMS)
        cmd_postprocess = self.postProcess(infile,
                                           out_fasta,
                                           out_gff,
                                           **PARAMS)

        statement = [x for x in [cmd_annotate, cmd_postprocess] if x]
        statement = " && ".join(statement)

        return statement, cluster_options 

    
class AnnotateWithMGM(AnnotateMetaAssembly):
    '''Run MetaGeneMark'''

    def annotate(self, infile, out_fasta, out_gff, **PARAMS):
        # mgm by default doesn't output compressed files
        out_fasta = out_fasta[:-3]
        out_gff = out_gff[:-3]
        
        # if the model file is missing from PARAMS, then assume it's in the
        # parent directory of the tool. Assumed to be MetaGeneMark_v1.mod
        df_modfile = os.path.dirname(spawn.find_executable('gmhmmp'))
        df_modfile = os.path.join(os.path.dirname(df_modfile),
                                  'MetaGeneMark_v1.mod')
        dfm = PARAMS.get('metagenemark_model_file', df_modfile)
        if dfm:
            df_modfile = dfm
            
        E.info(df_modfile)
        
        options = PARAMS['metagenemark_options']
        run_options = PARAMS['metagenemark_run_options']
        
        statement = ("gmhmmp "
                     " -m {df_modfile}"
                     " {options}"
                     " -f G -o {out_gff}"
                     " -d -D {out_fasta}"
                     " {infile}"
                     " &> {out_fasta}.log".format(**locals()))
        
        return statement, run_options

    def postProcess(self, infile, out_fasta, out_gff, **PARAMS):
        # mgm by default doesn't output compressed files
        out_fasta = out_fasta[:-3]
        out_gff = out_gff[:-3]

        statement = ("sed -i 's/[[:space:]]\+>/;/' {out_fasta} &&"
                     " gzip {out_fasta}; gzip {out_gff}".format(**locals()))

        return statement


class AnnotateWithProdigal(AnnotateWithMGM):
    '''Run Prodigal v2.*'''

    def annotate(self, infile, out_fasta, out_gff, **PARAMS):
        out_fasta = out_fasta[:-3]
        out_gff = out_gff[:-3]

        options = PARAMS['prodigal_options']
        run_options = PARAMS['prodigal_run_options']

        statement = ("prodigal "
                     " {options}"
                     " -i {infile}"
                     " -o {out_gff}"
                     " -f gff"
                     " -d {out_fasta}"
                     " -p meta"
                     " &> {out_fasta}.log".format(**locals()))

        return statement, run_options

    
    def postProcess(self, infile, out_fasta, out_gff, **PARAMS):
        '''An attempt to intelligently compress all the information in
        the fasta header so there is no whitespace'''
        out_fasta = out_fasta[:-3]
        out_gff = out_gff[:-3]


        out_complete = out_fasta + '_partial.txt'
        statement = (" cat {out_fasta} |"
                     " grep '>' |"
                     " awk -F \\\"partial=\\\" '{print $2}' |"
                     " awk -F \\\";\\\" '{print $1}' > {out_complete} &&"
                     " sed -i 's/[[:space:]]/;/' {out_fasta} &&"
                     " sed -i 's/#[[:space:]]/#/g' {out_fasta} &&"
                     " sed -i 's/[[:space:]]//g' {out_fasta} &&"
                     " gzip {out_fasta}; gzip {out_gff}".format(**locals()))

        return statement
    
    
###############################################################################
# Clustering Tools

def runClusteringTool(infile, out_fasta, tool, **PARAMS):
    '''Run clustering tool to collapse redundant sequences.
    Current options are:
                        CD-HIT
                        USEARCH
    '''

    if tool == 'cdhit':
        annotation_tool = ClusterWithCDHIT()
    elif tool == 'usearch':
        annotation_tool = ClusterWithUSEARCH()
    else:
        raise Exception("The clustering tool %s is not configured to"
                        " work with this pipeline" % tool)

    statement = annotation_tool.build(infile, out_fasta, **PARAMS)

    return statement


class ClusterMetaAssemblyAnnotations(object):
    '''Base class for simple handling of different annotation tools'''

    def clean(self, infile):
        '''There is whitespace in mgm output'''
        temp_file = P.getTempFilename('.')
        os.unlink(temp_file)
        temp_file = temp_file + '.fna'

        statement = "zcat {infile} | sed '/^$/d' > {temp_file}".format(**locals())

        return statement, temp_file
    
    def cluster(self, infile, out_fasta, **PARAMS):
        '''Cluster redundant sequences'''
        pass
    
    def postProcess(self, infile, temp_file, out_fasta, **PARAMS):
        '''Post-processing step to extract non-redundant sequences'''

        if out_fasta.endswith('.gz'):
            out_fasta = P.snip(out_fasta, '.gz')
        statement = ("rm -f {temp_file} &&"
                     " gzip {out_fasta}".format(**locals()))

        return statement
        
    def build(self, infile, out_fasta, **PARAMS):
        cmd_clean, temp_file = self.clean(infile)
        cmd_annotate, cluster_options = self.cluster(infile,
                                                     temp_file,
                                                     out_fasta,
                                                     **PARAMS)
        cmd_postprocess = self.postProcess(infile,
                                           temp_file,
                                           out_fasta,
                                           **PARAMS)

        statement = [x for x in [cmd_clean,
                                 cmd_annotate,
                                 cmd_postprocess] if x]
        statement = " && ".join(statement)

        return statement, cluster_options 


class ClusterWithCDHIT(ClusterMetaAssemblyAnnotations):
    '''Perform clustering with CDHIT'''

    def cluster(self, infile, temp_file, out_fasta, **PARAMS):
        # CD-HIT doesn't output compressed files
        out_fasta = out_fasta[:-3]
        
        run_options = PARAMS['cdhit_run_options']
        cdhit_options = PARAMS['cdhit_options']
        
        statement = ("cd-hit-est"
                     " -i {temp_file}"
                     " -o {out_fasta}"
                     " {cdhit_options}"
                     " -sc"
                     " -sf "
                     " &> {out_fasta}.log".format(**locals()))

        return statement, run_options

        # print '\n\n\n\n'
        # print statement
        # print run_options

class ClusterWithUSEARCH(ClusterMetaAssemblyAnnotations):
    '''Perform clustering with USEARCH
    Parameters are those suggested by Damian Rafal Plichta
    '''
    
    def cluster(self, infile, temp_file, out_fasta, **PARAMS):
        temp_file2 = P.getTempFilename('.')
        os.unlink(temp_file2)
        temp_file2 = temp_file2 + '.fna'

        centroids = P.snip(out_fasta, '.gz')
        consensus = re.sub('_geneset.', '_consensus.', centroids)
        uc_tab = P.snip(out_fasta, '.fasta.gz') + '.uc'

        run_options = PARAMS['usearch_run_options']
        
        statement = (" usearch"
                     "  -sortbylength {temp_file}"
                     "  -fastaout {temp_file2}"
                     "  -minseqlength {usearch__min_length}"
                     " &>> {outfile}.log &&"
                     # cluster annotations...
                     " usearch -cluster_fast {temp_file2}"
                     "  -id {usearch_id}"
                     "  {usearch_options}"
                     "  -centroids {centroids}"
                     "  -consout {consensus}"
                     "  -uc {uc_tab}"
                     "  &>> {out_fasta}.log &&"
                     " gzip {centroids} &&"
                     " gzip {consensus} &&"
                     " gzip {uc_tab} &&"
                     " rm -f {temp_file2}".format(**locals()))

        return statement, run_options

    
###############################################################################
# Mapping Tools

def mapReads(infile, outfile, index, tool, geneset=False, **PARAMS):
    '''Map reads in a manner that handles 1, 2, or 3 input files'''

    if tool == 'bowtie2':
        mapping_tool = Bowtie2Mapper()
    else:
        raise Exception("Currently only bowtie2 is supported for mapping")

    statement = mapping_tool.build(infile, outfile, index, geneset, **PARAMS)

    return statement


class Bowtie2Mapper(MetaAssembler):
    '''Inherits the same structure as meta-assembly class (as it handles
    the same input data).
    '''

    def preProcess(self, infiles, outfile, tmp_dir, **PARAMS):
        '''
        detect files containing reads, determine whether this
        is a hybrid assembly
        '''
        return ""
    
    
    def readMapper(self, infiles, outfile, tmp_dir, index, geneset, **PARAMS):
        '''Receive a list of up to three infiles, representing read1,
        read2 and singletons, generate a map command for bowtie2.
        '''

        # Begin building run command
        if geneset:
            # option for more lenient mapping to gene set
            run_cmd =  "bowtie2 {}".format(PARAMS["bowtie2_geneset_options"])
        else:
            run_cmd = "bowtie2 {}".format(PARAMS["bowtie2_mapping_options"])
        
        # Determine whether input files are fasta or fastq format
        f = self.getFormat(infiles[0])
        if f in ('fastq', 'fastq.gz', 'fastq.1.gz'):
            run_cmd += " -q"
        elif f in ('fasta', 'fasta.gz', 'fasta.1.gz'):
            run_cmd += " -f"
        else:
            raise ValueError("Unrecognised file format for mapper: %s" % f)

        # Location of index
        run_cmd += " -x " + index 
        
        # Fetch libraries to be mapped
        if len(infiles) == 1:
            # check whether input files are interlaced
            if self.checkPairs(infiles[0]):
                libraries =  ' --interleaved ' + infiles[0]
            else:
                libraries = ' -U ' + infiles[0]
        else:
            libraries = zip([' -1', '-2', '-U'], infiles)
            libraries = ' '.join([' '.join(x) for x in libraries])

        run_cmd += libraries

        # Log info is written to stderr
        run_cmd += ' 2> %s' % outfile + '.log'
        
        # By default outfiles are written to stdout.
        run_cmd += " | samtools view -Sbh - > %s" % \
                   os.path.join(tmp_dir, os.path.basename(outfile))

        return run_cmd

        
    def postProcess(self, infiles, outfile, tmp_dir, **PARAMS):
        '''Optionally remove non-unique reads, and strip sequence. 
        Automatically sort and index.
        '''
        tmp_f = os.path.join(tmp_dir, os.path.basename(outfile))
        
        cmd_postprocess = "cat {tmp_f} |"
        
        if PARAMS['mapping_remove_non_unique']:
            cmd_postprocess += (" cgat bam2bam"
                                "  --method=filter"
                                "  --filter-method=unique"
                                "  --log={outfile}.log |")

        if PARAMS['mapping_strip_sequence']:
            cmd_postprocess += (" cgat bam2bam"
                                "  --strip-method=all"
                                "  --method=strip-sequence"
                                "  --log={outfile}.log |")
            
        # cmd_postprocess += ("cgat bam2bam"
        #                     " --method=set-nh"
        #                     " --log={outfile}.log |"
        cmd_postprocess += (" samtools sort -o {outfile} &&"
                            " samtools index {outfile} &&"
                            " rm -rf {tmp_dir}")

        return cmd_postprocess.format(**locals())

    
    def build(self, infile, outfile, index, geneset, **PARAMS):
        '''Create run statement for assembler'''

        # Create a temporary directory for handling mapper output
        tmp_dir = os.path.join(os.path.dirname(outfile),
                              self.getTrack(infile) + '.tmp')
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        
        infiles = self.getPairedFiles(infile)
        cmds = []
        cmd_preprocess = self.preProcess(infiles, outfile, tmp_dir, **PARAMS)
        cmd_map = self.readMapper(infiles, outfile, tmp_dir, index, geneset, **PARAMS)
        cmd_postprocess = self.postProcess(infiles, outfile, tmp_dir, **PARAMS)

        cmds = [x for x in (cmd_preprocess, cmd_map, cmd_postprocess) if x]
        statement = ' && '.join(cmds)

        return statement

def mpa2taxtree(infile, outfile, map_out):
    '''Convert an MPA formatted file output by Kraken to a tax dump file
    used by RTG map and RTG species.
    '''

    # set up a taxonomy dictionary in which key=genome, value=genome_tax
    tax_dict = {}

    for line in IOTools.openFile(infile):
        genome, taxonomy = line.strip().split('\t')
        
        # Set up a taxonomy for this genome
        genome_tax = collections.OrderedDict()
        genome_tax['d'] = 'unclassified'
        genome_tax['p'] = 'unclassified'
        genome_tax['c'] = 'unclassified'
        genome_tax['o'] = 'unclassified'
        genome_tax['f'] = 'unclassified'
        genome_tax['g'] = 'unclassified'
        genome_tax['s'] = 'unclassified'
        
        # Taxonomy is listed domain -> species, unclassified levels are missing
        for tax_level in taxonomy.split('|'):
            level, tax = tax_level.split('__')
            genome_tax[level] = tax
            
        last = genome_tax['d']
        for k, v in genome_tax.items():
            if k == 'd' and v == 'unclassified':
                break
            elif v == 'unclassified':
                genome_tax[k] = 'unclassified' + '_' + last
            else:
                last = v
                    
        tax_dict[genome] = genome_tax

    # Create a dict with unique values to a set for each taxonomic level
    taxUniq = {'domain': set(), 'phylum': set(), 'class': set(),
               'order': set(), 'family': set(), 'genus': set(),
               'species': set()}

    for genome in tax_dict.keys():
        g = tax_dict[genome]
        taxUniq['domain'].add(g['d'])
        taxUniq['phylum'].add(g['p'])
        taxUniq['class'].add(g['c'])
        taxUniq['order'].add(g['o'])
        taxUniq['family'].add(g['f'])
        taxUniq['genus'].add(g['g'])
        taxUniq['species'].add(g['s'])   

    # Iterate over each unique level, create a unique taxID number
    taxIDs = {'domain': {}, 'phylum': {}, 'class': {},
              'order': {}, 'family': {}, 'genus': {}, 'species': {}}

    n = 1
    for level in ('domain', 'phylum', 'class',
                  'order', 'family', 'genus', 'species'):
        for tax in taxUniq[level]:
            n += 1
            taxIDs[level][tax] = str(n)

    # Create a tax tree
    out_tree = set()

    map_dict = {}
    
    for genome, t in tax_dict.items():
        n += 1
        map_dict[genome] = str(n)
        
        # taxID parentID rank name
        gn = (str(n), taxIDs['species'][t['s']], 'genome', genome)
        out_tree.add(gn)
        s = (taxIDs['species'][t['s']], taxIDs['genus'][t['g']], 'species', t['s'])
        out_tree.add(s)
        g = (taxIDs['genus'][t['g']], taxIDs['family'][t['f']], 'genus', t['g'])
        out_tree.add(g)
        f = (taxIDs['family'][t['f']], taxIDs['order'][t['o']], 'family', t['f'])
        out_tree.add(f)
        o = (taxIDs['order'][t['o']], taxIDs['class'][t['c']], 'order', t['o'])
        out_tree.add(o)
        c = (taxIDs['class'][t['c']], taxIDs['phylum'][t['p']], 'class', t['c'])
        out_tree.add(c)
        p = (taxIDs['phylum'][t['p']], taxIDs['domain'][t['d']], 'phylum', t['p'])
        out_tree.add(p)
        d = (taxIDs['domain'][t['d']], '1', 'domain', t['d'])
        out_tree.add(d)

    outf =  IOTools.openFile(outfile, 'w')
    outf.write('#RTG taxonomy version 1.0\n#taxID\tparentID\trank\tname\n')
    outf.write('1\t-1\tno rank\troot\n')
    for i in out_tree:
        outf.write('\t'.join(i) + '\n')
    outf.close()

    pickle.dump(map_dict, open(map_out,'wb'))


###############################################################################
# Miscellaneous Tools

def rarefyCountTable(df, minimum, drop=True):
    """Rarefy dataframe where rows are OTUs and columns are samples.
    If not minimum='min', then it must be an integer value
    WARNING: No random seed is applied"""

    if minimum != 'min':
        try:
            minimum = int(minimum)
        except ValueError:
            raise ValueError("Please ensure rarefying minimum is set to 'min',"
                             " or an integer value")
        df = df.iloc[:, list(df.sum() > minimum)]
    else:
        minimum = min(df.sum())

    E.info("Rarefying count data to minimum depth %i" % minimum)

    df = df.apply(lambda x: skbio.stats.subsample_counts(x, minimum))

    if drop:
        # Drop rows that are now all zeros
        df =  df[(df.T !=0).any()]
    
    return df


