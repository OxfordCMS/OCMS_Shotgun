'''
fastq2filteredfastq.py
=======================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python

Purpose
-------

filter a fastq file based on reads in a bam file.

Usage
-----

Example::

   zcat in.fastq.gz | python fastq2filteredfastq.py --bamfile=x.bam --outfile=filtered.fastq

Type::

   python fastq2filteredfastq.py --help

for command line help.

Command line options
--------------------

'''

import os
import sys
import pysam
import warnings
from cgatcore import pipeline as P
from cgatcore import experiment as E
from cgatcore import iotools as IOTools
from cgat import Fastq as Fastq
import glob

class bamFilter():
    def __init__(self, bamfile, paired, filter_mapping):
        self.bamfile = bamfile
        self.paired = paired
        self.filter_mapping = filter_mapping
        self.prefix = os.path.basename(bamfile.rstrip(".bam"))

    def filter_mapping_se(self, alignment):
        """Helper function pysam read object that are 
        mapped or unmapped for single end reads"""
        if self.filter_mapping == 'unmapped':
            condition = alignment.is_unmapped
        elif self.filter_mapping == 'mapped':
            condition = alignment.is_mapped
        
        # only returns something when condition is met
        if condition:
            return {alignment.qname: "read1"} 
        else:
            return

    def filter_mapping_pe(self, alignment, mate):
        """Helper function to return pysam alignment query name and a label
        of whether is forward, reverse, or singleton"""
        if self.filter_mapping == 'unmapped':
            condition1 = alignment.is_unmapped
            condition2 = mate.is_unmapped
        elif self.filter_mapping == 'mapped':
            condition1 = alignment.is_mapped
            condition2 = mate.is_mapped
        
        # both reads mapped (or unmapped)
        if condition1 and condition2:
            out = {alignment.qname: "read1", mate.qname: "read2"}
        # one read mapped (or unmapped)
        elif condition1 or condition2:
            out = {alignment.qname: "singleton"}
        # neither reads mapped (or unmapped)
        else:
            return

        return out

    def bam_filter_mapping(self):
        """
        Iterates through bam file and stores all mapped alignments in a set
        """
        # initialize dictionary of mapped alignments
        mapped = {}

        # temporary dictionary of unmated alignments
        unmated = {}

        # initialize dictionary to count reads
        self.counter = E.Counter()
        self.counter.input_alignments = 0
        self.counter.input_reads = 0
        self.counter.output_reads = 0
        self.counter.singleton_reads = 0

        with pysam.AlignmentFile(self.bamfile, "rb") as bamfile:
            for alignment in bamfile.fetch(until_eof=True):
                self.counter.input_alignments += 1
                qname = alignment.qname
                
                # paired end reads
                if self.paired:
                    # if both mates of pe reads found, delete from unmated and 
                    # process the pair to filter for mapped
                    if qname in unmated:
                        mate = unmated.pop(qname)
                        filtered_qname = self.filter_mapping_pe(
                            alignment, mate)
                    # if haven't found mate yet, keep in dictionary temporarily
                    else:
                        unmated[qname] = alignment
                        filtered_qname = {}
                # single end reads
                else:
                    filtered_qname = self.filter_mapping_se(alignment)

                # add mapped qname to dictionary if read is mapped
                if filtered_qname:
                    mapped.update(filtered_qname)
                    self.counter.output_reads += len(filtered_qname.keys())

            # Process remaining singletons at the end
            for remaining_read in unmated.values():
                filtered_singleton = self.filter_mapping_se(remaining_read)
                if filtered_singleton:
                    self.counter.singleton_reads += 1
                    mapped.update(filtered_singleton)
        
        return mapped
                
    def bamfiltered2fastq(self, in_fastq, outdir):
        """Iterate through bam file, use bam_filter_mapped to get a diciontary
        of alignments that have mapped. Iterates through fastq file and checks
        if fastq record is in the mapped dicitonary. Keeps or discards record
        depending on filter_mapping argument. If want unmapped reads, will keep
        fastq record if fastq identifier is not in mapped dictionary. If want 
        mapped reads, will keep fastq record if fastq identifier is in mapped 
        dictionary.

        Writes filtered compressed fastqs in specified directory. If working 
        with paired end reads, will write each read in respective fastq file, 
        with forward read in fastq.1.gz, reverse read in fastq.2.gz, singletons
        in fastq.3.gz"""
        summary = IOTools.open_file(f"{outdir}/{self.prefix}_filter_summary.tsv", 'w')

        # Open output fastq files for Read 1, Read 2, and 
        # Singletons (for paired-end reads)
        if self.paired:
            fastq1 = f"{outdir}/{self.prefix}_{self.filter_mapping}.fastq.1.gz"
            fastq2 = f"{outdir}/{self.prefix}_{self.filter_mapping}.fastq.2.gz"
            fastq3 = f"{outdir}/{self.prefix}_{self.filter_mapping}.fastq.3.gz"
            fastq1_out = IOTools.open_file(fastq1, 'w')
            fastq2_out = IOTools.open_file(fastq2, 'w')
            fastq3_out = IOTools.open_file(fastq3, 'w')
        else:
            fastq1 = f"{outdir}/{self.prefix}_{self.filter_mapping}.fastq.gz"
            fastq_out = IOTools.open_file(fastq1, 'w')
    
        # get dictionary of mapped reads
        filtered_qnames = self.bam_filter_mapping()
        
        # if no mapped reads found, symlink in fastq file to outdir
        if not filtered_qnames:
            indir = os.path.dirname(in_fastq)
            fastq_ins = os.glob.glob(f"{indir}/{self.prefix}*.fastq*.gz")
            fastq_outs = [fq.replace(indir, outdir) for fq in fastq_ins]
            for fq_in, fq_out in zip(fastq_ins, fastq_outs):
                os.symlink(fq_in, fq_out)
            E.info(f"No mapped reads for {in_fastq} in {self.bamfile}")

        # stream in fastq file
        for record in Fastq.iterate(IOTools.open_file(in_fastq)):
            self.counter.input_reads += 1
            # strip record header of spaces
            id = record.identifier.split(" ")[0]

            if id in filtered_qnames.keys():        
                # Handle single-end reads
                if not self.paired:
                    fastq_out.write(record.__str__()+ "\n")  
                # Handle paired-end reads
                else:
                    if filtered_qnames[id] == 'read1':
                        fastq1_out.write(record.__str__() + "\n")  
                    elif filtered_qnames[id] == 'read2':
                        fastq2_out.write(record.__str__() + "\n") 
                    elif filtered_qnames[id] == 'singleton':
                        fastq3_out.write(record.__str__() + "\n")
            else:
                continue

        # write summary file
        header, entry = zip(*self.counter.iteritems())
        entry = [str(x) for x in entry]
        header = '\t'.join(header)
        entry = '\t'.join(entry)
        summary.write(header + "\n")
        summary.write(entry + "\n")
        summary.close()
        # Close the fastq files
        if self.paired:
            fastq1_out.close()
            fastq2_out.close()
            fastq3_out.close()
        else:
            fastq_out.close()

def main(argv=None):
    """
    Writes filtered reads (mapped/unmapped) to separate FastQ files:
    - For paired-end reads: fastq.1.gz, fastq.2.gz, and fastq.3.gz for singletons.
    - For single-end reads: a single fastq.gz file.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.ArgumentParser(description = __doc__)

    parser.add_argument("-i", "--in_fastq", dest="in_fastq", type=str,
                        help="Fastq file to filter. Can be file or stdin.\
                            Should be the fastq used generate bam alignments. \
                            For paired end reads only provide fastq.1.gz")
    parser.add_argument("-o", "--outdir", dest="outdir", type=str,
                        help="Output directory for bam filtered fastqs"),
    parser.add_argument("-b", "--bamfile", dest="bamfile", type=str,
                        help="input bamfile to filter reads from")
    parser.add_argument("-f", "--filter_mapping", dest="filter_mapping", 
                        default='unmapped', const='unmapped', nargs='?', 
                        type=str, choices=['unmapped', 'mapped'],
                        help="type of read to keep")
    parser.add_argument("--paired", dest="paired", action="store_true", 
                        help="paired end data. if set to will give mapped \
                            reads as fastq.1 and fastq.2 of forward and reverse \
                            reads respectively and fastq.3 of singletons")

    # add common options (-h/--help, ...) and parse command line
    (args) = E.start(parser, argv=argv)

    E.info("iterating over bam file")
    E.info("filtering fastq file")

    # set up generator for bam file and alignments
    filter_reads = bamFilter(args.bamfile, args.paired, args.filter_mapping)
    
    # iterate through bam, filter, and, write filtered reads to fastq
    filter_reads.bamfiltered2fastq(args.in_fastq, args.outdir)
    
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))