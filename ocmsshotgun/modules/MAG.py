import re
import os
import ocmsshotgun.modules.Utility as utility
import cgatcore.pipeline as P

class samtools_align(utility.matchReference):
    def __init__(self, fastq1, outfile, reference_fasta, **params):
        """
        Initialize samtools_align with required arguments.
        """
        super().__init__(fastq1, outfile, **params)  # Inherit functionality from matchReference
        self.reference_fasta = reference_fasta  # Assign the reference FASTA

    def buildStatement(self):
        '''Align paired-end reads from FASTQ files to input FASTA reference and generate BAM files'''
        fastq1 = self.fastq1
        fastq2 = self.fastq2
        outfile1 = self.outfile
        reference_fasta = self.reference_fasta


        # Derived paths
        sample_out = P.snip(self.fastq1, self.fq1_suffix)
        logfile = sample_out + '.bam.log'
        tmp_fastq1 = P.get_temp_filename('.')
        tmp_fastq2 = P.get_temp_filename('.')

        # Parameters
        samtools_n_threads = self.PARAMS["samtools_job_threads"]
        samtools_options = self.PARAMS.get("samtools_options", "")

        # Build the alignment statement
        statement = (
            "zcat %(fastq1)s > %(tmp_fastq1)s && "
            "zcat %(fastq2)s > %(tmp_fastq2)s && "
            "bwa index  %(reference_fasta)s && "
            "bwa mem -t %(samtools_n_threads)s %(samtools_options)s %(reference_fasta)s "
            "%(tmp_fastq1)s %(tmp_fastq2)s | "
            "samtools view -Sb - | "
            "samtools sort -@ %(samtools_n_threads)s -o %(outfile1)s && "
            "rm -f %(tmp_fastq1)s %(tmp_fastq2)s && "
            "gzip %(logfile)s"
        ) % locals()

        return statement

