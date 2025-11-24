# Module: MetaRefinement
# Description:
#     Tools for metagenomic bin refinement steps:
#       1. ExtractRefinedBinReads — extract read IDs per refined bin.
#       2. FilterFastqByIds — filter original FASTQs using those read IDs.

import ocmstoolkit.modules.Utility as Utility
from cgatcore import pipeline as P
from pathlib import Path
from Bio import SeqIO
import gzip
import json
import pysam

# ---------------------------------------------------------------------------
# Class 1: ExtractRefinedBinReads
# ---------------------------------------------------------------------------

class ExtractRefinedBinReads(Utility.BaseTool):
    def __init__(self, inputs, **PARAMS):
        """
        Parameters
        ----------
        inputs : tuple
            (Sorted BAM file, contig→bin JSON file)
        PARAMS : dict
            Pipeline parameters (threads, memory, etc.)
        """
        dummy_outfile = f"03_refined_bin_reads.dir/{Path(inputs[0]).stem}_placeholder.txt"
        super().__init__(inputs[0], dummy_outfile, **PARAMS)

        self.bam_file, self.map_file = map(Path, inputs)
        self.sample_id = self.bam_file.stem.replace("_sorted", "")
        self.outdir = Path("03_refined_bin_reads.dir") / self.sample_id
        self.outdir.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------------
    # Core logic: extract read IDs from BAM
    # ---------------------------------------------------------------------
    def _extract_reads(self):
        """Perform read ID extraction from BAM."""
        print(f"[{self.sample_id}] Extracting read IDs from {self.bam_file.name}")

        # --- Load contig → bin mapping ---
        if not self.map_file.exists():
            raise FileNotFoundError(f"Mapping file not found for {self.sample_id}: {self.map_file}")

        with open(self.map_file) as f:
            contig_to_bin = json.load(f)

        # --- Group contigs by bin ---
        bin_to_contigs = {}
        for contig, bin_name in contig_to_bin.items():
            bin_to_contigs.setdefault(bin_name, []).append(contig)
        valid_contigs = set(contig_to_bin.keys())

        # --- Prepare output files ---
        bin_to_handles = {}
        for bin_name in bin_to_contigs:
            out_path = self.outdir / f"{bin_name}_read_ids.txt"
            bin_to_handles[bin_name] = open(out_path, "w")

        # Additional output files for unmapped/unassigned reads
        unmapped_path = self.outdir / "unmapped_read_ids.txt"
        unassigned_path = self.outdir / "unassigned_read_ids.txt"
        unmapped_handle = open(unmapped_path, "w")
        unassigned_handle = open(unassigned_path, "w")

        # --- Initialize counters ---
        read_counts = {bin_name: 0 for bin_name in bin_to_contigs}
        read_counts["unmapped"] = 0
        read_counts["unassigned"] = 0
        total_reads = 0

        # --- Iterate over BAM and write read IDs ---
        with pysam.AlignmentFile(self.bam_file, "rb") as bam_in:
            for read in bam_in.fetch(until_eof=True):
                total_reads += 1

                # Define suffix for directionality
                if read.is_read1:
                    suffix = "/1"
                elif read.is_read2:
                    suffix = "/2"
                else:
                    suffix = ""

                # Unmapped reads
                if read.is_unmapped:
                    unmapped_handle.write(read.query_name + suffix + "\n")
                    read_counts["unmapped"] += 1
                    continue

                # Reads mapped but not to valid contigs
                ref_name = bam_in.get_reference_name(read.reference_id)
                if ref_name not in valid_contigs:
                    unassigned_handle.write(read.query_name + suffix + "\n")
                    read_counts["unassigned"] += 1
                    continue

                # Mapped to a refined-bin contig
                bin_name = contig_to_bin[ref_name]
                handle = bin_to_handles.get(bin_name)
                if handle:
                    handle.write(read.query_name + suffix + "\n")
                    read_counts[bin_name] += 1

        # --- Close output handles ---
        for handle in bin_to_handles.values():
            handle.close()
        unmapped_handle.close()
        unassigned_handle.close()

        # --- Sanity check for missing contigs ---
        with pysam.AlignmentFile(self.bam_file, "rb") as bam_in:
            bam_contigs = set(bam_in.references)
        missing_contigs = [c for c in contig_to_bin if c not in bam_contigs]

        if missing_contigs:
            print(f"[WARN] {self.sample_id}: {len(missing_contigs)} contigs not in BAM header")
        else:
            print(f"[OK] {self.sample_id}: all refined-bin contigs found in BAM")

        print(f"[{self.sample_id}] Extracted read IDs for {len(bin_to_contigs)} refined bins.")

        # --- Write summary TSV report ---
        summary_path = self.outdir / f"{self.sample_id}_read_summary.tsv"

        # Calculate totals
        total_reads = read_counts.get("unmapped", 0) + read_counts.get("unassigned", 0)
        mapped_reads = 0
        for key, val in read_counts.items():
            if key not in ["unmapped", "unassigned"]:
                mapped_reads += val
        total_reads += mapped_reads

        # Compute summary lines
        with open(summary_path, "w") as s:
            s.write("bin_name\tread_count\t%_of_mapped_reads\t%_of_total_reads\n")

            # Each refined bin
            for bin_name, count in sorted(
                ((b, c) for b, c in read_counts.items() if b not in ["unmapped", "unassigned"]),
                key=lambda x: x[1],
                reverse=True,
            ):
                pct_mapped = (100 * count / mapped_reads) if mapped_reads > 0 else 0
                pct_total = (100 * count / total_reads) if total_reads > 0 else 0
                s.write(f"{bin_name}\t{count}\t{pct_mapped:.2f}\t{pct_total:.2f}\n")

            # Corrected combined totals
            pct_bins_mapped = 100 * mapped_reads / (mapped_reads + read_counts.get("unassigned", 0))
            pct_bins_total = 100 * mapped_reads / total_reads
            s.write(f"All_refined_bins_combined\t{mapped_reads}\t{pct_bins_mapped:.2f}\t{pct_bins_total:.2f}\n")

            pct_unassigned_mapped = 100 * read_counts.get("unassigned", 0) / (mapped_reads + read_counts.get("unassigned", 0))
            pct_unassigned_total = 100 * read_counts.get("unassigned", 0) / total_reads
            s.write(f"unassigned\t{read_counts.get('unassigned', 0)}\t{pct_unassigned_mapped:.2f}\t{pct_unassigned_total:.2f}\n")

            pct_unmapped_total = 100 * read_counts.get("unmapped", 0) / total_reads
            s.write(f"unmapped\t{read_counts.get('unmapped', 0)}\t\t{pct_unmapped_total:.2f}\n")

            s.write(f"total\t{total_reads}\t\t100.00\n")

        print(f"[{self.sample_id}] Read summary with percentages written to: {summary_path}")


    # ---------------------------------------------------------------------
    # Statement builder for CGAT-core P.run()
    # ---------------------------------------------------------------------
    def build_statement(self):
        """Return command string for pipeline execution."""
        statement = f"""
        python -c "from ocmsshotgun.modules.MetaRefinement import ExtractRefinedBinReads;
tool = ExtractRefinedBinReads(('{self.bam_file}', '{self.map_file}'));
tool._extract_reads()"
        """
        return statement.strip()


# ---------------------------------------------------------------------------
# Class 2: FilterFastqByIds
# ---------------------------------------------------------------------------

class FilterFastqByIds(Utility.BaseTool):
    """
    FilterFastqByIds tool class.

    Filters paired-end FASTQ files using read ID lists produced in previous
    steps of the pipeline. Produces filtered FASTQs for each refined bin.
    """

    def __init__(self, infile, outfile, **PARAMS):
        """
        Parameters
        ----------
        infile : str
            Path to read ID list (.txt)
        outfile : list[str]
            List of two output FASTQs [R1, R2]
        PARAMS : dict
            Pipeline parameters (threads, memory, etc.)
        """
        # BaseTool expects one output path → use R1
        super().__init__(infile, outfile[0], **PARAMS)

        self.read_id_file = Path(infile)
        self.out_r1, self.out_r2 = map(Path, outfile)

        # Derive sample info
        self.sample_id = self.out_r1.parent.name
        self.bin_name = self.out_r1.stem.replace("_R1", "")
        self.fq_dir = Path("input_fastqs.dir")

        # Build full paths to input FASTQs
        self.r1_path = self.fq_dir / f"{self.sample_id}.fastq.1.gz"
        self.r2_path = self.fq_dir / f"{self.sample_id}.fastq.2.gz"

        if not self.r1_path.exists() or not self.r2_path.exists():
            raise FileNotFoundError(
                f"FASTQ files for {self.sample_id} not found in {self.fq_dir}"
            )

    # ---------------------------------------------------------------------
    # Core logic: filter FASTQs
    # ---------------------------------------------------------------------
    def _filter_fastqs(self):
        """Perform FASTQ filtering."""
        self.out_r1.parent.mkdir(parents=True, exist_ok=True)
        print(f"[{self.sample_id}] Processing refined bin: {self.bin_name}")

        # Load read IDs
        with open(self.read_id_file) as f:
            read_ids = set(line.strip() for line in f if line.strip())
        print(f"  Loaded {len(read_ids)} read IDs")

        if not read_ids:
            print(f"  [WARN] No read IDs found, skipping bin")
            return

        # Helper for filtering one FASTQ file
        def filter_fastq(in_path, out_path):
            count_in, count_out = 0, 0
            with gzip.open(in_path, "rt") as in_fq, gzip.open(out_path, "wt") as out_fq:
                for rec in SeqIO.parse(in_fq, "fastq"):
                    count_in += 1
                    if rec.id in read_ids:
                        SeqIO.write(rec, out_fq, "fastq")
                        count_out += 1
            return count_in, count_out

        # Run filtering
        c1_in, c1_out = filter_fastq(self.r1_path, self.out_r1)
        c2_in, c2_out = filter_fastq(self.r2_path, self.out_r2)

        print(f"  R1: {c1_out}/{c1_in} reads written")
        print(f"  R2: {c2_out}/{c2_in} reads written")
        print(f"  [OK] Wrote FASTQs: {self.out_r1.name}, {self.out_r2.name}")

    # ---------------------------------------------------------------------
    # Statement builder for CGAT-core P.run()
    # ---------------------------------------------------------------------
    def build_statement(self):
        """Return shell command to run this class via Python inline execution."""
        statement = f"""
        python -c "from ocmsshotgun.modules.MetaRefinement import FilterFastqByIds;
tool = FilterFastqByIds('{self.read_id_file}', ['{self.out_r1}', '{self.out_r2}']);
tool._filter_fastqs()"
        """
        return statement.strip()

