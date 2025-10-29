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
    """
    ExtractRefinedBinReads tool class.

    Extracts read IDs for contigs belonging to refined bins
    using alignments in a BAM file. Creates one text file per
    refined bin listing all read IDs.
    """

    def __init__(self, inputs, **PARAMS):
        """
        Parameters
        ----------
        inputs : tuple
            (Sorted BAM file, contig→bin JSON file)
        PARAMS : dict
            Pipeline parameters (threads, memory, etc.)
        """
        # Pass a dummy outfile so Utility.BaseTool initializes cleanly, as this requires outfile to be valid filesystem paths (strings or Path objects).
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
            raise FileNotFoundError(
                f"Mapping file not found for {self.sample_id}: {self.map_file}"
            )

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
                if read.is_unmapped:
                    unmapped_handle.write(read.query_name + "\n")
                    read_counts["unmapped"] += 1
                    continue

                ref_name = bam_in.get_reference_name(read.reference_id) #retrive contig names this reads aligned to
                if ref_name not in valid_contigs: #skip if the reads not aligned to contigs resulted in refined bin
                # mapped but not to any refined-bin contig
                    unassigned_handle.write(read.query_name + "\n")
                    read_counts["unassigned"] += 1
                    continue
                    
                # mapped to a refined-bin contig
                bin_name = contig_to_bin[ref_name] #look up which refined bin that contig belong to 
                handle = bin_to_handles.get(bin_name) #retrieves the open file handle (e.g., the output file for that bin) from the dictionary bin_to_handles using the key bin_name
                if handle: #if a valid file handle exists for this bin, write the current read’s ID (name) to that bin’s output file
                    handle.write(read.query_name + "\n")
                    read_counts[bin_name] += 1

        # --- Close output handles ---
        for handle in bin_to_handles.values():
            handle.close()
        unmapped_handle.close()
        unassigned_handle.close()

        # --- Sanity check for missing contigs ---
        with pysam.AlignmentFile(self.bam_file, "rb") as bam_in:
            bam_contigs = set(bam_in.references) #gives you a set of unique contig names for quick, reference stores list of contigs/chromosome stored in a BAM file
        missing_contigs = [c for c in contig_to_bin if c not in bam_contigs] #identifies contigs that are present in your refined-bin mapping but missing from the sample’s actual BAM file

        if missing_contigs:
            print(f"[WARN] {self.sample_id}: {len(missing_contigs)} contigs not in BAM header")
        else:
            print(f"[OK] {self.sample_id}: all refined-bin contigs found in BAM")

        print(f"[{self.sample_id}] Extracted read IDs for {len(bin_to_contigs)} refined bins.")
        
        # --- Write summary TSV report ---
        summary_path = self.outdir / f"{self.sample_id}_read_summary.tsv"
        with open(summary_path, "w") as s:
            s.write("bin_name\tread_count\n")
            for bin_name, count in read_counts.items():
                s.write(f"{bin_name}\t{count}\n")
            s.write(f"total\t{total_reads}\n")

        print(f"[{self.sample_id}] Extracted read IDs for {len(bin_to_contigs)} bins.")
        print(f"[{self.sample_id}] Summary written to: {summary_path}")
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

