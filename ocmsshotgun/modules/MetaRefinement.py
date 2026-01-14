"""
MetaRefinement.py

Single-file module for refined-bin read extraction and paired FASTQ generation.

Key behaviours:
 - Option A: if either mate maps to a refined contig, recruit BOTH mates.
 - FASTQ extraction uses simultaneous paired-iteration (4-line records), enforcing strict 1:1 pairing.
 - Normalization is minimal and tuned to your data: only strip leading '@' from FASTQ headers.
"""

from pathlib import Path
from typing import Tuple, Set, Dict, Iterable
import gzip
import shlex
import json
import pysam
import subprocess
import os

# ----------------------------
# Utilities
# ----------------------------
def normalize_read_name(raw: str) -> str:
    """
    Minimal normalization tailored to your dataset:
      - Remove leading '@' from FASTQ header lines (BAM names do not have '@').
    """
    if not raw:
        return raw
    if raw.startswith("@"):
        raw = raw[1:]
    return raw

def ensure_gz_path(p: Path) -> Path:
    """Ensure an output path ends with .gz."""
    return p if str(p).endswith(".gz") else Path(str(p) + ".gz")

# ----------------------------
# Class: ExtractRefinedBinReads
# ----------------------------
class ExtractRefinedBinReads:
    """
    Extract base read names from a BAM file grouped by refined bins.

    Produces:
      03_refined_bin_reads.dir/<sample_id>/
        <bin_name>_read_ids.txt
        unmapped_read_ids.txt
        unassigned_read_ids.txt
        <sample_id>_read_summary.tsv
    """

    def __init__(self, inputs: Tuple[str, str], **PARAMS):
        """
        inputs: (bam_path, contig_to_bin_json)
        PARAMS: optional pipeline parameters (kept for API parity)
        """
        self.bam_path = Path(inputs[0])
        self.map_path = Path(inputs[1])
        self.PARAMS = PARAMS

        # Derive sample id from BAM filename; remove _sorted suffix if present
        self.sample_id = self.bam_path.stem.replace("_sorted", "")
        self.outdir = Path("03_refined_bin_reads.dir") / self.sample_id
        self.outdir.mkdir(parents=True, exist_ok=True)

    def _load_map(self) -> Dict[str, str]:
        if not self.map_path.exists():
            raise FileNotFoundError(f"Contig->bin JSON not found: {self.map_path}")
        with open(self.map_path) as fh:
            contig_to_bin = json.load(fh)
        return contig_to_bin

    def _extract_reads(self):
        contig_to_bin = self._load_map()

        # Precompute bin->contig mapping
        bin_to_contigs: Dict[str, Set[str]] = {}
        for contig, binname in contig_to_bin.items():
            bin_to_contigs.setdefault(binname, set()).add(contig)
        valid_contigs = set(contig_to_bin.keys())

        # Accumulate base read names per bin (Option A: base names only)
        bin_to_ids: Dict[str, Set[str]] = {b: set() for b in bin_to_contigs}
        unmapped_ids: Set[str] = set()
        unassigned_ids: Set[str] = set()
        total_reads = 0

        with pysam.AlignmentFile(str(self.bam_path), "rb") as bam_in:
            for read in bam_in.fetch(until_eof=True):
                total_reads += 1
                base = normalize_read_name(read.query_name)

                if read.is_unmapped:
                    unmapped_ids.add(base)
                    continue

                try:
                    ref_name = bam_in.get_reference_name(read.reference_id)
                except Exception:
                    ref_name = None

                if not ref_name or ref_name not in valid_contigs:
                    unassigned_ids.add(base)
                    continue

                binname = contig_to_bin[ref_name]
                bin_to_ids.setdefault(binname, set()).add(base)

        # Write per-bin id files
        for binname, ids in bin_to_ids.items():
            out_path = self.outdir / f"{binname}_read_ids.txt"
            with open(out_path, "w") as out_f:
                for rid in sorted(ids):
                    out_f.write(rid + "\n")

        # Write unmapped/unassigned lists
        with open(self.outdir / "unmapped_read_ids.txt", "w") as f:
            for rid in sorted(unmapped_ids):
                f.write(rid + "\n")

        with open(self.outdir / "unassigned_read_ids.txt", "w") as f:
            for rid in sorted(unassigned_ids):
                f.write(rid + "\n")

        # Write summary TSV
        mapped_reads = sum(len(s) for s in bin_to_ids.values())
        total_collected = mapped_reads + len(unassigned_ids) + len(unmapped_ids) or 1

        summary_path = self.outdir / f"{self.sample_id}_read_summary.tsv"
        with open(summary_path, "w") as s:
            s.write("bin_name\tread_count\t%_of_mapped_reads\t%_of_total_reads\n")
            for binname, ids in sorted(bin_to_ids.items(), key=lambda x: len(x[1]), reverse=True):
                count = len(ids)
                pct_mapped = (100.0 * count / mapped_reads) if mapped_reads else 0.0
                pct_total = 100.0 * count / total_collected
                s.write(f"{binname}\t{count}\t{pct_mapped:.2f}\t{pct_total:.2f}\n")

            pct_bins_mapped = 100.0 * mapped_reads / (mapped_reads + len(unassigned_ids)) if (mapped_reads + len(unassigned_ids)) else 0.0
            pct_bins_total = 100.0 * mapped_reads / total_collected
            s.write(f"All_refined_bins_combined\t{mapped_reads}\t{pct_bins_mapped:.2f}\t{pct_bins_total:.2f}\n")

            pct_unassigned_mapped = 100.0 * len(unassigned_ids) / (mapped_reads + len(unassigned_ids)) if (mapped_reads + len(unassigned_ids)) else 0.0
            pct_unassigned_total = 100.0 * len(unassigned_ids) / total_collected
            s.write(f"unassigned\t{len(unassigned_ids)}\t{pct_unassigned_mapped:.2f}\t{pct_unassigned_total:.2f}\n")

            pct_unmapped_total = 100.0 * len(unmapped_ids) / total_collected
            s.write(f"unmapped\t{len(unmapped_ids)}\t\t{pct_unmapped_total:.2f}\n")

            s.write(f"total\t{total_collected}\t\t100.00\n")

        print(f"[OK] {self.sample_id}: per-bin read-id lists -> {self.outdir}")
        print(f"[OK] summary -> {summary_path}")

    # ------------------------------
    # FIXED: build_statement at class level
    # ------------------------------
    def build_statement(self):
        cmd = (
            f'python3 -c "'
            f'from ocmsshotgun.modules.MetaRefinement import ExtractRefinedBinReads; '
            f"tool = ExtractRefinedBinReads(('{self.bam_path}', '{self.map_path}')); "
            f'tool._extract_reads()"'
        )
        return cmd

# ----------------------------
# Class: FilterFastqByIds (paired lockstep)
# ----------------------------
class FilterFastqByIds:
    """
    Given a read-id file (base names), extract paired FASTQs by iterating R1 and R2
    simultaneously (4-line FASTQ records). Only writes pairs where the base id is in the set.
    """

    def __init__(self, read_id_file: str, outfiles: Tuple[str, str], **PARAMS):
        self.read_id_file = Path(read_id_file)
        self.out_r1 = Path(outfiles[0])
        self.out_r2 = Path(outfiles[1])
        self.PARAMS = PARAMS

        if "general_input_fastqs_dir" not in PARAMS:
            raise ValueError("PARAMS must include 'general_input_fastqs_dir'")
        self.fastq_dir = Path(PARAMS["general_input_fastqs_dir"])

        parent = self.out_r1.parent
        self.sample_id = parent.name if parent.name else self.read_id_file.stem

        # expected input naming convention: <sample>.fastq.1.gz and .fastq.2.gz
        self.r1_in = self.fastq_dir / f"{self.sample_id}.fastq.1.gz"
        self.r2_in = self.fastq_dir / f"{self.sample_id}.fastq.2.gz"

        if not self.r1_in.exists() or not self.r2_in.exists():
            raise FileNotFoundError(f"Input FASTQs not found: {self.r1_in} {self.r2_in}")

        self.out_r1.parent.mkdir(parents=True, exist_ok=True)

    def _load_base_ids(self) -> Set[str]:
        ids: Set[str] = set()
        with open(self.read_id_file) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                base = normalize_read_name(line)
                if base:
                    ids.add(base)
        return ids

    def _iterate_paired_fastq(self, r1_path: Path, r2_path: Path) -> Iterable[Tuple[str, str]]:
        """Yield pairs of FASTQ record strings (rec1, rec2)."""
        with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2:
            while True:
                r1_lines = [f1.readline() for _ in range(4)]
                r2_lines = [f2.readline() for _ in range(4)]
                if not r1_lines[0] or not r2_lines[0]:
                    break
                yield ("".join(r1_lines), "".join(r2_lines))

    def run(self):
        print(f"[INFO] FilterFastqByIds (paired lockstep) sample={self.sample_id}")
        allowed_ids = self._load_base_ids()
        print(f"[INFO] Loaded {len(allowed_ids)} allowed IDs from {self.read_id_file}")

        out_r1_path = ensure_gz_path(self.out_r1)
        out_r2_path = ensure_gz_path(self.out_r2)

        written = 0
        processed = 0

        with gzip.open(out_r1_path, "wt") as out_r1_f, gzip.open(out_r2_path, "wt") as out_r2_f:
            for rec1, rec2 in self._iterate_paired_fastq(self.r1_in, self.r2_in):
                processed += 1
                header1 = rec1.split("\n", 1)[0]
                header2 = rec2.split("\n", 1)[0]
                id1 = normalize_read_name(header1.lstrip("@"))
                id2 = normalize_read_name(header2.lstrip("@"))

                # Strict validation: IDs must match
                if id1 != id2:
                    raise ValueError(f"FASTQ files out of sync at record {processed}: R1={id1} != R2={id2}")

                if id1 in allowed_ids:
                    out_r1_f.write(rec1)
                    out_r2_f.write(rec2)
                    written += 1

        print(f"[INFO] Processed paired records: {processed}")
        print(f"[INFO] Wrote paired records: {written}")
        if written == 0:
            print(f"[WARN] No reads written for sample {self.sample_id}. Check read-id file and FASTQ naming.")
        else:
            print(f"[OK] Output FASTQs (paired) written: {out_r1_path}, {out_r2_path} ({written} reads each)")

    def build_statement(self) -> str:
        # JSON-serialize PARAMS (JSON uses double-quotes so safe to embed in single-quoted Python literal)
        params_json = json.dumps(self.PARAMS)

        # Build a single-line python snippet. Use repr() for the file paths so they become valid Python literals.
        py_code = (
            "import json; "
            "from ocmsshotgun.modules.MetaRefinement import FilterFastqByIds; "
            f"tool = FilterFastqByIds({repr(str(self.read_id_file))}, "
            f"[{repr(str(self.out_r1))}, {repr(str(self.out_r2))}], "
            f"**json.loads({repr(params_json)})); "
            "tool.run()"
        )

        # Shell-quote the whole python snippet so it is passed as a single safe argument to the shell
        return "python3 -c " + shlex.quote(py_code)
# ----------------------------
# Module main guard
# ----------------------------
if __name__ == "__main__":
    print("MetaRefinement module.")
    print("Import ExtractRefinedBinReads, FilterFastqByIds")

