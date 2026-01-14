# BinAssembly.py
import os
import errno
from typing import List
import ocmsshotgun.modules.MetaAssembly as MA


def _safe_symlink(src: str, dst: str) -> None:
    """Create dst -> abs(src). If dst exists (file or symlink), remove it first."""
    src_abs = os.path.abspath(src)
    try:
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(src_abs, dst)
    except OSError as e:
        # race-safe recreate if exists concurrently
        if e.errno == errno.EEXIST:
            try:
                os.remove(dst)
            except OSError:
                pass
            os.symlink(src_abs, dst)
        else:
            raise


class runBinSpades(MA.MetaAssembler):
    """
    SPAdes wrapper:
      - create stable temp symlinks {bin}_R1.fastq.gz / _R2 / _R3
      - assembler() returns spades command only
      - postProcess() returns a shell command that (at runtime)
        checks contigs.fasta at top-level of {bin}.spades, creates top-level symlinks,
        removes temp links, and exits non-zero if contigs missing.
    """

    def _fetch_run_statement(self, lib_args: str, spades_out: str, **PARAMS) -> str:
        # threads/memory can be passed in PARAMS (fallbacks provided)
        threads = PARAMS.get("spades_threads", PARAMS.get("threads", 8))
        memory = PARAMS.get("spades_memory", PARAMS.get("memory", 32))
        # accept "32G" or int GB
        if isinstance(memory, str) and memory.upper().endswith("G"):
            mem_gb = int(memory[:-1])
        else:
            mem_gb = int(memory)
        return f"spades.py {lib_args} --only-assembler --threads {threads} --memory {mem_gb} -o {spades_out}"

    def assembler(self, infiles: List[str], outfile: str, **PARAMS) -> str:
        out_dir = os.path.dirname(outfile)
        bin_id = self.getTrack(infiles[0])
        spades_dir = os.path.join(out_dir, f"{bin_id}.spades")
        os.makedirs(spades_dir, exist_ok=True)

        def make_dst(i: int) -> str:
            return os.path.join(out_dir, f"{bin_id}_R{i}.fastq.gz")

        temp_links: List[str] = []
        if len(infiles) == 1:
            dst = make_dst(1)
            _safe_symlink(infiles[0], dst)
            temp_links.append(dst)
            lib_args = f"-s {dst}"
        elif len(infiles) == 2:
            for i, src in enumerate(infiles, 1):
                dst = make_dst(i)
                _safe_symlink(src, dst)
                temp_links.append(dst)
            lib_args = f"-1 {temp_links[0]} -2 {temp_links[1]}"
        elif len(infiles) == 3:
            for i, src in enumerate(infiles, 1):
                dst = make_dst(i)
                _safe_symlink(src, dst)
                temp_links.append(dst)
            lib_args = f"-1 {temp_links[0]} -2 {temp_links[1]} -s {temp_links[2]}"
        else:
            raise ValueError("runBinSpades: unexpected number of FASTQ inputs")

        # Build and return the SPAdes run command only. MetaAssembly.build() should
        # combine this with postProcess()'s command.
        spades_cmd = self._fetch_run_statement(lib_args, spades_dir, **PARAMS)
        return spades_cmd

    def postProcess(self, infiles: List[str], outfile: str, **PARAMS) -> str:
        """
        This shell fragment is run after the assembler command.
        It:
          - checks for contigs.fasta (non-empty) directly in spades_dir
          - creates absolute top-level symlinks binX.contigs.fasta and binX.scaffolds.fasta
          - makes the task outfile point to the contigs symlink
          - removes temp R1/R2/R3 links
          - exits non-zero if contigs missing
        """
        out_dir = os.path.dirname(outfile)
        bin_id = self.getTrack(infiles[0])
        spades_dir = os.path.join(out_dir, f"{bin_id}.spades")

        # absolute source paths (inside spades dir)
        contigs_src = os.path.abspath(os.path.join(spades_dir, "contigs.fasta"))
        scaffolds_src = os.path.abspath(os.path.join(spades_dir, "scaffolds.fasta"))

        # absolute output links one level up
        contigs_out = os.path.abspath(os.path.join(out_dir, f"{bin_id}.contigs.fasta"))
        scaffolds_out = os.path.abspath(os.path.join(out_dir, f"{bin_id}.scaffolds.fasta"))
        outfile_abs = os.path.abspath(outfile)

        # absolute temp link names
        temp_links_list = [
            os.path.abspath(os.path.join(out_dir, f"{bin_id}_R{i}.fastq.gz"))
            for i in (1, 2, 3)
        ]
        # create rm -f string (will ignore missing)
        temp_rm = " ".join(temp_links_list)
        temp_rm_cmd = f"rm -f {temp_rm}" if temp_rm else ""

        # shell to run AFTER spades finishes:
        # - require contigs_src exists & non-empty (-s)
        # - create/overwrite final symlinks using absolute paths
        # - remove temp links
        # - ensure outfile points to contigs_out
        # - otherwise fail (exit 1)
        shell = (
            f"if [ -s {contigs_src} ]; then "
            f"ln -sfn {contigs_src} {contigs_out} && "
            f"[ -s {scaffolds_src} ] && ln -sfn {scaffolds_src} {scaffolds_out} || true && "
            f"ln -sfn {contigs_out} {outfile_abs} && "
            f"{temp_rm_cmd}; "
            f"else echo 'ERROR: missing contigs.fasta in {spades_dir}' >&2; exit 1; fi"
        )

        return shell

