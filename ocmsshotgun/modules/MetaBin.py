import os
import glob

class DepthFileManager:
    def __init__(self, depth_dir="01_mapping.dir"):
        self.depth_dir = depth_dir

    def get_depth_files(self, mapping_mode):
        if mapping_mode == "many2one":
            return [os.path.join(self.depth_dir, "cumulative_metabat2_depth.txt")]
        elif mapping_mode == "one2one":
            depth_files = glob.glob(os.path.join(self.depth_dir, "*_metabat2_depth.txt"))
            # exclude cumulative marker files
            depth_files = [
                f for f in depth_files
                if "cumulative" not in os.path.basename(f)
            ]
            if not depth_files:
                raise FileNotFoundError(f"No per-sample depth files in {self.depth_dir}")
            return sorted(depth_files)
        else:
            raise ValueError(f"Unknown mapping mode: {mapping_mode}")


class MetaBAT2Runner:
    """
    Wraps MetaBAT2 command building for pooled (many2one) or per-sample (one2one) binning.
    """

    @staticmethod
    def run_for_sample(depth_file, output_dir, PARAMS, mapping_mode):
        fasta_dir = PARAMS["general"]["fasta_dir"]

        if mapping_mode == "many2one":
            # Expect exactly one pooled FASTA
            fasta_files = glob.glob(os.path.join(fasta_dir, "*.fasta"))
            if len(fasta_files) != 1:
                raise FileNotFoundError(
                    f"Expected 1 pooled FASTA in {fasta_dir}, found: {fasta_files}"
                )
            fasta_path = fasta_files[0]
            prefix = os.path.splitext(os.path.basename(fasta_path))[0]

        else:  # one2one
            prefix = os.path.basename(depth_file).replace("_metabat2_depth.txt", "")
            fasta_path = os.path.join(fasta_dir, f"{prefix}.fasta")
            if not os.path.exists(fasta_path):
                raise FileNotFoundError(
                    f"FASTA not found for {depth_file}: {fasta_path}"
                )

        # Build command
        threads = PARAMS["metabat2"].get("threads", 4)
        min_contig_len = PARAMS["metabat2"].get("min_contig_length", 2500)

        os.makedirs(output_dir, exist_ok=True)
        out_prefix = os.path.join(output_dir, f"{prefix}_bin")
        log_file = os.path.join(output_dir, f"{prefix}_metabat2.log")

        return (
            f"metabat2 -i {fasta_path} "
            f"-a {depth_file} "
            f"-o {out_prefix} "
            f"-m {min_contig_len} "
            f"-t {threads} "
            f"> {log_file} 2>&1 && gzip -f {out_prefix}.*.fa"
        )

    @staticmethod
    def run_all(output_dir, PARAMS):
        """
        Build MetaBAT2 commands for all samples (one2one) or pooled run (many2one).
        Returns a list of (prefix, command).
        """
        mapping_mode = PARAMS["mapfastq2fasta"]["mapping_mode"]
        depth_manager = DepthFileManager()
        depth_files = depth_manager.get_depth_files(mapping_mode)

        commands = []
        if mapping_mode == "many2one":
            prefix = "pooled"
            sample_dir = os.path.join(output_dir, prefix, "metabat2_bins")
            os.makedirs(sample_dir, exist_ok=True)
            cmd = MetaBAT2Runner.run_for_sample(depth_files[0], sample_dir, PARAMS, mapping_mode)
            commands.append((prefix, cmd))
        else:  # one2one
            for depth_file in depth_files:
                prefix = os.path.basename(depth_file).replace("_metabat2_depth.txt", "")
                sample_dir = os.path.join(output_dir, prefix, "metabat2_bins")
                os.makedirs(sample_dir, exist_ok=True)
                cmd = MetaBAT2Runner.run_for_sample(depth_file, sample_dir, PARAMS, mapping_mode)
                commands.append((prefix, cmd))

        return commands

class MaxBin2Runner:
    def __init__(self, assembly, prefix, output_dir, abundance_file=None, abundance_list=None, **PARAMS):
        """
        Initialize MaxBin2Runner using parameters from YAML config.

        :param assembly: Path to assembly FASTA file
        :param abundance_file: Path to abundance/depth file (for unpooled samples)
        :param abundance_list: Path to abundance list file (for pooled samples)
        :param output_dir: Output directory for MaxBin2 bins
        :param prefix: Prefix for output bin files
        :param PARAMS: Additional parameter like maxbin2_threads
        """
        self.assembly = assembly
        self.abundance_file = abundance_file
        self.abundance_list = abundance_list
        self.output_dir = output_dir
        self.prefix = prefix
        self.PARAMS = PARAMS

    def build_command(self):
        """
        Build statement for MaxBin2.
        """
        maxbin2_threads = self.PARAMS.get("maxbin2_threads", 4)

        output_prefix = os.path.join(self.output_dir, f"{self.prefix}_bin")

        statement = (
            f"MaxBin "
            f"-fasta {self.assembly} "
            f"-out {output_prefix} "
            f"-thread {maxbin2_threads}"
        )
        # Add abundance file or list
        if self.abundance_list:
            statement += f" -abund_list {self.abundance_list}"
        elif self.abundance_file:
            statement += f" -abund {self.abundance_file}"

        return statement.strip()

