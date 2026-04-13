import os
import glob

class DepthFileManager:
    def __init__(self, depth_dir="01_mapping.dir"):
        self.depth_dir = depth_dir

    def get_depth_files(self, mapping_mode, tool):
        """
        Get depth files for the given binning tool.
        
        Parameters
        ----------
        mapping_mode : str
            "many2one" or "one2one"
        tool : str
            Which tool's depth files to fetch ("metabat2" or "maxbin2")
        """
        if tool not in ("metabat2", "maxbin2"):
            raise ValueError(f"Unsupported tool: {tool}")

        if mapping_mode == "many2one":
            pattern = os.path.join(self.depth_dir, f"cumulative_{tool}_depth_*.txt")
            depth_files = sorted(glob.glob(pattern))
        
            if not depth_files:
                single = os.path.join(self.depth_dir, f"cumulative_{tool}_depth.txt")
                if os.path.exists(single):
                    depth_files = [single]
        
            if not depth_files:
                raise FileNotFoundError(f"No pooled cumulative {tool} depth files found in {self.depth_dir}")

            return depth_files

        elif mapping_mode == "one2one":
            pattern = os.path.join(self.depth_dir, f"*_{tool}_depth.txt")
            depth_files = glob.glob(pattern)
            
            # exclude cumulative marker files
            depth_files = [
                f for f in depth_files
                if "cumulative" not in os.path.basename(f)
            ]
            if not depth_files:
                raise FileNotFoundError(f"No per-sample {tool} depth files in {self.depth_dir}")
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
            name = os.path.basename(depth_file)
            if "_" in name:
                facility = os.path.splitext(name)[0].rsplit("_", 1)[-1]
            else:
                facility = os.path.splitext(name)[0]

            candidate = os.path.join(fasta_dir, f"{facility}_pooled.fasta")
            if os.path.exists(candidate):
                fasta_path = candidate
                prefix = os.path.splitext(os.path.basename(candidate))[0]
            else:
                fasta_files = glob.glob(os.path.join(fasta_dir, "*.fasta"))
                if len(fasta_files) != 1:
                    raise FileNotFoundError(
                        f"Expected 1 pooled FASTA in {fasta_dir}, found: {fasta_files}"
                    )
                fasta_path = fasta_files[0]
                prefix = os.path.splitext(os.path.basename(fasta_path))[0]

        else:
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
    def run_all(output_dir, PARAMS, tool="metabat2"):
        """
        Build binning commands for all samples (one2one) or pooled run (many2one).
        Returns a list of (prefix, command).
        """
        mapping_mode = PARAMS["mapfastq2fasta"]["mapping_mode"]
        depth_manager = DepthFileManager("01_mapping.dir")
        depth_files = depth_manager.get_depth_files(mapping_mode, tool=tool)

        if depth_files is None:
            raise RuntimeError(f"DepthFileManager.get_depth_files returned None for mode={mapping_mode}, tool={tool}")
        if not depth_files:
            raise FileNotFoundError(f"No depth files found for mode={mapping_mode}, tool={tool} in 01_mapping.dir")

        commands = []
        if mapping_mode == "many2one":
            for depth_file in depth_files:
                base = os.path.basename(depth_file)
                prefix = base.replace(f"cumulative_{tool}_depth_", "").replace(".txt", "")
                if not prefix:
                    prefix = "pooled"

                sample_dir = os.path.join(output_dir, prefix, f"{tool}_bins")
                os.makedirs(sample_dir, exist_ok=True)

                cmd = MetaBAT2Runner.run_for_sample(depth_file, sample_dir, PARAMS, mapping_mode)
                commands.append((prefix, cmd))

        else:  # one2one
            for depth_file in depth_files:
                prefix = os.path.basename(depth_file).replace(f"_{tool}_depth.txt", "")
                sample_dir = os.path.join(output_dir, prefix, f"{tool}_bins")
                os.makedirs(sample_dir, exist_ok=True)
                cmd = MetaBAT2Runner.run_for_sample(depth_file, sample_dir, PARAMS, mapping_mode)
                commands.append((prefix, cmd))

        return commands


class MaxBin2Runner:
    """
    Wraps MaxBin2 command building for pooled (many2one) or per-sample (one2one) binning.
    """

    @staticmethod
    def run_for_sample(depth_file, output_dir, PARAMS, mapping_mode):
        fasta_dir = PARAMS["general"]["fasta_dir"]

        if mapping_mode == "many2one":
            name = os.path.basename(depth_file)
            if "_" in name:
                facility = os.path.splitext(name)[0].rsplit("_", 1)[-1]
            else:
                facility = os.path.splitext(name)[0]

            candidate = os.path.join(fasta_dir, f"{facility}_pooled.fasta")
            if os.path.exists(candidate):
                fasta_path = candidate
                prefix = os.path.splitext(os.path.basename(candidate))[0]
            else:
                fasta_files = glob.glob(os.path.join(fasta_dir, "*.fasta"))
                if len(fasta_files) != 1:
                    raise FileNotFoundError(
                        f"Expected 1 pooled FASTA in {fasta_dir}, found: {fasta_files}"
                    )
                fasta_path = fasta_files[0]
                prefix = os.path.splitext(os.path.basename(fasta_path))[0]

        else:
            prefix = os.path.basename(depth_file).replace("_maxbin2_depth.txt", "")
            fasta_path = os.path.join(fasta_dir, f"{prefix}.fasta")
            if not os.path.exists(fasta_path):
                raise FileNotFoundError(
                    f"FASTA not found for {depth_file}: {fasta_path}"
                )

        # Params
        threads = PARAMS["maxbin2"].get("threads", 4)

        os.makedirs(output_dir, exist_ok=True)
        out_prefix = os.path.join(output_dir, f"{prefix}_bin")
        log_file = os.path.join(output_dir, f"{prefix}_maxbin2.log")

        # Build command
        return (
            f"run_MaxBin.pl -contig {fasta_path} "
            f"-abund {depth_file} "
            f"-out {out_prefix} "
            f"-thread {threads} "
            f"> {log_file} 2>&1 && gzip -f {out_prefix}.*.fasta"
        )

    @staticmethod
    def run_all(output_dir, PARAMS, tool="maxbin2"):
        mapping_mode = PARAMS["mapfastq2fasta"]["mapping_mode"]
        depth_manager = DepthFileManager("01_mapping.dir")
        depth_files = depth_manager.get_depth_files(mapping_mode, tool=tool)

        if depth_files is None:
            raise RuntimeError(f"DepthFileManager.get_depth_files returned None for mode={mapping_mode}, tool={tool}")
        if not depth_files:
            raise FileNotFoundError(f"No depth files found for mode={mapping_mode}, tool={tool} in 01_mapping.dir")

        commands = []
        if mapping_mode == "many2one":
            for depth_file in depth_files:
                base = os.path.basename(depth_file)
                prefix = base.replace(f"cumulative_{tool}_depth_", "").replace(".txt", "")
                if not prefix:
                    prefix = "pooled"

                sample_dir = os.path.join(output_dir, prefix, f"{tool}_bins")
                os.makedirs(sample_dir, exist_ok=True)

                cmd = MaxBin2Runner.run_for_sample(depth_file, sample_dir, PARAMS, mapping_mode)
                commands.append((prefix, cmd))
        else:
            for depth_file in depth_files:
                prefix = os.path.basename(depth_file).replace(f"_{tool}_depth.txt", "")
                sample_dir = os.path.join(output_dir, prefix, f"{tool}_bins")
                os.makedirs(sample_dir, exist_ok=True)

                cmd = MaxBin2Runner.run_for_sample(depth_file, sample_dir, PARAMS, mapping_mode)
                commands.append((prefix, cmd))

        return commands

