
import os

class MetaBAT2Runner:
    def __init__(self, assembly, depth_file, prefix, output_dir, **PARAMS):
        """
        Initialize MetaBAT2Runner using parameters from YAML config.
        :param assembly: Path to assembly FASTA file
        :param depth_file: Path to depth file
        :param output_dir: Output directory for MetaBAT2 bins
        :param prefix: Prefix for output bin files
        :param PARAMS: Additional parameters like metabat2_threads, metabat2_min_contig_length
        """
        self.assembly = assembly
        self.depth_file = depth_file
        self.output_dir = output_dir
        self.prefix = prefix
        self.PARAMS = PARAMS

    def build_command(self):
        """
        Build statement for MetaBAT2
        """
        metabat2_threads = self.PARAMS.get("threads", 4)
        metabat2_min_contig_length = self.PARAMS.get("min_contig_length", 2500)
        
        output_prefix = os.path.join(self.output_dir, f"{self.prefix}_bin")
        log_file = os.path.join(self.output_dir, f"{self.prefix}_metabat2.log")

        statement = (
            f"metabat2 "
            f"-i {self.assembly} "
            f"-a {self.depth_file} "
            f"-o {output_prefix} "
            f"-m {metabat2_min_contig_length} "
            f"-t {metabat2_threads} "
            f"> {log_file} 2>&1 && "  # Capture stdout+stderr into log
            f"gzip -f {output_prefix}.*.fa"  # Compress all resulting bin files 
        )
        return statement

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

