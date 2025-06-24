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
        metabat2_threads = self.PARAMS.get("metabat2_threads", 4)
        metabat2_min_contig_length = self.PARAMS.get("metabat2_min_contig_length", 2500)
        
        output_prefix = os.path.join(self.output_dir, f"{self.prefix}_bin")
        statement = (
            f"metabat2 "
            f"-i {self.assembly} "
            f"-a {self.depth_file} "
            f"-o {output_prefix} "
            f"-m {metabat2_min_contig_length} "
            f"-t {metabat2_threads}"
        )
        return statement

