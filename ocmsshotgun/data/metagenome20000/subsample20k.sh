#!/bin/bash

# load seqtk module
module load seqtk/1.4-GCC-12.2.0

DIR="/well/kir/collab/ocms_johnson_powrie_shared/test_data/metagenome"
# Check if the directory exists
if [ ! -d "$DIR" ]; then
    echo "Directory $DIR does not exist."
    exit 1
fi

for f1 in "$DIR"/*.fastq.1.gz; do

    # Extract the sample name by removing the suffix
    sample="${f1##*/}" # file basename
    sample="${sample%%.fastq.1.gz}" # remove suffix

    # get fastq2
    f2="${DIR}/${sample}.fastq.2.gz"
    
    # Construct the output file name
    out_f1="${sample}_short20000.fastq.1.gz"
    out_f2="${sample}_short20000.fastq.2.gz"
    
    # subsample 20k reads with seed set to 100
    seqtk sample -s100 "$f1" 20000 > "$out_f1"
    seqtk sample -s100 "$f2" 20000 > "$out_f2"
done
