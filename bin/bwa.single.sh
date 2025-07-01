#!/bin/bash

# Function to display help information
display_help() {
    echo "Usage: bash $0 reference.fasta input.fastq output_dir"
    echo
    echo "This script aligns reads to a reference genome using BWA and processes the results."
    echo
    echo "Arguments:"
    echo "  reference.fasta  Path to the reference genome file."
    echo "  input.fastq      Path to the input FASTQ file."
    echo "  output_dir       Directory to store output files."
    echo
    echo "Example:"
    echo "  bash $0 reference.fasta input.fastq output_directory"
    echo
    exit 1
}

# Check if all required arguments are provided
if [[ $# -lt 3 ]]; then
    display_help
fi

# Parse command line arguments
reference="$1"
input_fastq="$2"
output_dir="$3"
threads="$4"

echo "Running bwa.single.sh with the following parameters:"
echo "Reference: $reference"
echo "Input FASTQ: $input_fastq"
echo "Output Directory: $output_dir"

# Create output directory if it does not exist
echo "Creating output directory (if it doesn't exist): $output_dir"
mkdir -p "$output_dir"


# 检查输出目录中是否存在任何指定文件
files_exist=true
for file in coverage.txt repeat-spanning_read.txt spanning_reads_sorted.bam spanning_reads_sorted.bam.bai spanning_reads_sequencing_depth.pdf; do
    if [ ! -f "$output_dir/$file" ]; then
        files_exist=false
        break
    fi
done

# 如果文件不存在，则运行 minimap2 命令
if [ "$files_exist" = false ]; then
    echo "Aligning reads to reference genome..."
    echo "Running bwa."
    # Index reference genome
    echo "Indexing reference genome..."
    bwa index "$reference"
    bwa mem -t "$threads" "$reference" "$input_fastq" > "$output_dir/output.sam"

    # Remove BWA index files
    [ -f "${reference}.amb" ] && rm "${reference}.amb"
    [ -f "${reference}.ann" ] && rm "${reference}.ann"
    [ -f "${reference}.bwt" ] && rm "${reference}.bwt"
    [ -f "${reference}.pac" ] && rm "${reference}.pac"
    [ -f "${reference}.sa" ] && rm "${reference}.sa"
fi


