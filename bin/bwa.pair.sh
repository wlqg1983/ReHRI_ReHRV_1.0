#!/bin/bash

# Function to display help information
display_help() {
    echo "Usage: $0 <reference> <reads1> <reads2> <output_prefix>"
    echo
    echo "This script performs read mapping using minimap2 for paired-end reads and saves the results."
    echo
    echo "Arguments:"
    echo "  reference      Path to the reference genome file."
    echo "  reads1         Path to the first set of reads."
    echo "  reads2         Path to the second set of reads."
    echo "  output_prefix  Prefix for the output directory."
    echo
    echo "Example:"
    echo "  $0 reference.fasta reads1.fastq reads2.fastq output_directory"
    echo
    exit 1
}

# Check for the '-h' or '--help' option or if the number of arguments is incorrect
if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -ne 5 ]; then
    display_help
fi

# Read command line arguments
reference="$1"
reads1="$2"
reads2="$3"
output_dir="$4"
threads="$5"

# Create the output directory
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
    echo "Running bwa."
    bwa index "$reference"
    bwa mem -t "$threads" "$reference" "$reads1" "$reads2" > "$output_dir/output.sam"

    # Remove BWA index files
    [ -f "${reference}.amb" ] && rm "${reference}.amb"
    [ -f "${reference}.ann" ] && rm "${reference}.ann"
    [ -f "${reference}.bwt" ] && rm "${reference}.bwt"
    [ -f "${reference}.pac" ] && rm "${reference}.pac"
    [ -f "${reference}.sa" ] && rm "${reference}.sa"
fi


