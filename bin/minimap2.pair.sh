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
output_prefix="$4"
thread="$5"

# Append "_results" to the output directory name and create directory
#output_dir=$output_prefix
mkdir -p "$output_prefix"
#echo "Results will be stored in '$output_dir'."

# Copy reference genome to output directory
cp "$reference" "$output_prefix/"

# Run minimap2 to perform mapping
#minimap2 -ax sr -t "$thread" "$reference" "$reads1" "$reads2" > "$output_prefix/output.sam"


# 检查输出目录中是否存在任何指定文件
files_exist=true
for file in coverage.txt repeat-spanning_read.txt spanning_reads_sorted.bam spanning_reads_sorted.bam.bai spanning_reads_sequencing_depth.pdf; do
    if [ ! -f "$output_prefix/$file" ]; then
        files_exist=false
        break
    fi
done

# 如果文件不存在，则运行 minimap2 命令
if [ "$files_exist" = false ]; then
    echo "Running minimap2."
    minimap2 -ax sr -t "$thread" "$reference" "$reads1" "$reads2" > "$output_prefix/output.sam"
fi


