#!/bin/bash

# Function to display help information
display_help() {
    echo "Usage: $0 <reference> <reads> <output_dir>"
    echo
    echo "This script performs read mapping using minimap2 and saves the results in a specified output directory."
    echo
    echo "Arguments:"
    echo "  reference    Path to the reference genome file."
    echo "  reads        Path to the reads file."
    echo "  output_dir   Directory to store output files."
    echo
    echo "Example:"
    echo "  $0 reference.fasta reads.fastq output_directory"
    echo
    exit 1
}

# Check for the '-h' or '--help' option or if the number of arguments is incorrect
if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -ne 4 ]; then
    display_help
fi

# Read command line arguments
reference="$1"
reads="$2"
output_dir="$3"
thread="$4"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"
#echo "Results will be stored in '$output_dir'."

# Run minimap2 to perform mapping
#minimap2 -ax sr -t "$thread" "$reference" "$reads" > "$output_dir/output.sam"


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
    echo "Running minimap2."
    minimap2 -ax sr -t "$thread" "$reference" "$reads" > "$output_dir/output.sam"
fi
