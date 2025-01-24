#!/bin/bash

# Function to display help information
display_help() {
    echo "Usage: $0 <reference> <long_reads> <output_prefix> <threads> <seqdepth_type>"
    echo
    echo "This script performs read mapping using minimap2."
    echo
    echo "Arguments:"
    echo "  reference       Path to the reference genome file."
    echo "  long_reads      Path to the long reads file."
    echo "  output_prefix   Prefix for the output directory to store results."
    echo "  threads         Number of threads to use for minimap2."
    echo "  seqdepth_type   Type of sequencing ('ont' for Oxford Nanopore or 'pacbio' for Pacific Biosciences)."
    echo
    exit 1
}

# Check for the '-h' or '--help' option or if the number of arguments is incorrect
if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -ne 5 ]; then
    display_help
fi

# Read command line arguments
reference="$1"
long_reads="$2"
output_prefix="$3"
threads="$4"
seqdepth_type="$5"

# Determine the minimap2 option based on seqdepth_type
minimap2_option=""
if [ "$seqdepth_type" == "ont" ]; then
    minimap2_option="map-ont"
elif [ "$seqdepth_type" == "pacbio" ]; then
    minimap2_option="map-pb"
else
    echo "Invalid seqdepth_type. Please specify 'ont' or 'pacbio'."
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_prefix"
echo "Storing results in directory '$output_prefix'."

# Run minimap2 to perform mapping
#echo "Running minimap2 with $seqdepth_type sequencing..."
#minimap2 -ax "$minimap2_option" -t "$threads" "$reference" "$long_reads" > "$output_prefix/output.sam"


# 检查输出目录中是否存在任何指定文件
# 初始化变量  
all_files_exist=true  
  
# 循环遍历文件列表  
for file in coverage.txt repeat-spanning_read.txt spanning_reads_sorted.bam spanning_reads_sorted.bam.bai spanning_reads_sequencing_depth.pdf; do  
    if [ ! -f "$output_prefix/$file" ]; then  
        all_files_exist=false  
        break  # 如果找到一个不存在的文件就跳出循环  
    fi  
done  
  
# 如果有任何文件不存在，则运行 minimap2 命令  
if [ "$all_files_exist" = false ]; then  
    echo "At least one file is missing, running minimap2."  
    minimap2 -ax "$minimap2_option" -t "$threads" "$reference" "$long_reads" > "$output_prefix/output.sam"  
fi


