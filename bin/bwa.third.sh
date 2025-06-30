#!/bin/bash

# Function to display help information
display_help() {
    echo "Usage: bash $0 reference.fasta input.fastq output_dir thread seqdepth_type"
    echo
    echo "This script aligns reads to a reference genome using BWA."
    echo
    echo "Arguments:"
    echo "  reference.fasta  Path to the reference genome file."
    echo "  input.fastq      Path to the input FASTQ file."
    echo "  output_dir       Directory to store output files."
    echo "  thread           Number of threads to use."
    echo "  seqdepth_type    Sequencing depth type: 'ont' or 'pacbio'."
    echo
    exit 1
}

# Check if all required arguments are provided
if [[ $# -lt 5 ]]; then
    display_help
fi

# Parse command line arguments
reference="$1"
input_fastq="$2"
output_dir="$3"
thread="$4"
seqdepth_type="$5"

# Create output directory if it does not exist
mkdir -p "$output_dir"

# Index reference genome
#bwa index "$reference"

# Determine the correct preset for the sequencing depth type
preset_option=""
if [ "$seqdepth_type" == "ont" ]; then
    preset_option="ont2d"
elif [ "$seqdepth_type" == "pacbio" ]; then
    preset_option="pacbio"
else
    echo "Error: Invalid sequencing depth type. Please specify 'ont' or 'pacbio'."
    exit 1
fi


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
    # Index reference genome
    bwa index "$reference"
    bwa mem -t "$thread" -x "$preset_option" "$reference" "$input_fastq" > "$output_dir/output.sam"

    # Remove BWA index files
    [ -f "${reference}.amb" ] && rm "${reference}.amb"
    [ -f "${reference}.ann" ] && rm "${reference}.ann"
    [ -f "${reference}.bwt" ] && rm "${reference}.bwt"
    [ -f "${reference}.pac" ] && rm "${reference}.pac"
    [ -f "${reference}.sa" ] && rm "${reference}.sa"
fi








