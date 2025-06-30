#!/usr/bin/env python
"""
DNA Segment Annotation Tool

This script annotates DNA segments based on their positions relative to genes in a GenBank file.
It identifies whether segments are within genes, partially overlapping genes, or in intergenic regions.

Usage:
    python annotate_location_from_gbf.py -g <genbank_file> -i <positions_file> [-o output_prefix]

Output:
    Generates '<prefix>annotation.position.txt' containing the annotation results.
    When -o is not provided, uses input file name as prefix.
"""

import argparse
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import chardet
import codecs
import os

def detect_and_convert_to_utf8(file_path):
    """
    Detect file encoding and convert content to UTF-8
    
    Args:
        file_path (str): Path to the input file
        
    Returns:
        list: List of lines from the file in UTF-8 encoding
        
    Raises:
        UnicodeDecodeError: If the file cannot be decoded with detected or fallback encoding
    """
    # Detect encoding by reading first 10KB of file
    with open(file_path, 'rb') as f:
        raw_data = f.read(10000)
        result = chardet.detect(raw_data)
        encoding = result['encoding']
        print(f"Detected encoding: {encoding} (confidence: {result['confidence']:.2f})")

    try:
        # Try reading with detected encoding
        with codecs.open(file_path, 'r', encoding=encoding) as f:
            return f.read().splitlines()
    except UnicodeDecodeError:
        # Fallback to latin-1 if detection fails
        print("Falling back to latin-1 encoding")
        with codecs.open(file_path, 'r', encoding='latin-1') as f:
            return f.read().splitlines()

def parse_location(location_str):
    """
    Parse BioPython location string into start/end positions
    
    Args:
        location_str (str): Location string from GenBank feature
        
    Returns:
        list: List of (start, end) position tuples
    """
    location_str = location_str.strip()
    positions = []

    # Handle complex locations with joins/complements
    location_parts = location_str.split('(')
    for part in location_parts:
        sub_loc_str = part.split(')')[0].strip() if ')' in part else part.strip()
        
        # Extract all numbers from location string
        import re
        numbers = re.findall(r'\d+', sub_loc_str)
        if len(numbers) == 2:
            positions.append((int(numbers[0]), int(numbers[1])))
        else:
            positions.append((None, None))

    return positions

def annotate_DNA_segment(segment_start, segment_end, gene_positions):
    """
    Annotate a DNA segment based on gene positions
    
    Args:
        segment_start (int): Start position of DNA segment
        segment_end (int): End position of DNA segment
        gene_positions (list): List of (gene_name, start, end) tuples
        
    Returns:
        list: List of annotation strings for the segment
    """
    annotations = []
    segment_start, segment_end = sorted([segment_start, segment_end])
    
    # Create FeatureLocation object for easy position comparison
    segment_location = FeatureLocation(start=segment_start, end=segment_end)

    if not gene_positions:
        return ["No genes found in genome"]

    # Check if segment is before first gene
    first_gene_name, first_start, first_end = gene_positions[0]
    if segment_end < first_start:
        return [f"Before({first_gene_name})"]

    # Check if segment is after last gene
    last_gene_name, last_start, last_end = gene_positions[-1]
    if segment_start > last_end:
        return [f"After({last_gene_name})"]

    # Check each gene and intergenic region
    for i, (gene_name, gene_start, gene_end) in enumerate(gene_positions):
        # Check if segment is completely within current gene
        if segment_start >= gene_start and segment_end <= gene_end:
            annotations.append(f"Within({gene_name})")
            continue
        
        # Check for partial overlaps with current gene
        if (segment_start < gene_end and segment_end > gene_start):
            if segment_start < gene_start:
                annotations.append(f"Partial({gene_name})")
            elif segment_end > gene_end:
                annotations.append(f"Partial({gene_name})")

        # Check intergenic regions between current and next gene
        if i < len(gene_positions) - 1:
            next_gene_name, next_start, next_end = gene_positions[i+1]
            inter_start = gene_end + 1
            inter_end = next_start - 1
            
            # Check if segment is completely within intergenic region
            if segment_start >= inter_start and segment_end <= inter_end:
                annotations.append(f"Within({gene_name},{next_gene_name})")
            # Check for partial overlaps with intergenic region
            elif segment_start < inter_end and segment_end > inter_start:
                annotations.append(f"Partial({gene_name},{next_gene_name})")

    return annotations if annotations else ["No annotations found"]

def get_output_prefix(input_file, user_prefix=None):
    """
    Determine output prefix based on user input or input filename
    
    Args:
        input_file (str): Path to input file
        user_prefix (str, optional): User-specified prefix
        
    Returns:
        str: Output prefix to use
    """
    if user_prefix:
        return user_prefix
    
    # Extract base filename without extension
    base_name = os.path.basename(input_file)
    file_name, _ = os.path.splitext(base_name)
    return file_name

def main():
    # Set up command line argument parser with detailed help
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""DNA Segment Annotation Tool

Annotates DNA segments based on their positions relative to genes in a GenBank file.
Outputs annotation results to '<prefix>annotation.position.txt'.""",
        epilog="""Annotation Categories:
  Within(gene): Segment completely within a gene
  Partial(gene): Segment partially overlaps a gene
  Within(gene1,gene2): Segment in intergenic region
  Before(gene)/After(gene): Segment at genome ends

Examples:
  Basic usage:
    python annotate_location_from_gbf.py -g genome.gb -i segments.txt
  
  With custom output prefix:
    python annotate_location_from_gbf.py -g genome.gb -i segments.txt -o my_results
  
  The positions file should contain one segment per line with start and end positions:
    1000 2000
    3000 4000
    or
    1000,2000
    3000,4000""")

    # Define command line arguments
    parser.add_argument('-g', '--genbank', required=True,
                        help='GenBank format file containing gene annotations')
    parser.add_argument('-i', '--input', required=True,
                        help='File containing DNA segment positions (one segment per line)')
    parser.add_argument('-o', '--output', required=False,
                        help='Prefix for output filename (default: input filename)')
    
    args = parser.parse_args()

    # Read and parse GenBank file
    print(f"Reading GenBank file: {args.genbank}")
    try:
        genome_record = SeqIO.read(args.genbank, "genbank")
    except Exception as e:
        print(f"Error reading GenBank file: {e}")
        return

    # Extract and sort gene positions
    gene_positions = []
    for feature in genome_record.features:
        if feature.type == "gene":
            gene_name = feature.qualifiers.get('gene', ['Unknown'])[0]
            locations = parse_location(str(feature.location))
            for loc_start, loc_end in locations:
                if loc_start is not None and loc_end is not None:
                    gene_positions.append((gene_name, loc_start, loc_end))
    
    # Sort genes by genomic position
    gene_positions.sort(key=lambda x: (x[1], x[2]))
    print(f"Found {len(gene_positions)} repeats in GenBank file")

    # Read and process positions file
    print(f"Reading segment positions from: {args.input}")
    try:
        lines = detect_and_convert_to_utf8(args.input)
        segment_positions = []
        for line_num, line in enumerate(lines, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Handle various delimiters (space, comma, tab)
            parts = line.replace(',', ' ').split()
            if len(parts) >= 2:
                try:
                    # Convert to float first, then to int (handles numbers with .00)
                    start = int(float(parts[0]))
                    end = int(float(parts[1]))
                    segment_positions.append((start, end))
                except ValueError as e:
                    print(f"Warning: Invalid numbers in line {line_num}: {line}")
            else:
                print(f"Warning: Skipping malformed line {line_num}: {line}")
    except Exception as e:
        print(f"Error reading positions file: {e}")
        return

    print(f"Processing {len(segment_positions)} DNA segments")

    # Determine output filename
    output_prefix = get_output_prefix(args.input, args.output)
    output_file = f"{output_prefix}.pos.anno.txt"
    print(f"Writing annotations to: {output_file}")
    
    with open(output_file, 'w', encoding='utf-8') as f:
        # Write header with column descriptions
        f.write("# DNA Segment Position Annotations\n")
        f.write("# Columns: Segment_Start Segment_End Annotations\n")
        f.write("#\n")
        
        for start, end in segment_positions:
            annotations = annotate_DNA_segment(start, end, gene_positions)
            f.write(f"{start}\t{end}\t{', '.join(annotations)}\n")

    print("Annotation completed successfully!\n")
    #print("Note: Please verify segments marked 'Partial' or with no annotations.")
    #print("      These may require manual inspection.\n")

if __name__ == "__main__":
    main()
    