#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DNA Segment Annotation Tool (version=1.0)

Description:
Annotates DNA segments using clinical genetics classification standards:
  - exonic:    Overlaps protein-coding exons (highest priority)
  - intronic:  Within gene but not in exons
  - genic:     Fully within gene boundaries (but not intronic or exonic)
  - partial:   Partially overlaps gene region
  - intergenic: Between annotated genes
  - terminal:  Before first/after last gene (linear genomes only)

Output format:
  start_pos  end_pos  annotation_type  gene_names
"""

import argparse
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import chardet
import codecs
import os, sys, re
import time
from datetime import datetime
import logging
import traceback 

################################################################################
# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)

################################################################################
def check_file_exists(file_path):
    """Check if a file exists at the specified path."""
    if not os.path.exists(file_path):
        print(f"\nError: The specified file does not exist at path: {os.path.abspath(file_path)}")
        print("Please verify:")
        print("  1. The file name is correct")
        print("  2. The file is in the expected directory")
        print("  3. You have proper read permissions for the file")
        sys.exit(1)

################################################################################
def convert_file_to_utf8(file_path):
    """Convert file encoding to UTF-8 with proper error handling."""
    # [原有代码保持不变...]

################################################################################
def parse_location(location_str):
    """Parse GenBank location string with strand awareness."""
    # [原有代码保持不变...]

################################################################################
def get_gene_name(feature):
    """Safely extract gene name from feature qualifiers."""
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    elif 'gene_synonym' in feature.qualifiers:
        return feature.qualifiers['gene_synonym'][0]
    elif 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    elif 'product' in feature.qualifiers:
        return feature.qualifiers['product'][0]
    return 'unnamed'

################################################################################
def load_genes_and_exons(genome_record):
    """Load genes and exons with improved feature handling."""
    genes = []
    exons = []
    noncoding_genes = set()
    rna_types = {'rRNA', 'tRNA', 'ncRNA', 'misc_RNA'}

    for feature in genome_record.features:
        try:
            gene_name = get_gene_name(feature)
            
            # Handle gene features
            if feature.type == 'gene':
                locations = []
                if hasattr(feature.location, 'parts'):
                    for part in feature.location.parts:
                        locations.append({
                            'start': min(part.start, part.end) + 1,
                            'end': max(part.start, part.end),
                            'strand': part.strand
                        })
                else:
                    locations.append({
                        'start': min(feature.location.start, feature.location.end) + 1,
                        'end': max(feature.location.start, feature.location.end),
                        'strand': feature.location.strand
                    })
                
                genes.append({
                    'name': gene_name,
                    'locations': locations,
                    'strand': feature.location.strand
                })
            
            # Mark non-coding RNA genes
            elif feature.type in rna_types:
                noncoding_genes.add(gene_name)
            
            # Handle CDS and exon features
            elif feature.type in ['CDS', 'exon']:
                if hasattr(feature.location, 'parts'):
                    for part in feature.location.parts:
                        exons.append({
                            'gene': gene_name,
                            'start': min(part.start, part.end) + 1,
                            'end': max(part.start, part.end),
                            'strand': part.strand
                        })
                else:
                    exons.append({
                        'gene': gene_name,
                        'start': min(feature.location.start, feature.location.end) + 1,
                        'end': max(feature.location.start, feature.location.end),
                        'strand': feature.location.strand
                    })
        
        except Exception as e:
            logger.warning(f"Skipping malformed feature: {str(e)}")
            if logger.getEffectiveLevel() == logging.DEBUG:
                logger.debug(f"Feature causing error: {feature}")

    return genes, exons, noncoding_genes

################################################################################
def annotate_segment(segment, genes, exons, noncoding_genes, min_intergenic=50):
    """Annotate DNA segments with robust feature handling."""
    seg_start, seg_end = sorted(segment)
    detailed_annotations = []
    overlapping_genes = set()
    
    # Get all gene boundaries
    gene_boundaries = []
    for gene in genes:
        for loc in gene['locations']:
            gene_boundaries.append((loc['start'], loc['end'], gene['name']))
    gene_boundaries.sort()
    
    # Check terminal regions
    if gene_boundaries and seg_end < gene_boundaries[0][0]:
        detailed_annotations.append(f"before({gene_boundaries[0][2]})[{seg_start}-{seg_end}]")
        return detailed_annotations, []
    
    if gene_boundaries and seg_start > gene_boundaries[-1][1]:
        detailed_annotations.append(f"after({gene_boundaries[-1][2]})[{seg_start}-{seg_end}]")
        return detailed_annotations, []
    
    # Prepare all features for annotation
    all_features = []
    for gene in genes:
        for loc in gene['locations']:
            all_features.append({
                'type': 'gene',
                'name': gene['name'],
                'start': loc['start'],
                'end': loc['end']
            })
    for exon in exons:
        all_features.append({
            'type': 'exon',
            'gene': exon['gene'],
            'start': exon['start'],
            'end': exon['end']
        })
    all_features.sort(key=lambda x: x['start'])
    
    # Perform annotation
    current_pos = seg_start
    while current_pos <= seg_end:
        best_annotation = None
        best_feature = None
        
        # Find overlapping features
        for feat in all_features:
            if current_pos <= feat.get('end', 0) and current_pos >= feat.get('start', 0):
                if not best_annotation or feat.get('type') == 'exon':
                    best_annotation = feat.get('type')
                    best_feature = feat
        
        if best_feature and 'name' in best_feature:
            overlap_start = max(current_pos, best_feature.get('start', 0))
            overlap_end = min(seg_end, best_feature.get('end', 0))
            
            # Special handling for non-coding RNAs and genes without exons
            if (best_feature['name'] in noncoding_genes or 
                not any(e.get('gene') == best_feature['name'] for e in exons)):
                detailed_annotations.append(
                    f"exonic({best_feature['name']})[{overlap_start}-{overlap_end}]"
                )
                overlapping_genes.add(best_feature['name'])
            else:
                # Standard gene with exons
                is_intronic = True
                for exon in exons:
                    if (exon.get('gene') == best_feature['name'] and 
                        not (overlap_end < exon.get('start', 0) or 
                             overlap_start > exon.get('end', 0))):
                        is_intronic = False
                        break
                
                if is_intronic:
                    detailed_annotations.append(
                        f"intronic({best_feature['name']})[{overlap_start}-{overlap_end}]"
                    )
                else:
                    detailed_annotations.append(
                        f"exonic({best_feature['name']})[{overlap_start}-{overlap_end}]"
                    )
                
                overlapping_genes.add(best_feature['name'])
            
            current_pos = overlap_end + 1
        else:
            # Handle intergenic regions
            next_gene_start = min(
                (f.get('start', seg_end+1) for f in all_features if f.get('start', 0) > current_pos),
                default=seg_end + 1
            )
            inter_end = min(seg_end, next_gene_start - 1)
            
            if inter_end >= current_pos:
                left_genes = [f.get('name', '') for f in all_features 
                             if f.get('type') == 'gene' and f.get('end', 0) < current_pos]
                right_genes = [f.get('name', '') for f in all_features 
                              if f.get('type') == 'gene' and f.get('start', 0) > inter_end]
                
                left = left_genes[-1] if left_genes else "START"
                right = right_genes[0] if right_genes else "END"
                
                if (inter_end - current_pos + 1) >= min_intergenic:
                    detailed_annotations.append(
                        f"intergenic({left}-{right})[{current_pos}-{inter_end}]"
                    )
            
            current_pos = inter_end + 1
    
    return detailed_annotations, sorted(overlapping_genes)

################################################################################
def main():
    """Main execution function with enhanced error handling."""
    parser = argparse.ArgumentParser(
        description='DNA Segment Annotation Tool (version=1.2)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-g', '--genbank', required=True,
                       help='Input GenBank file (.gb, .gbk)')
    parser.add_argument('-i', '--input', required=True,
                       help='File with genomic positions (start end per line)')
    parser.add_argument('-o', '--output',
                       help='Output filename prefix')
    parser.add_argument('-d', '--debug', action='store_true',
                       help='Enable detailed debug output')
    parser.add_argument('--min-intergenic', type=int, default=1,
                       help='Minimum length to show intergenic region details')

    try:
        args = parser.parse_args()
        
        check_file_exists(args.input)
        convert_file_to_utf8(args.input)
        
        if args.debug:
            logger.setLevel(logging.DEBUG)
            logger.debug("Debug mode enabled with detailed trace")
        else:
            logger.setLevel(logging.INFO)

        # Load genome
        logger.info(f"Loading genome from: {args.genbank}")
        try:
            genome_record = SeqIO.read(args.genbank, 'genbank')
            genome_length = len(genome_record.seq)
            logger.info(f"Genome loaded. Length: {genome_length} bp")
            if args.debug:
                logger.debug(f"First 3 features: {[f.type for f in genome_record.features[:3]]}")
        except Exception as e:
            logger.error(f"GenBank parsing failed: {str(e)}")
            if args.debug:
                logger.error(traceback.format_exc())
            return 1

        # Load genes and exons with non-coding RNA info
        genes, exons, noncoding_genes = load_genes_and_exons(genome_record)
        logger.info(f"Loaded {len(genes)} genes and {len(exons)} exons")
        if args.debug:
            logger.debug(f"Non-coding RNA genes: {sorted(noncoding_genes)}")

        # Load segments
        segments = []
        logger.info(f"Loading segments from: {args.input}")
        with open(args.input, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                try:
                    parts = re.split(r'[ \t!@#$%^&*()_+=,<>?/\|;:~-]+', line.strip())
                    if len(parts) != 2:
                        logger.warning(f"Line {line_num}: Invalid column number. Skip!")
                        continue
    
                    start = int(float(parts[0]))
                    end = int(float(parts[1]))
                    segments.append((min(start, end), max(start, end)))
                except ValueError as e:
                    logger.warning(f"Line {line_num}: Invalid numerical format - {str(e)}")
                    continue

        if not segments:
            logger.error("No valid segments found in input file")
            return 1
        logger.info(f"Loaded {len(segments)} segments for annotation")

        # Perform annotation
        output_prefix = args.output or os.path.splitext(os.path.basename(args.input))[0]
        output_file = f"{output_prefix}_annotations.txt"
        
        logger.info(f"Annotating segments...")
        start_time = time.time()
        
        with open(output_file, 'w') as f:
            f.write("# DNA Segment Annotation Report (v1.2)\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Genome: {os.path.basename(args.genbank)}\n")
            f.write(f"# Genome Length: {genome_length}\n")
            f.write(f"# Genes Loaded: {len(genes)}\n")
            f.write(f"# Exons Loaded: {len(exons)}\n")
            f.write(f"# Non-coding RNAs: {len(noncoding_genes)}\n")
            f.write("#\n# Annotation Key:\n")
            f.write("#  exonic(gene)[X-Y]    - Coding or non-coding RNA region\n")
            f.write("#  intronic(gene)[X-Y]  - Within intron of coding gene\n")
            f.write("#  intergenic(A-B)[X-Y] - Between specified genes\n")
            f.write("#  before(gene)[X-Y]    - Before first gene\n")
            f.write("#  after(gene)[X-Y]     - After last gene\n#\n")
            f.write("Start\tEnd\tAnnotation_Details\tGene_List\n")
    
            for seg_start, seg_end in segments:
                try:
                    annotations, gene_names = annotate_segment(
                        (seg_start, seg_end), 
                        genes, 
                        exons,
                        noncoding_genes,
                        min_intergenic=args.min_intergenic
                    )
            
                    f.write(f"{seg_start}\t{seg_end}\t")
                    f.write(f"{'; '.join(annotations) if annotations else 'unannotated'}\t")
                    f.write(f"{','.join(gene_names) if gene_names else '.'}\n")
                except Exception as e:
                    logger.error(f"Failed to annotate segment {seg_start}-{seg_end}: {str(e)}")
                    if args.debug:
                        logger.error(traceback.format_exc())
                    f.write(f"{seg_start}\t{seg_end}\tannotation_failed\t.\n")

        elapsed = time.time() - start_time
        logger.info(f"Annotations completed in {elapsed:.2f} seconds")
        logger.info(f"Results saved to: {output_file}")
        
        return 0

    except Exception as e:
        logger.error(f"Fatal error: {str(e)}")
        if args.debug:
            logger.error(traceback.format_exc())
        return 1

################################################################################
if __name__ == "__main__":
    sys.exit(main())
    