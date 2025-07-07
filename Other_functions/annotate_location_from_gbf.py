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
import time  # 用于性能计时
from datetime import datetime  # 用于生成时间戳
import logging
import traceback 

################################################################################
# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)

################################################################################################################################################################
def check_file_exists(file_path):
    """
    Check if a file exists at the specified path.
    
    Args:
        file_path (str): Path to the file to be checked
        
    Returns:
        None: Exits program with error code 1 if file doesn't exist
    """
    # Check if file exists before proceeding
    if not os.path.exists(file_path):
        print(f"\nError: The specified file does not exist at path: {os.path.abspath(file_path)}")
        print("Please verify:")
        print("  1. The file name is correct")
        print("  2. The file is in the expected directory")
        print("  3. You have proper read permissions for the file")
        sys.exit(1)  # Exit with error code 1

################################################################################
def convert_file_to_utf8(file_path):
    """
    Convert file encoding to UTF-8 (no return value, directly modifies the source file)
    Automatically handles BOM headers, line endings, and encoding exceptions
    
    Args:
        file_path (str): Path to the file to be converted
        
    Raises:
        UnicodeDecodeError: When the file cannot be decoded with any encoding
        IOError: When file read/write errors occur
    """
    # Phase 1: Detect file encoding
    with open(file_path, 'rb') as f:
        raw_data = f.read(20480)
        result = chardet.detect(raw_data)
        original_encoding = result['encoding']
        #print(f"Detected original encoding: {original_encoding} (confidence: {result['confidence']:.2f})")

    # Phase 2: Read file content (try multiple encodings)
    content = None
    encodings_to_try = [
        'utf-8-sig',  # Try UTF-8 with BOM first
        original_encoding,  # Then try the detected encoding
        'gb18030',     # Common encoding in Chinese environments
        'latin-1'      # Final fallback
    ]
    
    for enc in encodings_to_try:
        try:
            with open(file_path, 'r', encoding=enc) as f:
                content = f.read()
            #print(f"Successfully read file with encoding [{enc}]")
            original_encoding = enc
            break
        except UnicodeDecodeError:
            continue
    
    if content is None:
        raise UnicodeDecodeError(f"Unable to decode file: {file_path}")

    # Phase 3: Convert non-UTF-8 files
    if original_encoding.lower() not in ['utf-8', 'utf8', 'utf-8-sig']:
        print("\nWarning: Detected non-UTF-8 encoded file!")
        print(f"File path: {os.path.abspath(file_path)}")
        print(f"Current encoding: {original_encoding}")
        print("This operation will:")
        print("1. Permanently convert the file encoding to UTF-8")
        print("2. Overwrite the original file")
        print("\nPlease ensure you have backed up the original file!")
        
        # 10-second countdown confirmation
        try:
            for i in range(10, 0, -1):
                print(f"\rAuto-proceeding in {i} seconds... (Press Ctrl+C to cancel)", end='')
                sys.stdout.flush()
                time.sleep(10)
            print("\nStarting conversion...")
        except KeyboardInterrupt:
            print("\nOperation cancelled by user")
            return

        # Phase 4: Safely write to a temporary file
        temp_file = file_path + '.tmp'
        try:
            # Normalize line endings to Linux format (\n)
            normalized_content = content.replace('\r\n', '\n').replace('\r', '\n')
            
            with open(temp_file, 'w', encoding='utf-8', newline='\n') as f:
                f.write(normalized_content)
                
            # Atomically replace the original file
            os.replace(temp_file, file_path)
            print(f"File successfully converted to UTF-8 encoding (no BOM)")
            
        except Exception as e:
            if os.path.exists(temp_file):
                os.remove(temp_file)
            raise IOError(f"File conversion failed: {str(e)}")

################################################################################
def parse_location(location_str):
    """Parse GenBank location string with strand awareness"""
    try:
        # Handle complement strand
        if 'complement(' in location_str:
            parts = re.split(r'complement\(([^)]+)\)', location_str)
            locations = []
            for part in parts[1:-1]:  # Skip empty parts
                ranges = re.findall(r'(\d+)\.\.(\d+)', part)
                for start, end in ranges:
                    locations.append((int(start), int(end), -1))  # -1 for complement strand
            return locations
        
        # Handle normal strand
        ranges = re.findall(r'(\d+)\.\.(\d+)', location_str)
        return [(int(s), int(e), 1) for s, e in ranges]  # 1 for forward strand
    
    except Exception as e:
        logger.warning(f"Location parsing error: {str(e)}")
        return [(None, None, 1)]

################################################################################
def load_genes_and_exons(genome_record):
    """Load genes and exons with strand information"""
    genes = []
    exons = []
    
    for feature in genome_record.features:
        # Extract gene features
        if feature.type == 'gene':
            gene_name = feature.qualifiers.get('gene', ['unnamed'])[0]
            # Get the actual location (handles complement and joins)
            locations = []
            if hasattr(feature.location, 'parts'):
                # Handle compound locations (joins, complements)
                for part in feature.location.parts:
                    locations.append({
                        'start': min(part.start, part.end) + 1,  # Biopython uses 0-based
                        'end': max(part.start, part.end),
                        'strand': part.strand
                    })
            else:
                # Simple location
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
        
        # Extract CDS and exon features
        elif feature.type in ['CDS', 'exon']:
            gene_name = feature.qualifiers.get('gene', ['unnamed'])[0]
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
    
    return genes, exons

################################################################################
def annotate_segment(segment, genes, exons, min_intergenic=50):
    seg_start, seg_end = sorted(segment)
    detailed_annotations = []
    overlapping_genes = set()
    
    # 获取所有基因的边界和名字
    gene_boundaries = []
    for gene in genes:
        for loc in gene['locations']:
            gene_boundaries.append((loc['start'], loc['end'], gene['name']))
    gene_boundaries.sort()
    
    # 获取第一个和最后一个基因的名字
    first_gene_name = gene_boundaries[0][2] if gene_boundaries else None
    last_gene_name = gene_boundaries[-1][2] if gene_boundaries else None
    
    # Before first gene
    if gene_boundaries and seg_end < gene_boundaries[0][0]:
        detailed_annotations.append(f"before({first_gene_name})[{seg_start}-{seg_end}]")
        return detailed_annotations, []
    
    # After last gene
    if gene_boundaries and seg_start > gene_boundaries[-1][1]:
        detailed_annotations.append(f"after({last_gene_name})[{seg_start}-{seg_end}]")
        return detailed_annotations, []
    
    # 按位置排序所有特征
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
    
    # 分段注释
    current_pos = seg_start
    while current_pos <= seg_end:
        # 查找当前区域的最佳注释
        best_annotation = None
        best_feature = None
        
        for feat in all_features:
            if current_pos <= feat['end'] and current_pos >= feat['start']:
                # 优先选择exon，其次是gene
                if not best_annotation or feat['type'] == 'exon':
                    best_annotation = feat['type']
                    best_feature = feat
        
        if best_feature:
            overlap_start = max(current_pos, best_feature['start'])
            overlap_end = min(seg_end, best_feature['end'])
            
            if best_feature['type'] == 'exon':
                detailed_annotations.append(
                    f"exonic({best_feature['gene']})[{overlap_start}-{overlap_end}]"
                )
                overlapping_genes.add(best_feature['gene'])
            else:
                # 检查是否是intronic区域（在基因内但不在任何exon中）
                is_intronic = True
                for exon in exons:
                    if exon['gene'] == best_feature['name'] and \
                       not (overlap_end < exon['start'] or overlap_start > exon['end']):
                        is_intronic = False
                        break
                
                if is_intronic:
                    detailed_annotations.append(
                        f"intronic({best_feature['name']})[{overlap_start}-{overlap_end}]"
                    )
                else:
                    detailed_annotations.append(
                        f"genic({best_feature['name']})[{overlap_start}-{overlap_end}]"
                    )
                
                overlapping_genes.add(best_feature['name'])
            
            current_pos = overlap_end + 1
        else:
            # 处理intergenic区域
            next_gene_start = min(
                (f['start'] for f in all_features if f['start'] > current_pos),
                default=seg_end + 1
            )
            inter_end = min(seg_end, next_gene_start - 1)
            
            if inter_end >= current_pos:
                # 找出相邻基因
                left_genes = [
                    f['name'] for f in all_features 
                    if f['type'] == 'gene' and f['end'] < current_pos
                ]
                right_genes = [
                    f['name'] for f in all_features 
                    if f['type'] == 'gene' and f['start'] > inter_end
                ]
                
                left = left_genes[-1] if left_genes else "START"
                right = right_genes[0] if right_genes else "END"
                
                # 只有当intergenic区域足够长时才显示
                if (inter_end - current_pos + 1) >= min_intergenic:
                    detailed_annotations.append(
                        f"intergenic({left}-{right})[{current_pos}-{inter_end}]"
                    )
                
            current_pos = inter_end + 1
    
    return detailed_annotations, sorted(overlapping_genes)

################################################################################################################################################################
def main():
    parser = argparse.ArgumentParser(
        description='DNA Segment Annotation Tool (version=1.0)',
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
        
        check_file_exists(args.input)          # 检测文件是否存在
        convert_file_to_utf8(args.input)       # 检测是否是utf-8格式，如果不是则进行格式转换，并替换源文件
        
        # 配置日志级别
        if args.debug:
            logger.setLevel(logging.DEBUG)
            logger.debug("Debug mode enabled with detailed trace")
        else:
            logger.setLevel(logging.INFO)

        # 加载GenBank文件
        logger.info(f"Loading genome from: {args.genbank}")
        try:
            genome_record = SeqIO.read(args.genbank, 'genbank')  # 修复变量名
            genome_length = len(genome_record.seq)  # 使用正确的变量名
            logger.info(f"Genome loaded. Length: {genome_length} bp")
            logger.debug(f"First 3 features: {[f.type for f in genome_record.features[:3]]}")
        except Exception as e:
            logger.error(f"GenBank parsing failed: {str(e)}")
            if args.debug:
                logger.error(traceback.format_exc())
            return 1

        # 加载基因和外显子
        genes, exons = load_genes_and_exons(genome_record)
        logger.info(f"Loaded {len(genes)} genes and {len(exons)} exons")
        
        # 调试输出基因结构
        if args.debug and genes:
            sample_gene = genes[0]
            logger.debug(f"Sample gene '{sample_gene['name']}':")
            logger.debug(f"  Locations: {sample_gene['locations']}")
            logger.debug(f"  Strand: {'+' if sample_gene['strand'] > 0 else '-'}")
            logger.debug(f"  Associated exons: {[e for e in exons if e['gene'] == sample_gene['name']]}")

        # 加载待注释片段（支持多种分隔符）
        segments = []
        logger.info(f"Loading segments from: {args.input}")
        with open(args.input, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                try:
                    # 支持多种分隔符：空格、逗号、分号、横线、下划线（排除句号避免干扰小数）
                    pattern = r'[ \t!@#$%^&*()_+=,<>?/\|;:~-]+'
                    parts = re.split(pattern, line.strip())
    
                    # 验证列数
                    if len(parts) != 2:
                        logger.warning(f"Line {line_num}: Invalid column number. Skip!")
                        continue
    
                    # 预校验非数字字符（支持科学计数法如1e5）
                    def is_valid_number(s):
                        return re.fullmatch(r'^-?\d+\.?\d*([eE][-+]?\d+)?$', s.strip()) is not None
    
                    if not (is_valid_number(parts[0]) and is_valid_number(parts[1])):
                        logger.warning(f"Line {line_num}: Location information contains non numeric characters, (start={parts[0]}, end={parts[1]}). Skip!")
                        continue
    
                    # 转换为整数并确保升序
                    start = int(float(parts[0]))  # 处理科学计数法
                    end = int(float(parts[1]))
                    segments.append((min(start, end), max(start, end)))

                except ValueError as e:
                    logger.warning(f"Line {line_num}: Invalid numerical format - {str(e)}")
                    continue

        if not segments:
            logger.error("No valid segments found in input file")
            return 1
        logger.info(f"Loaded {len(segments)} segments for annotation")

        # 执行注释
        output_prefix = args.output or os.path.splitext(os.path.basename(args.input))[0]
        output_file = f"{output_prefix}_annotations.txt"
        
        logger.info(f"Annotating segments...")
        start_time = time.time()
        
        with open(output_file, 'w') as f:
            # 写入文件头注释
            f.write("# DNA Segment Annotation Report\n")
            f.write(f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Genome: {os.path.basename(args.genbank)}\n")
            f.write(f"# Genome Length: {genome_length}\n")
            f.write(f"# Genes Loaded: {len(genes)}\n")
            f.write(f"# Exons Loaded: {len(exons)}\n")
            f.write("#\n# Annotation Key:\n")
            f.write("#  exonic(gene)[X-Y]    - Overlaps coding exon\n")
            f.write("#  intronic(gene)[X-Y]  - Within intron region\n")
            f.write("#  genic(gene)[X-Y]     - Within gene but non-coding\n")
            f.write("#  intergenic(A-B)[X-Y] - Between specified genes\n")
            f.write("#  before(gene)[X-Y]    - Before first gene\n")
            f.write("#  after(gene)[X-Y]     - After last gene\n")
            f.write("#\n")
    
            # 写入表头作为第一行数据
            f.write("Start\tEnd\tAnnotation_Details\tGene_List\n")
    
            # 处理每个片段
            for seg_start, seg_end in segments:
                annotations, gene_names = annotate_segment(
                    (seg_start, seg_end), 
                    genes, 
                    exons,
                    min_intergenic=args.min_intergenic
                )
        
                # 格式化输出行
                f.write(f"{seg_start}\t{seg_end}\t")
                f.write(f"{'; '.join(annotations) if annotations else 'unannotated'}\t")
                f.write(f"{','.join(gene_names) if gene_names else '.'}\n")

        elapsed = time.time() - start_time
        logger.info(f"Annotations completed in {elapsed:.2f} seconds")
        logger.info(f"Results saved to: {output_file}")
        
        return 0

    except Exception as e:
        logger.error(f"Fatal error: {str(e)}")
        if args.debug:
            logger.error(traceback.format_exc())
        return 1

################################################################################################################################################################
if __name__ == "__main__":
    import sys
    sys.exit(main())
