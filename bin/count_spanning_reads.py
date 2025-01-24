#!/usr/bin/env python

import argparse
import os
import sys

#########################################################################################################################################################################
def read_spanning_read_name(file_path, output_file):
    data = []
    
    if not os.path.isfile(file_path):    # 检查验repeat-spanning_read.txt是否存在, 不存在就终止
        return
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        if len(lines) < 3:    # Assuming a valid file has at least 3 lines
            return     #False    # Indicating no content generated

        ref_genome_name = lines[0].split(":")[1].strip()
        ref_genome_length = lines[1].split(":")[1].strip()
        
        fully_covering_reads_info = lines[2].split(":")[1].strip()
        fully_covering_reads_count = int(''.join(filter(str.isdigit, fully_covering_reads_info)))
        data = [ref_genome_name, ref_genome_length, str(fully_covering_reads_count)]

    with open(output_file, 'a') as f:  # Append mode
        f.write("\t".join(data) + '\n')

#########################################################################################################################################################################
def main():
    parser = argparse.ArgumentParser(description="Extract and save key information from a 'spanning_read_name.txt' file.")
    parser.add_argument("-i", "--input", required=True, help="Input file path for 'spanning_read_name.txt'.")
    parser.add_argument("-o", "--output", required=True, help="Output file path to save the results.")
    
    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # 在main函数中处理输出文件逻辑
    output_file = args.output
    
    # Process input file
    if os.path.isfile(args.input):
        read_spanning_read_name(args.input, output_file)
    
    #print("************ One repeat unit has been processed completely! ************")

#########################################################################################################################################################################
if __name__ == "__main__":
    main()


