#!/usr/bin/env python

import argparse
import os
import sys

################################################################################################################################################################
def read_spanning_read_name(file_path, output_file, spanning_length, read_number):
    data = []
    
    if not os.path.isfile(file_path):    # 检验文件是否存在，不存在就终止
        print(f"Error: The file {file_path} does not exist.")
        return
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        if len(lines) < 3:    # 假设一个有效的文件至少有3行
            print("Error: The input file is too short to process.")
            return

        # 使用制表符分割字段
        ref_genome_name = lines[0].strip().split("\t")[1]
        ref_genome_length = lines[1].strip().split("\t")[1]
        fully_covering_reads_info = lines[2].strip().split("\t")[1]
        
        # 对截取出来的序列中的重复序列，进行长度的筛选
        #repeat_length = ref_genome_length - trimmed_length*2 

        # 读取第四行至最后的内容
        spanning_reads = []
        for line in lines[3:]:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                _, min_spanning_length = parts
                min_spanning_length = int(min_spanning_length)
                if min_spanning_length >= spanning_length:         # 这里筛选的是 reads跨越重复序列的 核苷酸的数量
                    spanning_reads.append(min_spanning_length)

        # 检查满足条件的行数是否大于等于read_number
        if len(spanning_reads) >= read_number:                   # 这里筛选 reads 数量 和 重复序列的长度 reset_repeat_length 为一个数组列表
            data = [ref_genome_name, ref_genome_length, str(len(spanning_reads))]
            with open(output_file, 'a') as f:  # 追加模式
                f.write("\t".join(data) + '\n')
        else:
            print(f"Skipped: The number of repeat-spanning reads {len(spanning_reads)} is less than the required number of {read_number}.")

################################################################################################################################################################
def main():
    parser = argparse.ArgumentParser(description="Extract and save key information from a 'spanning_read_name.txt' file.")
    parser.add_argument("-i", "--input", required=True, help="Input file path for 'spanning_read_name.txt'.")
    parser.add_argument("-o", "--output", required=True, help="Output file path to save the results.")
    parser.add_argument("-cl", "--spanning_length", type=int, required=True, help="The length of a read that spanning repeats.")
    parser.add_argument("-rn", "--read_number", type=int, required=True, help="Number of reads spanning repeats.")

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # 在main函数中处理输出文件逻辑    count_subconfig_spanning_reads.py
    output_file = args.output
    
    # Process input file
    if os.path.isfile(args.input):
        read_spanning_read_name(args.input, output_file, args.spanning_length, args.read_number)    
    else:
        print(f"Error: The input file {args.input} does not exist.")


################################################################################################################################################################
if __name__ == "__main__":
    main()
    
    