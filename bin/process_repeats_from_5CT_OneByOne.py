#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main():
    # 创建解析器
    parser = argparse.ArgumentParser(description="Process repeat sequences from 5CT and extract relevant sequences.")
    
    # 添加参数
    parser.add_argument('-p', '--project_id', required=True, help='Project ID.')
    parser.add_argument('-s', '--seq_file', required=True, help='Genome sequence file.')
    parser.add_argument('-r', '--rep_5CT', required=True, help='Repeat information with 5CT.')
    parser.add_argument('-n', '--chr_num', help='Number of chromosome.')
    
    # 解析参数
    args = parser.parse_args()
    
    # 获取当前工作路径
    current_path = os.getcwd()
    
    # 读取重复序列位置信息数据
    data = pd.read_csv(args.rep_5CT, sep='\t')
    
    # 去掉含ctg的行
    data = data[~data["fragment_id"].str.contains("ctg")]
    
    # 提取数字和字母部分进行排序
    data['fragment_id_num'] = data['fragment_id'].str.extract(r'(\d+)').astype(int)
    data['fragment_id_suffix'] = data['fragment_id'].str.extract(r'(\D)$')
    
    # 按照数字和字母部分排序
    data = data.sort_values(by=['fragment_id_num', 'fragment_id_suffix'])
    
    # 去重，保留第一次出现的行
    data = data.drop_duplicates(subset='fragment_id', keep='first')
    
    # 去除添加的辅助列
    data = data.drop(columns=['fragment_id_num', 'fragment_id_suffix'])
    
    # 取前三列
    data = data.iloc[:, :3]
    
    # 计算重复序列的长度
    data["length"] = abs(data["end"] - data["start"]) + 1
    
    # 标记正负
    data["direction"] = data.apply(lambda row: "minus" if row["end"] < row["start"] else "plus", axis=1)
    
    # 调整列顺序
    data = data[["fragment_id", "length", "start", "end", "direction"]]
    
    
    # 读取基因组文件  格式为FASTA文件
    fasta_file = args.seq_file
    if os.path.isfile(fasta_file):
        # Since the file contains only one sequence, directly parse the first record
        sequence_record = next(SeqIO.parse(fasta_file, "fasta"))
        # 计算序列的长度
        sequence_length = len(sequence_record.seq)
    else:
        print(f"Genome sequence {fasta_file} file not found.")
        sys.exit(1)

    # 截取序列并保存为新的FASTA文件
    extracted_sequences = []
    record_length = 0
    for index, row in data.iterrows():
        fragment_id = row["fragment_id"]
        start = row["start"]
        end = row["end"]
    
        if start > end:
            start, end = end, start  # 先保证起点 < 终点
    
        # 调整索引以适应FASTA格式
        start = start - 1
    
        # 更新最大结束位置
        if end > record_length:
            record_length = end
    
        if sequence_length < record_length:
            print("ERROR: The repeat position exceeds the length of the genome.")
            exit(1)
    
        # 边界检查
        if start < 0 or end > sequence_length:
            print(f"ERROR: Index out of range for fragment {fragment_id}. start={start}, end={end}, sequence_length={sequence_length}")
            continue
    
        # Access the sequence from the single SeqRecord object
        sequence = sequence_record.seq[start:end]
    
        if row["direction"] == "minus":
            sequence = sequence.reverse_complement()
    
        extracted_sequences.append(SeqRecord(sequence, id=fragment_id, description=""))

    if args.chr_num:
        output_fasta_file = os.path.join(current_path, f"repeat_sequences_{args.project_id}_chr{args.chr_num}.fasta")
        SeqIO.write(extracted_sequences, output_fasta_file, "fasta")
    else:
        output_fasta_file = os.path.join(current_path, f"repeat_sequences_{args.project_id}.fasta")
        SeqIO.write(extracted_sequences, output_fasta_file, "fasta")


if __name__ == "__main__":
    main()


