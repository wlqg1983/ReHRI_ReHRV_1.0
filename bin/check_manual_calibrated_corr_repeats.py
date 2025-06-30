#!/usr/bin/env python

import argparse
import pandas as pd
from Bio import SeqIO
import sys, time, re

################################################################################################################################################################
def clean_chromosome_number(chromosome):
    # 去除所有不可见字符（包括空格、制表符等）
    cleaned_chromosome = re.sub(r'\s+', '', chromosome)
    return cleaned_chromosome

################################################################################
def extract_chromosome_number(chromosome):
    cleaned_chromosome = clean_chromosome_number(chromosome)
    # 检查编号是否符合chr后跟数字的格式
    match = re.search(r'^chr\d+$', cleaned_chromosome)
    if not match:
        raise ValueError("The chromosome number format should be chr1, chr2,..., chr10,..., chr100, etc.")

    # 提取编号中的数字部分，用于排序
    return int(re.search(r'\d+', cleaned_chromosome).group())

################################################################################
def validate_manual_calibrated_list(manual_calibrated_list, count):
    # 检查列数
    num_columns = manual_calibrated_list.shape[1]

    if count == 1:
        if num_columns != 10:
            print("ERROR: The format of the input repetitive sequence information is incorrect.")
            print("ATTENTION: When count=1, the correct header should be: fragment_id, length, start, end, direction, paired_id, paired_length, paired_start, paired_end, paired_direction. (TSV format)\n")
            exit(1)
        return None, None
    else:
        if num_columns != 12:
            print("ERROR: The format of the input repetitive sequence information is incorrect.")
            print("ATTENTION: When count>1, the correct header should be: fragment_id, length, start, end, direction, chromosome, paired_id, paired_length, paired_start, paired_end, paired_direction, paired_chromosome. (TSV format)\n")
            exit(1)

        # 验证主片段和配对片段的染色体编号
        def validate_chromosomes(chrom_column):
            chromosomes = manual_calibrated_list[chrom_column].unique()
            cleaned_chromosomes = [clean_chromosome_number(chrom) for chrom in chromosomes]
            for chromosome in cleaned_chromosomes:
                extract_chromosome_number(chromosome)
            return sorted(cleaned_chromosomes, key=extract_chromosome_number)

        main_chromosomes = validate_chromosomes("chromosome")
        paired_chromosomes = validate_chromosomes("paired_chromosome")
        
        return main_chromosomes, paired_chromosomes

################################################################################
def calculate_chromosome_lengths(fasta_file, count):
    seq_lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_lengths.append(len(record.seq))
    if count == 1:
        return [seq_lengths[0]]  # 仅返回一条染色体的长度
    elif count > 1:
        cumulative_lengths = [sum(seq_lengths[:i+1]) for i in range(len(seq_lengths))]
        return cumulative_lengths

################################################################################
def validate_chromosome_headers(fasta_file, chromosomes, count):
    if count == 1:
        return  # 单染色体不需要验证
    
    fasta_headers = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
    if sorted(fasta_headers) != sorted(chromosomes):
        print(f"Warning: The chromosomes in the {fasta_file} will be processed in the given order, treated as chr1, chr2, and so on.")
        time.sleep(5)  # 等待用户终止程序

################################################################################
def adjust_start_end(manual_calibrated_list, cumulative_lengths, count):
    if count == 1:
        return manual_calibrated_list  # 单染色体不需要调整
    
    # 调整主片段的start和end
    manual_calibrated_list["start"] = manual_calibrated_list.apply(
        lambda row: row["start"] + cumulative_lengths[int(row["chromosome"][3:]) - 2] 
        if int(row["chromosome"][3:]) > 1 else row["start"], axis=1)
    manual_calibrated_list["end"] = manual_calibrated_list.apply(
        lambda row: row["end"] + cumulative_lengths[int(row["chromosome"][3:]) - 2] 
        if int(row["chromosome"][3:]) > 1 else row["end"], axis=1)
    
    # 调整配对片段的start和end
    manual_calibrated_list["paired_start"] = manual_calibrated_list.apply(
        lambda row: row["paired_start"] + cumulative_lengths[int(row["paired_chromosome"][3:]) - 2] 
        if int(row["paired_chromosome"][3:]) > 1 else row["paired_start"], axis=1)
    manual_calibrated_list["paired_end"] = manual_calibrated_list.apply(
        lambda row: row["paired_end"] + cumulative_lengths[int(row["paired_chromosome"][3:]) - 2] 
        if int(row["paired_chromosome"][3:]) > 1 else row["paired_end"], axis=1)
    
    return manual_calibrated_list

################################################################################
def check_start_end_max_length(manual_calibrated_list, sorted_chromosomes, cumulative_lengths, count):
    if count == 1:
        # 单染色体只需检查是否超过总长度
        max_start = manual_calibrated_list["start"].max()
        max_paired_start = manual_calibrated_list["paired_start"].max()
        max_end = manual_calibrated_list["end"].max()
        max_paired_end = manual_calibrated_list["paired_end"].max()
        chromosome_length = cumulative_lengths[0]
        
        if (max_start > chromosome_length or max_end > chromosome_length or 
            max_paired_start > chromosome_length or max_paired_end > chromosome_length):
            raise ValueError(f"The repetitive sequence position exceeds the chromosome length.\nMaximum start: {max_start}, maximum end: {max_end}, maximum paired start: {max_paired_start}, maximum paired end: {max_paired_end}, chromosome length: {chromosome_length}.\n")
        return
    
    # 多染色体检查
    def check_fragment(df, chrom_column, start_column, end_column, fragment_type):
        for i, chromosome in enumerate(sorted_chromosomes):
            chrom_rows = df[df[chrom_column] == chromosome]
            if chrom_rows.empty:
                continue

            max_start = chrom_rows[start_column].max()
            max_end = chrom_rows[end_column].max()
            chromosome_length = cumulative_lengths[i] if i == 0 else cumulative_lengths[i] - cumulative_lengths[i - 1]

            if max_start > chromosome_length or max_end > chromosome_length:
                raise ValueError(f"The {fragment_type} repetitive sequence position of Chromosome {chromosome} exceeds the length of this chromosome.\nMaximum start: {max_start}, maximum end: {max_end}, chromosome length: {chromosome_length}.\n")
    
    # 检查主片段和配对片段
    check_fragment(manual_calibrated_list, "chromosome", "start", "end", "main")
    check_fragment(manual_calibrated_list, "paired_chromosome", "paired_start", "paired_end", "paired")

################################################################################
def check_direction(df, count):
    # 检查主片段和配对片段的direction
    for index, row in df.iterrows():
        # 主片段检查
        direction = row['direction']
        start = row['start']
        end = row['end']
        if direction == 'plus' and end <= start:
            raise ValueError(f"Error in row {index + 1}: When main fragment direction is 'plus', 'end' should be > 'start'.")
        elif direction == 'minus' and end >= start:
            raise ValueError(f"Error in row {index + 1}: When main fragment direction is 'minus', 'end' should be < 'start'.")
            
        # 配对片段检查
        paired_direction = row['paired_direction']
        paired_start = row['paired_start']
        paired_end = row['paired_end']
        if paired_direction == 'plus' and paired_end <= paired_start:
            raise ValueError(f"Error in row {index + 1}: When paired fragment direction is 'plus', 'paired_end' should be > 'paired_start'.")
        elif paired_direction == 'minus' and paired_end >= paired_start:
            raise ValueError(f"Error in row {index + 1}: When paired fragment direction is 'minus', 'paired_end' should be < 'paired_start'.")

######################################################################################################################################################################
################################################################################
def main():
    parser = argparse.ArgumentParser(description="Detect the data format of input files and process chromosome sequences")
    parser.add_argument("-f", "--fasta_file", required=True, help="FASTA file storing chromosome sequences")
    parser.add_argument("-l", "--list_file", required=True, help="manual_calibrated_list file (TSV format)")
    parser.add_argument("-c", "--count", type=int, required=True, help="Chromosome count (1 or >1)")
    
    args = parser.parse_args()
    
    if not isinstance(args.count, int) or args.count < 1:
        parser.error("The number of chromosomes must be a positive integer (1 or >1).")
        exit(1)

    # 读取 manual_calibrated_list 文件
    manual_calibrated_list = pd.read_csv(args.list_file, sep='\t')
    
    # 调用 check_direction 函数进行检查
    check_direction(manual_calibrated_list, args.count)

    # 计算染色体长度及连加矩阵
    cumulative_lengths = calculate_chromosome_lengths(args.fasta_file, args.count)

    if args.count == 1:
        # 验证列数
        main_chromosomes, paired_chromosomes = validate_manual_calibrated_list(manual_calibrated_list, args.count)
        
        # 检查位置是否超过染色体长度
        check_start_end_max_length(manual_calibrated_list, None, cumulative_lengths, args.count)
        
        # 输出结果
        manual_calibrated_list.to_csv("adjusted_manual_calibrated_list.tsv", sep='\t', index=False)
    else:
        # 验证染色体编号
        main_chromosomes, paired_chromosomes = validate_manual_calibrated_list(manual_calibrated_list, args.count)
        validate_chromosome_headers(args.fasta_file, main_chromosomes, args.count)
        
        # 排序数据框
        manual_calibrated_list["chromosome"] = pd.Categorical(manual_calibrated_list["chromosome"], categories=main_chromosomes, ordered=True)
        manual_calibrated_list = manual_calibrated_list.sort_values(by="chromosome")
        
        # 检查位置是否超过染色体长度
        check_start_end_max_length(manual_calibrated_list, main_chromosomes, cumulative_lengths, args.count)
        
        # 调整start和end
        adjusted_list = adjust_start_end(manual_calibrated_list, cumulative_lengths, args.count)
        adjusted_list.to_csv("adjusted_manual_calibrated_list.tsv", sep='\t', index=False)

if __name__ == "__main__":
    main()


