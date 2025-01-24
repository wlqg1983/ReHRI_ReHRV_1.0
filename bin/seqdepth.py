#!/usr/bin/env python

import argparse
import os, sys, csv
import subprocess
import shutil, re, glob
from Bio import SeqIO
import pandas as pd

################################################################################
def draw_seq_depth_next(alignment, reference, sequencing_type, output_prefix, threads, input_fastq): 
    # Get the absolute path of the current working directory
    abs_work_dir = os.path.dirname(os.path.realpath(__file__))

    # Call the alignment script based on sequencing type and alignment type
    if sequencing_type == "NGS single end":
        if alignment == "minimap2":
            subprocess.run(["bash", os.path.join(abs_work_dir, "minimap2.single.sh"), reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_results"), threads], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True) 
        elif alignment == "bwa":
            subprocess.run(["bash", os.path.join(abs_work_dir, "bwa.single.sh"), reference, input_fastq,
                             os.path.join(os.getcwd(), output_prefix + "_results"), threads], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True) 
    elif sequencing_type == "NGS pair end":
        if alignment == "minimap2":
            subprocess.run(["bash", os.path.join(abs_work_dir, "minimap2.pair.sh"), reference, input_fastq[0], input_fastq[1],
                             os.path.join(os.getcwd(), output_prefix + "_results"), threads], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True) 
        elif alignment == "bwa":
            subprocess.run(["bash", os.path.join(abs_work_dir, "bwa.pair.sh"), reference, input_fastq[0], input_fastq[1],
                             os.path.join(os.getcwd(), output_prefix + "_results"), threads], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True) 

################################################################################
def draw_seq_depth_third(alignment, reference,output_prefix, threads, input_fastq, seqdepth_type): 
    # Get the absolute path of the current working directory
    abs_work_dir = os.path.dirname(os.path.realpath(__file__))

    # Call the alignment script based on sequencing type and alignment type
    if alignment == "minimap2":
        subprocess.run(["bash", os.path.join(abs_work_dir, "minimap2.third.sh"), reference, input_fastq,
                         os.path.join(os.getcwd(), output_prefix + "_results"), threads, seqdepth_type], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True) 
    elif alignment == "bwa":
        subprocess.run(["bash", os.path.join(abs_work_dir, "bwa.third.sh"), reference, input_fastq,
                         os.path.join(os.getcwd(), output_prefix + "_results"), threads, seqdepth_type], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True) 

################################################################################
def extract_trimmed_length_from_filename(filename):
    # 修改正则表达式以匹配新的文件名格式
    match = re.search(r'_(\d+)\.fasta$', filename)
    if match:
        return int(match.group(1))
    else:
        return None

################################################################################
def calculate_fasta_length(input_fasta):
    length = 0
    with open(input_fasta, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue  # 跳过标题行
            length += len(line.strip())  # 计算序列长度，去除空白字符
    return length

################################################################################
def draw_seq_depth_blast(reference, sequencing_type, output_prefix, threads, input_fasta, seqdepth_type, check_spanning_length, evalue): 
    # 获取输入数据文件的前缀
    ref_prefix = os.path.splitext(os.path.basename(reference))[0]
    command = f"/usr/bin/makeblastdb -in {reference} -dbtype nucl -parse_seqids -out {output_prefix}_results/{ref_prefix}"
    # 运行命令
    os.system(command)

    # 构建命令,为后面添加不同的参数预留接口。
    command = f"/usr/bin/blastn -query {input_fasta} -out {output_prefix}_results/{ref_prefix} -db {output_prefix}_results/{ref_prefix} -outfmt 6 -evalue {evalue} -num_alignments 10000 -num_threads {threads}"

    # 运行命令
    os.system(command)
    
    # 计算reference长度
    ref_len = calculate_fasta_length(reference)

    # 提取ref文件的文件名，并提取trimmed length
    ref_name = os.path.basename(reference)    # 提取截取序列文件名
    extract_trimmed_length = extract_trimmed_length_from_filename(reference)    # 提取截取序列的长度
    
    # 修正spanning_length 
    if extract_trimmed_length < check_spanning_length:
        check_spanning_length = extract_trimmed_length

    # repeat unit position，start end
    ref_start = extract_trimmed_length+1          
    ref_end = ref_len - extract_trimmed_length

    blast_result = f"{output_prefix}_results/{ref_prefix}"    # 假设你的 TSV 文件名遵循这个格式
    output_fasta = f"{output_prefix}_results/selected_reads.fasta"

    # 检查文件是否存在且非空，对blast比对的结果进行筛选
    if os.path.exists(blast_result) and os.stat(blast_result).st_size > 0:
        # 读取CSV文件的第1, 9, 10列
        df = pd.read_csv(blast_result, sep='\t', header=None, usecols=[0, 8, 9], names=['qseqid', 'sstart', 'send'])
    
        # 如果第9列（sstart）大于第10列（send），则交换两列的值
        df.loc[df['sstart'] > df['send'], ['sstart', 'send']] = df.loc[df['sstart'] > df['send'], ['send', 'sstart']].values
    
        # 过滤数据，获取满足条件的 reads 名称，######## 核心的比对算法
        filtered_df = df[(df['sstart'] <= ref_start - check_spanning_length) & (df['send'] >= ref_end + check_spanning_length)]
        
        # 将满足条件的 reads 名称添加到集合中
        selected_reads = set(filtered_df['qseqid'])
    else:
        selected_reads = set()
    
    # 从输入的 FASTA 文件中提取对应的序列并保存到新的 FASTA 文件中
    with open(output_fasta, 'w') as outfile:
        if selected_reads:
            with open(input_fasta, 'r') as infile:
                for record in SeqIO.parse(infile, 'fasta'):
                    if record.id in selected_reads:
                        SeqIO.write(record, outfile, 'fasta')
        else:
            # 如果没有满足条件的 reads，创建一个空的 FASTA 文件
            pass

    # Get the absolute path of the current working directory
    abs_work_dir = os.path.dirname(os.path.realpath(__file__))
    
    # Call the alignment script based on sequencing type and alignment type
    if sequencing_type == "NGS single end" or sequencing_type == "NGS pair end":          ####将blast筛选的read比对到reference上
        subprocess.run(["bash", os.path.join(abs_work_dir, "minimap2.single.sh"), reference, output_fasta,
                             os.path.join(os.getcwd(), output_prefix + "_results"), threads], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True) 
    
    if sequencing_type == "TGS":                   ####将blast筛选的read比对到reference上
        subprocess.run(["bash", os.path.join(abs_work_dir, "minimap2.third.sh"), reference, output_fasta,
                         os.path.join(os.getcwd(), output_prefix + "_results"), threads, seqdepth_type], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True) 

    # 删除blast创建的中间结果
    database_files = glob.glob(f"{output_prefix}_results/{ref_prefix}.n*")
    for database_file in database_files:
        if os.path.exists(database_file):
            os.remove(database_file) 
    # 删除其他中间结果
    if os.path.exists(f"{output_prefix}_results/{ref_prefix}"):
        os.remove(f"{output_prefix}_results/{ref_prefix}") 
    if os.path.exists(f"{output_prefix}_results/selected_reads.fasta"):
        os.remove(f"{output_prefix}_results/selected_reads.fasta")

################################################################################
def check_file_format_efficient(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()

        # 删除空行和去除行首尾的空格
        while not first_line:
            first_line = file.readline().strip()

        # 去除行首尾的空格后再进行格式检查
        first_line = first_line.strip()
        if first_line.startswith('>'):
            return 'FASTA'
        elif first_line.startswith('@'):
            return 'FASTQ'
        else:
            return 'Unknown'
            
################################################################################
def merge_fasta_files(file1, file2, output_file):
    with open(output_file, 'w') as outfile:
        with open(file1, 'r') as infile1:
            outfile.write(infile1.read())
        with open(file2, 'r') as infile2:
            outfile.write(infile2.read())
    return output_file
    
########################################################################################################################################################################
########################################################################################################################################################################
if __name__ == "__main__":
    # Create argparse object and set up command line arguments
    parser = argparse.ArgumentParser(prog="seqdepth.py", description="Draw sequence depth.")

    parser.add_argument("-alignment", choices=["minimap2", "bwa", "blast"], default="minimap2", help="Specifies the alignment software to be used. Default is minimap2.")
    parser.add_argument("-reference", type=str, help="Specifies the reference genome or a specific segment of the genome.")
    parser.add_argument("-single", metavar="FASTQ", help="Specifies single-end sequencing data.")
    parser.add_argument("-pair", metavar=("FASTQ1", "FASTQ2"), nargs=2, help="Specifies paired-end sequencing data.")
    parser.add_argument("-third", metavar="FASTQ", help="Specifies third-generation sequencing data.")
    parser.add_argument("-output", type=str, help="Specify the output prefix for the image file.")
    parser.add_argument("-threads", type=str, help="threads for the mapping softwares.")
    parser.add_argument("-seqdepth_type", type=str, choices=['ont', 'pacbio', None], help="Type for third-generation sequencing data. Choices: 'ont', 'pacbio'.")
    parser.add_argument("-spanning_length", default=1, type=int, help="Spanning length of reads supporting mitogenome recombination.")
    parser.add_argument("-evalue", default="1e-5", type=str, help="Spanning length of reads supporting mitogenome recombination.")

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    alignment = args.alignment
    reference = args.reference
    threads = args.threads
    seqdepth_type = args.seqdepth_type
    check_spanning_length = args.spanning_length

    sequencing_type = ""
    if args.single:
        sequencing_type = "NGS single end"
        single_fastq = args.single
    elif args.pair:
        sequencing_type = "NGS pair end"
        pair_fastq1, pair_fastq2 = args.pair
    elif args.third:
        sequencing_type = "TGS"
        third_fastq = args.third

    output_prefix = "default_output_prefix"
    if args.output:
        output_prefix = args.output

    input_fastq = ""
    if sequencing_type == "NGS single end":
        input_fastq = single_fastq
    elif sequencing_type == "NGS pair end":
        input_fastq = (pair_fastq1, pair_fastq2)
    elif sequencing_type == "TGS":
        input_fastq = third_fastq

    # Create the output folder if it doesn't exist
    output_folder = output_prefix + "_results"
    # Check if the folder exists
    if not os.path.exists(output_folder):
        # Create a new folder
        os.makedirs(output_folder)

    # Copy the reference file to the output folder
    shutil.copy(reference, output_folder)
    
    # 构建完整的文件路径
    output_file = os.path.join(output_folder, 'output.sam')

    # 调用 draw_seq_depth 函数
    if not args.seqdepth_type and alignment != "blast":
        draw_seq_depth_next(alignment, reference, sequencing_type, output_prefix, threads, input_fastq)
    elif args.seqdepth_type and alignment != "blast":
        draw_seq_depth_third(alignment, reference, output_prefix, threads, input_fastq, args.seqdepth_type)
    
    if alignment == "blast":
        if sequencing_type == "NGS single end" or sequencing_type == "TGS":
            if check_file_format_efficient(input_fastq) != "FASTA":         # 检查输入的reads的格式是否为fasta，因为要使用blast软件
                print(f"Error: The format of reads should be 'fasta' when 'blast' was used. But now is {check_file_format_efficient(input_fastq)}.\n")
                exit(1)
            draw_seq_depth_blast(reference, sequencing_type, output_prefix, threads, input_fastq, seqdepth_type, check_spanning_length, args.evalue)
        if sequencing_type == "NGS pair end":
            if check_file_format_efficient(input_fastq[0]) != "FASTA":      # 检查输入的reads的格式是否为fasta，因为要使用blast软件
                print(f"Error: The format of reads should be 'fasta' when 'blast' was used. But now is {check_file_format_efficient(input_fastq[0])}.\n")
                exit(1)
            if check_file_format_efficient(input_fastq[1]) != "FASTA":      # 检查输入的reads的格式是否为fasta，因为要使用blast软件
                print(f"Error: The format of reads should be 'fasta' when 'blast' was used. But now is {check_file_format_efficient(input_fastq[1])}.\n")
                exit(1)
            pair_merged_fasta = f"pair_merged_fasta.fasta"
            pair_merged = merge_fasta_files(input_fastq[0], input_fastq[1], pair_merged_fasta)    # 双端数据合并到一起
            draw_seq_depth_blast(reference, sequencing_type, output_prefix, threads, pair_merged, seqdepth_type, check_spanning_length, args.evalue)
    
            # 删除fastq2fasta的中间结果
            if os.path.exists(pair_merged_fasta):
                os.remove(pair_merged_fasta)
    
