#!/usr/bin/env python

from Bio import SeqIO
import sys

def replace_spaces_in_headers(fasta_file):
    """
    Replace spaces with underscores in the headers of a FASTA file.
    
    Parameters:
    fasta_file (str): Path to the input FASTA file.
    """
    with open(fasta_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # 替换以 '>' 开头的行中的空格为下划线
    modified_lines = [line.replace(' ', '_') if line.startswith('>') else line for line in lines]

    with open(fasta_file, 'w') as f:
        f.writelines(modified_lines)

def modify_headers_in_place(fasta_file):
    """
    Modify the headers of each sequence in the FASTA file by adding chr{i} after '>'.
    
    Parameters:
    fasta_file (str): Path to the input FASTA file.
    """
    records = []
    
    # 读取并修改所有序列的header
    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta"), start=1):
        # Modify the header to include chr{i}
        new_id = f'chr{i}_{record.id}'.replace(' ', '_')
        record.id = new_id
        record.description = ''
        records.append(record)
    
    # 将修改后的序列写回原文件
    with open(fasta_file, 'w', encoding='utf-8') as out_handle:
        SeqIO.write(records, out_handle, "fasta")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python modify_headers.py <fasta_file>")
        sys.exit(1)

    inputfasta = sys.argv[1]  # 从外部参数获取输入文件路径
    
    replace_spaces_in_headers(inputfasta)
    modify_headers_in_place(inputfasta)


