#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os, re
import sys

###################################################################################################################################
def reverse_complement_dna(dna_sequence):
    # 检查dna_sequence是否为Seq对象，如果是，则直接调用reverse_complement方法
    if isinstance(dna_sequence, Seq):
        return dna_sequence.reverse_complement()
    # 如果dna_sequence是字符串，则处理简并碱基并考虑大小写
    elif isinstance(dna_sequence, str):
        # 将序列转换为大写以统一处理
        dna_sequence_upper = dna_sequence.upper()
        complement = {
            'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
            'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
            'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N'
        }
        # 验证序列中的每个碱基
        for base in dna_sequence_upper:
            if base not in complement:
                raise ValueError("Fasta file contains illegal base characters.")
                sys.exit(0)
        
        # 生成互补序列
        return ''.join(complement[base] for base in reversed(dna_sequence_upper))
    else:
        raise TypeError("Input must be a Bio.Seq.Seq object or a string.")

########################################################################################################################################################################
def trim_IR_sequence_circle(genome, fields, output_dir, trim_len, id1, id2, rev_compl_IR): 

    repeat1_start = int(float(fields[id1]['start']))  
    repeat1_end = int(float(fields[id1]['end']))  
    repeat1_rev = False
    repeat2_rev = False
    # Check if repeat1_start is greater than repeat1_end
    if repeat1_start > repeat1_end:
        repeat1_rev = True
        # Swap the values of repeat1_start and repeat1_end
        repeat1_start, repeat1_end = repeat1_end, repeat1_start
  
    repeat2_start = int(float(fields[id2]['start']))  
    repeat2_end = int(float(fields[id2]['end']))
    # Check if repeat2_start is greater than repeat2_end
    if repeat2_start > repeat2_end:
        repeat2_rev = True
        # Swap the values of repeat1_start and repeat1_end
        repeat2_start, repeat2_end = repeat2_end, repeat2_start
        
    # 重新排序以确保repeat1始终在repeat2之前
    if repeat1_start > repeat2_start:
        # 交换repeat1和repeat2的所有信息
        repeat1_start, repeat2_start = repeat2_start, repeat1_start
        repeat1_end, repeat2_end = repeat2_end, repeat1_end
        repeat1_rev, repeat2_rev = repeat2_rev, repeat1_rev
        
    #检查重复单元是否重合
    if repeat1_end >= repeat2_start:
        print(f"Disregard two overlapping repeats marked with {id1, id2}!!")
        return  # Exit the function early
    

    # 提取 rep1.fasta 和 rep2.fasta
    rep1_sequence = genome.seq[repeat1_start-1:repeat1_end]
    rep2_sequence = genome.seq[repeat2_start-1:repeat2_end]
    repeat1_length = len(rep1_sequence)
    repeat2_length = len(rep2_sequence)

    trim_len = int(trim_len)
    if len(genome.seq) - len(rep1_sequence) >= trim_len:
        #upstream_length = repeat1_start
        downstream_length = len(genome.seq) - repeat1_end
        upstream_seq = genome.seq[:repeat1_start-1]
        downstream_seq = genome.seq[repeat1_end:]
        
        genome_reset = downstream_seq + genome.seq + upstream_seq  # 是假设的基因组，只为方便截取序列

        # 重置重复序列单元的start和end，最后再截取
        # Update repeat1_start to reflect the change
        repeat1_new_start = downstream_length + repeat1_start
        repeat1_new_end = repeat1_new_start + repeat1_length

        # 从重塑后的基因组中提取AB
        ab_sequence = genome_reset[repeat1_new_start-trim_len-1 : repeat1_new_end + trim_len]
        if not repeat1_rev:
            ab_record = SeqRecord(ab_sequence, id=f"IR_AB_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
        if repeat1_rev:
            ab_sequence = reverse_complement_dna(ab_sequence)
            ab_record = SeqRecord(ab_sequence, id=f"IR_AB_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_IR and not repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"IR_AB_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_IR and repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"IR_AB_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
    else:
        lasso_trim_len = len(genome.seq) - len(rep1_sequence)  #trim_len，计算截取序列的长度
        upstream_seq = genome.seq[:repeat1_start-1]
        downstream_seq = genome.seq[repeat1_end:]
        
        # 从重塑后的基因组中提取AB
        ab_sequence = downstream_seq + genome.seq + upstream_seq  # 假设read已经套圈基因组
        if not repeat1_rev:
            ab_record = SeqRecord(ab_sequence, id=f"IR_AB_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        if repeat1_rev:
            ab_sequence = reverse_complement_dna(ab_sequence)
            ab_record = SeqRecord(ab_sequence, id=f"IR_AB_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_IR and not repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"IR_AB_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_IR and repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"IR_AB_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        sys.stderr.write(f"ATTENTION: The length of the genome containing {id1, id2} repeat is too short ({len(genome.seq)} bp) and the length of the trimmed sequence is too long ({trim_len} bp), so there are not enough sequences in the genome to trim AB sequence. The maximum length of the genome allowed to be trimmed is {len(genome.seq) - len(rep1_sequence)} bp. So, it is assumed here that reads can cover the genome multiple times and reset the length of the trimmed sequence (-tl) to {lasso_trim_len} bp automatically.") 
        
################################################################################
    if len(genome.seq) - len(rep2_sequence) >= trim_len:
        #upstream_length = repeat2_start
        downstream_length = len(genome.seq) - repeat2_end
        upstream_seq = genome.seq[:repeat2_start-1]
        downstream_seq = genome.seq[repeat2_end:]
        #提取CD
        genome_reset = downstream_seq + genome.seq + upstream_seq  # 是假设的基因组，只为方便截取序列
        
        # Update repeat2_start to reflect the change
        repeat2_new_start = downstream_length + repeat2_start
        repeat2_new_end = repeat2_new_start + repeat2_length

        # 从重塑后的基因组中提取CD
        cd_sequence = genome_reset[repeat2_new_start-trim_len-1 : repeat2_new_end+trim_len]
        if not repeat2_rev:
            cd_record = SeqRecord(cd_sequence, id=f"IR_CD_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
        if repeat2_rev:
            cd_sequence = reverse_complement_dna(cd_sequence)
            cd_record = SeqRecord(cd_sequence, id=f"IR_CD_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_IR and not repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"IR_CD_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_IR and repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"IR_CD_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
    else:
        lasso_trim_len = len(genome.seq) - len(rep2_sequence)
        upstream_seq = genome.seq[:repeat2_start]
        downstream_seq = genome.seq[repeat2_end-1:]
        
        # 从重塑后的基因组中提取CD
        cd_sequence = downstream_seq + genome.seq + upstream_seq  # 假设read已经套圈基因组
        if not repeat2_rev:
            cd_record = SeqRecord(cd_sequence, id=f"IR_CD_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        if repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record = SeqRecord(cd_sequence, id=f"IR_CD_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_IR and not repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"IR_CD_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_IR and repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"IR_CD_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        print(f"ATTENTION: The length of the genome containing {id1, id2} repeat is too short ({len(genome.seq)} bp) and the length of the trimmed sequence is too long ({trim_len} bp), so there are not enough sequences in the genome to trim CD sequence. The maximum length of the genome allowed to be trimmed is {len(genome.seq) - len(rep2_sequence)} bp. So, it is assumed here that long reads can cover the genome multiple times and reset the length of the trimmed sequence (-tl) to {lasso_trim_len} bp automatically.")

########################################################################################################################################################################
def trim_IR_sequence_linear(genome, fields, output_dir, trim_len, id1, id2, rev_compl_IR): 

    repeat1_start = int(float(fields[id1]['start']))  
    repeat1_end = int(float(fields[id1]['end']))  
    repeat1_rev = False
    repeat2_rev = False
    # Check if repeat1_start is greater than repeat1_end
    if repeat1_start > repeat1_end:
        repeat1_rev = True
        # Swap the values of repeat1_start and repeat1_end
        repeat1_start, repeat1_end = repeat1_end, repeat1_start
  
    repeat2_start = int(float(fields[id2]['start']))  
    repeat2_end = int(float(fields[id2]['end']))
    # Check if repeat2_start is greater than repeat2_end
    if repeat2_start > repeat2_end:
        repeat2_rev = True
        # Swap the values of repeat1_start and repeat1_end
        repeat2_start, repeat2_end = repeat2_end, repeat2_start
        
    # 重新排序以确保repeat1始终在repeat2之前
    if repeat1_start > repeat2_start:
        # 交换repeat1和repeat2的所有信息
        repeat1_start, repeat2_start = repeat2_start, repeat1_start
        repeat1_end, repeat2_end = repeat2_end, repeat1_end
        repeat1_rev, repeat2_rev = repeat2_rev, repeat1_rev
        
    #检查重复单元是否重合
    if repeat1_end >= repeat2_start:
        print(f"Disregard two overlapping repeats marked with {id1, id2}!!")
        return  # Exit the function early
    

    # 提取 rep1.fasta 和 rep2.fasta
    rep1_sequence = genome.seq[repeat1_start-1:repeat1_end]
    rep2_sequence = genome.seq[repeat2_start-1:repeat2_end]
    repeat1_length = len(rep1_sequence)
    repeat2_length = len(rep2_sequence)

    trim_len = int(trim_len)
    
    #计算repeat1的上下游序列长度
    upstream_seq = genome.seq[:repeat1_start-1]
    downstream_seq = genome.seq[repeat1_end:]
    upstream_length = len(upstream_seq)
    downstream_length = len(downstream_seq)

    if upstream_length >= trim_len and downstream_length >= trim_len:
        # 从重塑后的基因组中提取AB
        ab_sequence = genome.seq[repeat1_start-trim_len-1 : repeat1_end+trim_len]
        if not repeat1_rev:
            ab_record = SeqRecord(ab_sequence, id=f"IR_AB_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
        if repeat1_rev:
            ab_sequence = reverse_complement_dna(ab_sequence)
            ab_record = SeqRecord(ab_sequence, id=f"IR_AB_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_IR and not repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"IR_AB_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_IR and repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"IR_AB_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
    else:
        lasso_trim_len = min(upstream_length, downstream_seq)    #计算截取序列的长度
        
        # 从重塑后的基因组中提取AB
        ab_sequence = genome.seq[repeat1_start-lasso_trim_len-1 : repeat1_end + lasso_trim_len]
        if not repeat1_rev:
            ab_record = SeqRecord(ab_sequence, id=f"IR_AB_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        if repeat1_rev:
            ab_sequence = reverse_complement_dna(ab_sequence)
            ab_record = SeqRecord(ab_sequence, id=f"IR_AB_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_IR and not repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"IR_AB_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_IR and repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"IR_AB_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"IR_AB_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        print(f"ATTENTION: The length of the subgenome1 containing {id1, id2} repeat is too short ({len(genome.seq)} bp) and the length of the trimmed sequence is too long ({trim_len} bp), so there are not enough sequences in the genome to trim AB sequence. Then the length of trimmed sequence (-tl) is forced to set {lasso_trim_len} bp.") 
        
################################################################################
    
    #计算repeat1的上下游序列长度
    upstream_seq = genome.seq[:repeat2_start-1]
    downstream_seq = genome.seq[repeat2_end:]
    upstream_length = len(upstream_seq)
    downstream_length = len(downstream_seq)
    
    if upstream_length >= trim_len and downstream_length >= trim_len:
        # 从重塑后的基因组中提取CD
        cd_sequence = genome.seq[repeat2_start-trim_len-1 : repeat2_end+trim_len]
        if not repeat2_rev:
            cd_record = SeqRecord(cd_sequence, id=f"IR_CD_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
        if repeat2_rev:
            cd_sequence = reverse_complement_dna(cd_sequence)
            cd_record = SeqRecord(cd_sequence, id=f"IR_CD_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_IR and not repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"IR_CD_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_IR and repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"IR_CD_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
    else:
        lasso_trim_len = min(upstream_length, downstream_length)

        # 从重塑后的基因组中提取CD
        cd_sequence = genome.seq[repeat2_start-lasso_trim_len-1 : repeat2_end + lasso_trim_len]
        if not repeat2_rev:
            cd_record = SeqRecord(cd_sequence, id=f"IR_CD_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        if repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record = SeqRecord(cd_sequence, id=f"IR_CD_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_IR and not repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"IR_CD_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_IR and repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"IR_CD_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"IR_CD_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        print(f"ATTENTION: The length of the subgenome1 containing {id1, id2} repeat is too short ({len(genome.seq)} bp) and the length of the trimmed sequence is too long ({trim_len} bp), so there are not enough sequences in the genome to trim CD sequence. Then the length of trimmed sequence (-tl) is forced to set {lasso_trim_len} bp.") 
    
#######################################################################################################################################################################

def trim_DR_sequence_circle(genome, fields, output_dir, trim_len, id1, id2, rev_compl_DR):

    repeat1_start = int(float(fields[id1]['start'])) 
    repeat1_end = int(float(fields[id1]['end'])) 
    repeat1_rev = False

    # Check if repeat1_start is greater than repeat1_end
    if repeat1_start > repeat1_end:
        repeat1_rev = True
        # Swap the values of repeat1_start and repeat1_end
        repeat1_start, repeat1_end = repeat1_end, repeat1_start
  
    repeat2_start = int(float(fields[id2]['start'])) 
    repeat2_end = int(float(fields[id2]['end']))
    repeat2_rev = False
    
    # Check if repeat2_start is greater than repeat2_end
    if repeat2_start > repeat2_end:
        repeat2_rev = True
        # Swap the values of repeat1_start and repeat1_end
        repeat2_start, repeat2_end = repeat2_end, repeat2_start
        
    # 重新排序以确保repeat1始终在repeat2之前
    if repeat1_start > repeat2_start:
        # 交换repeat1和repeat2的所有信息
        repeat1_start, repeat2_start = repeat2_start, repeat1_start
        repeat1_end, repeat2_end = repeat2_end, repeat1_end
        repeat1_rev, repeat2_rev = repeat2_rev, repeat1_rev
        
    #检查重复单元是否重合
    if repeat1_end >= repeat2_start:
        print(f"Disregard two overlapping repeats marked with {id1, id2}!!")
        return  # Exit the function early
    

    # 提取 rep1.fasta 和 rep2.fasta
    rep1_sequence = genome.seq[repeat1_start-1:repeat1_end]
    rep2_sequence = genome.seq[repeat2_start-1:repeat2_end]
    repeat1_length = len(rep1_sequence)
    repeat2_length = len(rep2_sequence)
    
    trim_len = int(float(trim_len))
        
    if len(genome.seq) - len(rep1_sequence) >= trim_len:
        # 提取原理：以一个重复序列单元为中心，左右两侧是重复序列单元以外的所有序列，即左右两侧的序列是相同的
        # 检测 rep1_sequence 上游序列和 rep2_sequence 下游序列的长度
        #upstream_length = repeat1_start
        downstream_length = len(genome.seq) - repeat1_end
        upstream_seq = genome.seq[:repeat1_start-1]
        downstream_seq = genome.seq[repeat1_end:]
        
        genome_reset = downstream_seq + genome.seq + upstream_seq  # 是假设的基因组，只为方便截取序列
        
        # 重置重复序列单元的start和end，最后再截取
        # Update repeat1_start to reflect the change
        repeat1_new_start = downstream_length + repeat1_start
        repeat1_new_end = repeat1_new_start + repeat1_length

        # 从重塑后的基因组中提取AB
        ab_sequence = genome_reset[repeat1_new_start-trim_len-1 : repeat1_new_end + trim_len]
        if not repeat1_rev:
            ab_record = SeqRecord(ab_sequence, id=f"DR_AB_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
        if repeat1_rev:
            ab_sequence = reverse_complement_dna(ab_sequence)
            ab_record = SeqRecord(ab_sequence, id=f"DR_AB_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_DR and not repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"DR_AB_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_DR and repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"DR_AB_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
    else:
        lasso_trim_len = len(genome.seq) - len(rep1_sequence)  #trim_len，计算截取序列的长度
        upstream_seq = genome.seq[:repeat1_start-1]
        downstream_seq = genome.seq[repeat1_end:]
        
        # 从重塑后的基因组中提取AB
        ab_sequence = downstream_seq + genome.seq + upstream_seq  # 假设read已经套圈基因组
        if not repeat1_rev:
            ab_record = SeqRecord(ab_sequence, id=f"DR_AB_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        if repeat1_rev:
            ab_sequence = reverse_complement_dna(ab_sequence)
            ab_record = SeqRecord(ab_sequence, id=f"DR_AB_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_DR and not repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"DR_AB_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_DR and repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"DR_AB_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        print(f"ATTENTION: The length of the genome containing {id1, id2} repeat is too short ({len(genome.seq)} bp) and the length of the trimmed sequence is too long ({trim_len} bp), so there are not enough sequences in the genome to trim AB sequence. The maximum length of the genome allowed to be trimmed is {len(genome.seq) - len(rep1_sequence)} bp. So, it is assumed here that reads can cover the genome multiple times and reset the length of the trimmed sequence (-tl) to {lasso_trim_len} bp automatically.") 

################################################################################
    if len(genome.seq) - len(rep2_sequence) >= trim_len:
        #upstream_length = repeat2_start
        downstream_length = len(genome.seq) - repeat2_end
        upstream_seq = genome.seq[:repeat2_start-1]
        downstream_seq = genome.seq[repeat2_end:]
        #提取CD
        genome_reset = downstream_seq + genome.seq + upstream_seq  # 是假设的基因组，只为方便截取序列
        
        # Update repeat2_start to reflect the change
        repeat2_new_start = downstream_length + repeat2_start
        repeat2_new_end = repeat2_new_start + repeat2_length

        # 从重塑后的基因组中提取CD
        cd_sequence = genome_reset[repeat2_new_start-trim_len-1 : repeat2_new_end+trim_len]
        if not repeat2_rev:
            cd_record = SeqRecord(cd_sequence, id=f"DR_CD_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
        if repeat2_rev:
            cd_sequence = reverse_complement_dna(cd_sequence)
            cd_record = SeqRecord(cd_sequence, id=f"DR_CD_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_DR and not repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"DR_CD_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_DR and repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"DR_CD_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
    else:
        lasso_trim_len = len(genome.seq) - len(rep2_sequence)
        upstream_seq = genome.seq[:repeat2_start-1]
        downstream_seq = genome.seq[repeat2_end:]
        
        # 从重塑后的基因组中提取CD
        cd_sequence = downstream_seq + genome.seq + upstream_seq  # 假设read已经套圈基因组
        if not repeat2_rev:
            cd_record = SeqRecord(cd_sequence, id=f"DR_CD_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        if repeat2_rev:
            cd_sequence = reverse_complement_dna(cd_sequence)
            cd_record = SeqRecord(cd_sequence, id=f"DR_CD_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_DR and not repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"DR_CD_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_DR and repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"DR_CD_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        print(f"ATTENTION: The length of the genome containing {id1, id2} repeat is too short ({len(genome.seq)} bp) and the length of the trimmed sequence is too long ({trim_len} bp), so there are not enough sequences in the genome to trim CD sequence. The maximum length of the genome allowed to be trimmed is {len(genome.seq) - len(rep1_sequence)} bp. So, it is assumed here that reads can cover the genome multiple times and reset the length of the trimmed sequence (-tl) to {lasso_trim_len} bp automatically.") 

#######################################################################################################################################################################
def trim_DR_sequence_linear(genome, fields, output_dir, trim_len, id1, id2, rev_compl_DR):

    repeat1_start = int(float(fields[id1]['start'])) 
    repeat1_end = int(float(fields[id1]['end'])) 
    repeat1_rev = False

    # Check if repeat1_start is greater than repeat1_end
    if repeat1_start > repeat1_end:
        repeat1_rev = True
        # Swap the values of repeat1_start and repeat1_end
        repeat1_start, repeat1_end = repeat1_end, repeat1_start
  
    repeat2_start = int(float(fields[id2]['start'])) 
    repeat2_end = int(float(fields[id2]['end']))
    repeat2_rev = False
    
    # Check if repeat2_start is greater than repeat2_end
    if repeat2_start > repeat2_end:
        repeat2_rev = True
        # Swap the values of repeat1_start and repeat1_end
        repeat2_start, repeat2_end = repeat2_end, repeat2_start
        
    # 重新排序以确保repeat1始终在repeat2之前
    if repeat1_start > repeat2_start:
        # 交换repeat1和repeat2的所有信息
        repeat1_start, repeat2_start = repeat2_start, repeat1_start
        repeat1_end, repeat2_end = repeat2_end, repeat1_end
        repeat1_rev, repeat2_rev = repeat2_rev, repeat1_rev
        
    #检查重复单元是否重合
    if repeat1_end >= repeat2_start:
        print(f"Disregard two overlapping repeats marked with {id1, id2}!!")
        return  # Exit the function early
    
    # 提取 rep1.fasta 和 rep2.fasta
    rep1_sequence = genome.seq[repeat1_start-1:repeat1_end]
    rep2_sequence = genome.seq[repeat2_start-1:repeat2_end]
    repeat1_length = len(rep1_sequence)
    repeat2_length = len(rep2_sequence)
    
    trim_len = int(float(trim_len))

    # 计算重复单元1的上下游序列的长度
    upstream_seq = genome.seq[:repeat1_start-1]
    downstream_seq = genome.seq[repeat1_end:]
    upstream_length = len(upstream_seq)
    downstream_length = len(downstream_seq)

    if upstream_length >= trim_len and downstream_length >= trim_len:
        # 从重塑后的基因组中提取AB
        ab_sequence = genome.seq[repeat1_start-trim_len-1 : repeat1_end+trim_len]
        if not repeat1_rev:
            ab_record = SeqRecord(ab_sequence, id=f"DR_AB_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
        if repeat1_rev:
            ab_sequence = reverse_complement_dna(ab_sequence)
            ab_record = SeqRecord(ab_sequence, id=f"DR_AB_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_DR and not repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"DR_AB_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_DR and repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"DR_AB_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
    else:
        lasso_trim_len = min(upstream_length, downstream_length)  #计算截取序列的长度
        
        # 从重塑后的基因组中提取AB
        ab_sequence = genome.seq[repeat1_start-lasso_trim_len-1 : repeat1_end+lasso_trim_len]
        if not repeat1_rev:
            ab_record = SeqRecord(ab_sequence, id=f"DR_AB_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        if repeat1_rev:
            ab_sequence = reverse_complement_dna(ab_sequence)
            ab_record = SeqRecord(ab_sequence, id=f"DR_AB_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_DR and not repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"DR_AB_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_DR and repeat1_rev:
            ab_sequence_c = reverse_complement_dna(ab_sequence)
            ab_record_c = SeqRecord(ab_sequence_c, id=f"DR_AB_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(ab_record_c, os.path.join(output_dir, f"DR_AB_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        print(f"ATTENTION: The length of the subgenome1 containing {id1, id2} repeat is too short ({len(genome.seq)} bp) and the length of the trimmed sequence is too long ({trim_len} bp), so there are not enough sequences in the genome to trim AB sequence. Then the length of trimmed sequence (-tl) is forced to set {lasso_trim_len} bp.") 

################################################################################
    # 计算重复单元2的上下游序列的长度
    upstream_seq = genome.seq[:repeat2_start-1]
    downstream_seq = genome.seq[repeat2_end:]
    upstream_length = len(upstream_seq)
    downstream_length = len(downstream_seq)
    
    if upstream_length >= trim_len and downstream_length >= trim_len:
        #提取CD
        cd_sequence = genome.seq[repeat2_start-trim_len-1 : repeat2_end+trim_len]
        if not repeat2_rev:
            cd_record = SeqRecord(cd_sequence, id=f"DR_CD_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
        if repeat2_rev:
            cd_sequence = reverse_complement_dna(cd_sequence)
            cd_record = SeqRecord(cd_sequence, id=f"DR_CD_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_DR and not repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"DR_CD_{id1}_{id2}_minus_{trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_minus_{trim_len}.fasta"), "fasta")
        if rev_compl_DR and repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"DR_CD_{id1}_{id2}_plus_{trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_plus_{trim_len}.fasta"), "fasta")
    else:
        lasso_trim_len = min(upstream_length, downstream_length)
        
        # 从重塑后的基因组中提取CD
        cd_sequence = genome.seq[repeat2_start-lasso_trim_len-1 : repeat2_end+lasso_trim_len]
        if not repeat2_rev:
            cd_record = SeqRecord(cd_sequence, id=f"DR_CD_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        if repeat2_rev:
            cd_sequence = reverse_complement_dna(cd_sequence)
            cd_record = SeqRecord(cd_sequence, id=f"DR_CD_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_DR and not repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"DR_CD_{id1}_{id2}_minus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_minus_{lasso_trim_len}.fasta"), "fasta")
        if rev_compl_DR and repeat2_rev:
            cd_sequence_c = reverse_complement_dna(cd_sequence)
            cd_record_c = SeqRecord(cd_sequence_c, id=f"DR_CD_{id1}_{id2}_plus_{lasso_trim_len}", description="")
            SeqIO.write(cd_record_c, os.path.join(output_dir, f"DR_CD_{id1}_{id2}_plus_{lasso_trim_len}.fasta"), "fasta")
        print(f"ATTENTION: The length of the subgenome1 containing {id1, id2} repeat is too short ({len(genome.seq)} bp) and the length of the trimmed sequence is too long ({trim_len} bp), so there are not enough sequences in the genome to trim CD sequence. Then the length of trimmed sequence (-tl) is forced to set {lasso_trim_len} bp.") 

#######################################################################################################################################################################
def find_all_inverted_repeat_pairs(new_data):
    inverted_repeat_pairs = []  # Create a list to store Inverted_repeat pairs
    inverted_repeat_count = 0
    inverted_repeat_found = False

    for row in new_data[1:]:
        fragment_id = row[0]
        sequence_type = row[3].lower()
        paired_id = row[4]

        if sequence_type == 'inverted_repeat':
            if paired_id != '-':
                inverted_repeat_pairs.append((fragment_id, paired_id))
                inverted_repeat_found = True
                inverted_repeat_count += 1

    if inverted_repeat_count % 2 != 0:
        print(f"Inverted Repeat {inverted_repeat_pairs} are odd and not paired. ")
        sys.exit(1)

    # 7. Sort the list of tuples based on the first element (fragment_id)
    inverted_repeat_pairs = sorted(inverted_repeat_pairs, key=lambda x: x[0])

    # 8. Sort the list after removing duplicates while considering case and sorting elements within each tuple
    inverted_repeat_pairs = [tuple(sorted(pair)) for pair in inverted_repeat_pairs]
    inverted_repeat_pairs = sorted(inverted_repeat_pairs)
    original_pairs_count = len(inverted_repeat_pairs)

    # 9. Remove duplicates by converting the list to a set and back to a list
    inverted_repeat_pairs = list(set(inverted_repeat_pairs))
    sorted_pairs_count = len(inverted_repeat_pairs)
    
    if sorted_pairs_count != original_pairs_count/2:
        print(f"Not all repetitive sequences have two repeat units. Check {inverted_repeat_pairs}.")
        sys.exit(1)  # Terminate the programme and return status code 1. 

    return inverted_repeat_pairs

##########################################################################################################################
def find_all_direct_repeat_pairs(new_data):
    error_count = 0
    direct_repeat_found = False  # Flag to track if direct repeat pairs are found
    id_pairs = []

    while error_count <= 3:
        try:
            # 6. Find Direct Repeat sequences
            direct_repeat_pairs = []  # Create a list to store Direct Repeat pairs

            direct_repeat_count = 0
            for row in new_data[1:]:
                fragment_id = row[0]
                sequence_type = row[3].lower()
                paired_id = row[4]

                # Check if the sequence is a Direct Repeat
                if sequence_type == 'direct_repeat':
                    # Check if paired_id is not '-'
                    if paired_id != '-':
                        # Append a tuple with (fragment_id, paired_id)
                        direct_repeat_pairs.append((fragment_id, paired_id))
                        direct_repeat_found = True  # Set the flag to True
                        direct_repeat_count += 1  # Increment the count

            if direct_repeat_count % 2 != 0:
                print(f"Direct Repeat {direct_repeat_pairs} are odd and not paired.")
                sys.exit(1)  # Terminate the program with an error code

            # 7. Sort the list of tuples based on the first element (fragment_id)
            direct_repeat_pairs = sorted(direct_repeat_pairs, key=lambda x: x[0])

            # 8. Sort the list after removing duplicates while considering case and sorting elements within each tuple
            direct_repeat_pairs = [tuple(sorted(pair)) for pair in direct_repeat_pairs]
            direct_repeat_pairs = sorted(direct_repeat_pairs)
            original_pairs_count = len(direct_repeat_pairs)

            # 9. Remove duplicates by converting the list to a set and back to a list
            direct_repeat_pairs = list(set(direct_repeat_pairs))
            sorted_pairs_count = len(direct_repeat_pairs)
            # Check if the number of pairs after sorting is twice the number of original pairs
            if sorted_pairs_count != original_pairs_count / 2:
                print("IDs may be written incorrectly and the number of paired Direct Repeats does not match.")
                sys.exit(1)  # Terminate the program and return a status code of 1.

            # 10. Iterate through the found ID pairs
            for pair in direct_repeat_pairs:
                id1, id2 = pair
                #find_indices_and_split(new_data, id1, id2, output_file)
                id_pairs.append((id1, id2))

            # Exit the loop
            break

        except ValueError:
            print("Invalid input. Please ensure both IDs are present in the original data.")
            error_count += 1  # Increment error count
    return direct_repeat_pairs

##########################################################################################################################
def find_rows_with_id(data, id1, id2):  
    result = []  
    for row in data:  
        if row[0] == id1 and row[4] == id2:  
            result.append(row)  
        elif row[0] == id2 and row[4] == id1:  
            result.append(row)  
    return result  
    
###################################################################################################################################
def filter_unique_repeats(data):
    # 创建一个字典来存储每个 fragment_id 和 paired_id 组合的出现次数
    repeat_counts = {}

    # 遍历数据，统计每个 fragment_id 和 paired_id 组合的出现次数
    for row in data:
        fragment_id = row[0]
        paired_id = row[4]
        repeat_type = row[3].lower()

        if repeat_type in ['inverted_repeat', 'direct_repeat']:
            pair = tuple(sorted((fragment_id, paired_id)))  # 对 fragment_id 和 paired_id 进行排序

            if pair in repeat_counts:
                repeat_counts[pair] += 1
            else:
                repeat_counts[pair] = 1

    # 创建一个新的数据集来存储仅出现一次的行
    unique_repeats = []

    # 创建一个新的数据集来存储剩余的行
    new_data = []

    for row in data:
        fragment_id = row[0]
        paired_id = row[4]
        repeat_type = row[3].lower()

        if repeat_type in ['inverted_repeat', 'direct_repeat']:
            pair = tuple(sorted((fragment_id, paired_id)))  # 对 fragment_id 和 paired_id 进行排序

            # 如果这个组合仅出现一次，则加入 unique_repeats
            if repeat_counts[pair] == 1:
                unique_repeats.append(row)
            else:
                new_data.append(row)
        else:
            # 如果不是重复序列类型，直接加入 new_data
            new_data.append(row)

    return new_data, unique_repeats

########################################################################################################################################################################
def parse_extramaincon_compl_chain(value):
    if isinstance(value, str):
        normalized_value = value.strip().lower()
        if normalized_value == "yes" or normalized_value == "y":
            return True
        elif normalized_value == "no" or normalized_value == "n":
            return False
        else:
            raise ValueError("Invalid input: must be 'Yes', 'Y', 'No' or 'N'.")
    else:
        raise TypeError("Input must be a string")

########################################################################################################################################################################
def read_last_line(filename):
    with open(filename, 'rb') as file:
        file.seek(-2, os.SEEK_END)  # 跳到文件的倒数第二个字节
        while file.read(1) != b'\n':  # 向后查找第一个换行符
            file.seek(-2, os.SEEK_CUR)
        return file.readline().decode()
        
########################################################################################################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process repeat sequences and extract flanking sequences.")

    parser.add_argument("-r", "--reference", required=True, dest="genome_file", help="Path to the genome FASTA file.")
    parser.add_argument('-rf', '--repeattable', required=True, help="Path to the input data file. A five-colume table.")
    parser.add_argument("-o", "--output_dir_prefix", required=True, help="Prefix for the output directory and files.")
    parser.add_argument("-tl", "--trim_length", default=1000, type=int, help="Length to trim the sequences.")
    parser.add_argument("-rc", "--reverse_complement", required=True, type=str, help="If to consider complementary chains (Yes/No).")
    parser.add_argument("-gt", "--genome_type", required=True, help="The type of the input genome. Circle (C) or linear (L).")

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # 读取基因组FASTA文件
    genome = SeqIO.read(args.genome_file, "fasta")
    # 读取基因组的长度
    genome_length = len(genome.seq)

    # 检测基因组的长度与提供的列表是否一致。
    try:
        last_line = read_last_line(args.repeattable)  # 读取最后一行
        columns = last_line.strip().split('\t')       # 分割最后一行以获取列
        if len(columns) >= 5:  # 假设 start 和 end 是第 2 和第 3 列（基于从0开始的索引）
            # 解析 start 和 end 值
            start = int(float(columns[1]))  # 转换 start 值
            end = int(float(columns[2]))    # 转换 end 值
            max_value = max(start, end)     # 找到最大值
            if not max_value == genome_length:
                print(f"The genome length is {genome_length}, which is different the length ({max_value}) in the five-colume table.")
                sys.exit(1)
    except OSError as e:
        print(f"Error reading file: {e}")

    # 检查输出目录是否存在
    if not os.path.isdir(args.output_dir_prefix):
        # 创建新的输出目录
        os.makedirs(args.output_dir_prefix)

    # 1. Process the data from the input file
    with open(args.repeattable, 'r') as file:
        lines = file.readlines()
    
    # 2. Process the original data, lowercase only the headers and type column
    data = [lines[0].strip().lower().split('\t')]  # Process the header
    data += [line.strip().split('\t') for line in lines[1:]]  # Process the data rows
    
    # 现在 new_data 是更新后的数据集，unique_repeats 包含重复单元仅出现一次的重复序列
    data, unique_repeats = filter_unique_repeats(data)
    
    if unique_repeats:
        print(f"Group {unique_repeats} does not have enough items for pairing. Each group must have at least 2 items.")
        sys.exit(1)
    
    inverted_repeat_pairs_ID = find_all_inverted_repeat_pairs(data)
    sequences = {}  
    for pair in inverted_repeat_pairs_ID:  
        id1, id2 = pair  
        rows = find_rows_with_id(data, id1, id2)  
        if rows:  
            sequences[id1] = {'fragment_id': rows[0][0], 'start': rows[0][1], 'end': rows[0][2], 'type': rows[0][3], 'paired_id': rows[0][4]}  
            sequences[id2] = {'fragment_id': rows[1][0], 'start': rows[1][1], 'end': rows[1][2], 'type': rows[1][3], 'paired_id': rows[1][4]}  
            
        rev_compl_IR = parse_extramaincon_compl_chain(args.reverse_complement)
        
        if args.genome_type == "C":
            trim_IR_sequence_circle(genome, sequences, args.output_dir_prefix, args.trim_length, id1, id2, rev_compl_IR)
        if args.genome_type == "L":
            trim_IR_sequence_linear(genome, sequences, args.output_dir_prefix, args.trim_length, id1, id2, rev_compl_IR)

    direct_repeat_pairs_ID = find_all_direct_repeat_pairs(data)
    sequences = {}  
    for pair in direct_repeat_pairs_ID:  
        id1, id2 = pair  
        rows = find_rows_with_id(data, id1, id2)  
        if rows:  
            sequences[id1] = {'fragment_id': rows[0][0], 'start': rows[0][1], 'end': rows[0][2], 'type': rows[0][3], 'paired_id': rows[0][4]}  
            sequences[id2] = {'fragment_id': rows[1][0], 'start': rows[1][1], 'end': rows[1][2], 'type': rows[1][3], 'paired_id': rows[1][4]}  
            
        rev_compl_DR = parse_extramaincon_compl_chain(args.reverse_complement)

        if args.genome_type == "C":
            trim_DR_sequence_circle(genome, sequences, args.output_dir_prefix, args.trim_length, id1, id2, rev_compl_DR)
        if args.genome_type == "L":
            trim_DR_sequence_linear(genome, sequences, args.output_dir_prefix, args.trim_length, id1, id2, rev_compl_DR)

