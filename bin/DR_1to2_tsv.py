#!/usr/bin/env python

import argparse
import csv
import sys
import pandas as pd
from Bio import SeqIO
import warnings

#######################################################################################################################################################################
# Filter out FutureWarning for DataFrame concatenation
warnings.filterwarnings("ignore", category=FutureWarning)

##########################################################################################################################################
def find_indices_and_split(new_data, id1, id2, output_file):

    if not isinstance(new_data, pd.DataFrame):
        new_data = pd.DataFrame(new_data, columns=['fragment_id', 'length', 'type', 'paired_id'])
        
    # Find the row numbers for ID1 and ID2 in new_data
    id1_index = None
    id2_index = None
    for i, row in new_data.iterrows():
        if row['fragment_id'] == id1:
            id1_index = i
        elif row['fragment_id'] == id2:
            id2_index = i
            
    # 确保id1在id2之前，因为序列是以线性的方式进行处理的
    if id1_index > id2_index:
        id1_index, id2_index = id2_index, id1_index
        id1, id2 = id2, id1

    # Ensure both IDs were found
    if id1_index is not None and id2_index is not None:
        # Subtable 1: From the first row of new_data to the row containing ID1
        subtable_1 = new_data[:id1_index +1].values.tolist() 

        # Subtable 2: From the row after ID1 to the row containing ID2
        subtable_2 = new_data[id1_index + 1:id2_index + 1]
        header = ["fragment_id", "length", "type", "paired_id"]
        subtable_2 = [header] + subtable_2.values.tolist()

        # Subtable 3: From the row after ID2 to the last row of new_data
        subtable_3 = new_data[id2_index + 1:].values.tolist()

        # Append Subtable 3 to Subtable 1 to create Subtable 4
        subtable_4 = subtable_1 + subtable_3
        subtable_4 = [header] + subtable_4

        # Convert Subtable 2 and Subtable 4 to 5-column tables
        subtable_2_5_column = convert_length_to_start_end(subtable_2)
        subtable_4_5_column = convert_length_to_start_end(subtable_4)

        # Save the converted tables to separate files with the "output" prefix
        output_file_1 = f"{output_file}_DR_{id1}_{id2}_Chr1_1to2.tsv"
        output_file_2 = f"{output_file}_DR_{id1}_{id2}_Chr2_1to2.tsv"
        
        # Convert the 'start' and 'end' columns to integers, skipping the header row
        subtable_2_5_column = [row if row[0] == 'fragment_id' else [row[0], int(row[1]), int(row[2]), row[3], row[4]] for row in subtable_2_5_column]
        subtable_4_5_column = [row if row[0] == 'fragment_id' else [row[0], int(row[1]), int(row[2]), row[3], row[4]] for row in subtable_4_5_column]

        # Define a function to save the table to a file
        def save_table_to_file(table, filename):
            with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')  # Set the delimiter to tab
                writer.writerows(table)

        # Save Subtable 2 as Subconfiguration 1, 线性或者环状，5CT
        save_table_to_file(subtable_4_5_column, output_file_1)
        
        # Save Subtable 4 as Subconfiguration 2, 必须是环状，5CT
        save_table_to_file(subtable_2_5_column, output_file_2)

##########################################################################################################################################
# Function to convert sequence length to start and end positions
def convert_length_to_start_end(table):
    # Define the header for the new 5-column tables
    new_header = ["fragment_id", "start", "end", "type", "paired_id"]
    converted_table = [new_header]
    end = 0

    for row in table[1:]:
        fragment_id, original_length, sequence_type, paired_id = row

        if original_length != 'length':
            original_length = int(original_length)

            # Calculate start and end positions based on the assumption of positive length
            new_start = end + 1
            new_end = end + abs(original_length)
            end = new_end

            # Adjust positions based on the sign of original_length
            if original_length < 0:
                new_start, new_end = new_end, new_start

            converted_table.append([fragment_id, new_start, new_end, sequence_type, paired_id])

    return converted_table

##########################################################################################################################################
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

##########################################################################################################################################
def update_inputfasta(data_5CT, id1, id2, inputfasta, output_file, genome_type): 

    #先检查是不是fasta格式
    file_format = check_file_format_efficient(inputfasta)
    if file_format != 'FASTA':
        print(f"ERROR: The infile {inputfasta} is not in FASTA format.")
        sys.exit(1)
    
    sequences = []  
    with open(inputfasta, 'r', encoding='utf-8') as fasta_file:  
        current_sequence = ""  
        for line in fasta_file:  
            line = line.strip()  
            if line.startswith(">"):  
                if current_sequence:  
                    sequences.append(current_sequence)  
                    current_sequence = ""  
            else:  
                current_sequence += line  
        if current_sequence:  
            sequences.append(current_sequence)  
  
    if len(sequences) > 1:  
        print(f"ERROR: Only one sequence should be in: {inputfasta}.")
        sys.exit(1) 
    elif len(sequences) == 0:  
        print(f"ERROR: At least one sequence should be in: {inputfasta}.")
        sys.exit(1)  
        
        
    genome = SeqIO.read(inputfasta, "fasta")    # 读取fasta格式的文件
    
    if genome_type == "C":
        genome_type = "Circle"
    elif genome_type == "L":
        genome_type = "Linear"
    else:
        print("ERROR: The provided genome type was invalid! It should be C/L.")
        sys.exit(1)

    # 假设data_5CT是一个文件路径，我们需要先从文件中读取数据并构建一个字段映射  
    fields = {}  
    with open(data_5CT, 'r', encoding='utf-8') as file:  
        for line in file:  
            parts = line.strip().split('\t')  
            if len(parts) == 5:  
                fragment_id, start, end, type, paired_id = parts  
                fields[fragment_id] = {  
                    'start': start,  
                    'end': end,  
                    'type': type,  
                    'paired_id': paired_id  
                }  

    # 提取repeat1和repeat2的信息  
    repeat1_info = fields.get(id1)  
    repeat2_info = fields.get(id2)  
      
    if not repeat1_info or not repeat2_info:  
        print(f"One or both of the fragment IDs {id1} and {id2} were not found in the data.")  
        return  
   
    # 转换开始和结束位置为整数  
    repeat1_start_str = repeat1_info['start'] 
    repeat1_start = int(float(repeat1_start_str.replace('bp', '')))
    repeat1_end_str = repeat1_info['end'] 
    repeat1_end = int(float(repeat1_end_str.replace('bp', '')))
    repeat2_start_str = repeat2_info['start'] 
    repeat2_start = int(float(repeat2_start_str.replace('bp', '')))
    repeat2_end_str = repeat2_info['end'] 
    repeat2_end = int(float(repeat2_end_str.replace('bp', '')))
    
    # Check if repeat1_start is greater than repeat1_end
    if repeat1_start > repeat1_end:
        # Swap the values of repeat1_start and repeat1_end
        repeat1_start, repeat1_end = repeat1_end, repeat1_start

    # Check if repeat2_start is greater than repeat2_end
    if repeat2_start > repeat2_end:
        # Swap the values of repeat1_start and repeat1_end
        repeat2_start, repeat2_end = repeat2_end, repeat2_start
        
    # 重新排序以确保repeat1始终在repeat2之前
    if repeat1_start > repeat2_start:
        # 交换repeat1和repeat2的所有信息
        repeat1_start, repeat2_start = repeat2_start, repeat1_start
        repeat1_end, repeat2_end = repeat2_end, repeat1_end
        
    #检查重复单元是否重合
    if repeat1_end >= repeat2_start:
        print(f"Disregard two overlapping repeats marked with {id1, id2}!!")
        return  # Exit the function early

    # 检测 rep1_sequence 上游序列和 rep2_sequence 下游序列的长度
    upstream_length = repeat1_start
    downstream_length = len(genome.seq) - repeat2_end

    # 提取 rep1.fasta 和 rep2.fasta
    rep1_sequence = genome.seq[repeat1_start-1:repeat1_end]
    rep2_sequence = genome.seq[repeat2_start-1:repeat2_end]

    #将基因组拆分为两个亚型
    genome_sub1 = rep1_sequence + genome.seq[repeat2_end:] + genome.seq[:repeat1_start-1]
    genome_sub2 = genome.seq[repeat1_end : repeat2_end] 

    # 构造输出文件名  
    output_file_1 = f"{output_file}_DR_{id1}_{id2}_Chr1_1to2.fasta"
    output_file_2 = f"{output_file}_DR_{id1}_{id2}_Chr2_1to2.fasta"

    # 将更新后的序列写入到新的fasta文件中  
    if genome_type == "Circle":
        with open(output_file_1, 'w', encoding='utf-8') as fasta_file:  
            fasta_header = f">DR_{id1}_{repeat1_start}_{repeat1_end}_{id2}_{repeat2_start}_{repeat2_end}_Chromosome1_{len(genome_sub1)}_{genome_type}_1to2"  
            fasta_file.write(f"{fasta_header}\n")
            for i in range(0, len(genome_sub1), 100):  
                fasta_file.write(str(genome_sub1[i : i+100]) + "\n") 

        with open(output_file_2, 'w', encoding='utf-8') as fasta_file:  
            fasta_header = f">DR_{id1}_{repeat1_start}_{repeat1_end}_{id2}_{repeat2_start}_{repeat2_end}_Chromosome2_{len(genome_sub2)}_{genome_type}_1to2"  
            fasta_file.write(f"{fasta_header}\n")  
            for i in range(0, len(genome_sub2), 100):  
                fasta_file.write(str(genome_sub2[i : i+100]) + "\n") 

    elif genome_type == "Linear":
        with open(output_file_1, 'w', encoding='utf-8') as fasta_file:  
            fasta_header = f">DR_{id1}_{repeat1_start}_{repeat1_end}_{id2}_{repeat2_start}_{repeat2_end}_Chromosome1_{len(genome_sub1)}_{genome_type}_1to2"  
            fasta_file.write(f"{fasta_header}\n")
            for i in range(0, len(genome_sub1), 100):  
                fasta_file.write(str(genome_sub1[i : i+100]) + "\n") 

        with open(output_file_2, 'w', encoding='utf-8') as fasta_file:  
            fasta_header = f">DR_{id1}_{repeat1_start}_{repeat1_end}_{id2}_{repeat2_start}_{repeat2_end}_Chromosome2_{len(genome_sub2)}_Circle_1to2"  
            fasta_file.write(f"{fasta_header}\n")  
            for i in range(0, len(genome_sub2), 100):  
                fasta_file.write(str(genome_sub2[i : i+100]) + "\n") 

################################################################################
def drop_duplicates_in_8CT(merged_file):
    #### 去除8CT中的重复行
    
    # 步骤1: 去除重复行  去除重复的表头
    merged_file = merged_file.drop_duplicates()    # 初级去重
    merged_file = merged_file[~(merged_file.astype(str).apply(lambda row: row.astype(str).str.contains('ctg').any(), axis=1))]    # 去掉间区，以ctg为标记

    # 分割数据为含有空值和不含有空值的两部分  
    nan_rows = merged_file[merged_file.isnull().any(axis=1)]  
    non_nan_rows = merged_file[~merged_file.isnull().any(axis=1)]  

    # 创建一个新列，其中包含排序后的fragment_id和paired_id的元组  
    def create_sorted_id_tuple(row):  
        fragment_id = row['fragment_id']  
        paired_id = row['paired_id']  
        # 对fragment_id和paired_id进行排序  
        sorted_ids = sorted([(fragment_id, paired_id), (paired_id, fragment_id)])  
        # 返回排序后的第一个元组（因为两个元组是排序后的，所以它们中的任何一个都可以作为唯一标识符）  
        return sorted_ids[0]  

    # 假设 non_nan_rows 是一个DataFrame，且其索引是有效的  
    non_nan_rows = non_nan_rows.copy()
    non_nan_rows.loc[:, 'sorted_id_tuple'] = non_nan_rows.apply(create_sorted_id_tuple, axis=1)

    # 查找相同的行并记录index  
    duplicate_dict = {}  
    duplicate_indices = []  
    for index, row in non_nan_rows.iterrows():  
        sorted_id_tuple = row['sorted_id_tuple']  
        if sorted_id_tuple in duplicate_dict:  
            duplicate_indices.append(index)  
        else:  
            duplicate_dict[sorted_id_tuple] = index  

    # 保留每组重复项的第一项（即duplicate_dict中记录的索引对应的行）  
    unique_non_nan_rows = non_nan_rows.drop(duplicate_indices)  
    
    # 最终的DataFrame，包含含有空值的行和去重后的不含空值的行  
    final_merged_file = pd.concat([nan_rows, unique_non_nan_rows], ignore_index=True)
    merged_file = final_merged_file.iloc[:, :8]
    
    return merged_file
    
##########################################################################################################################################
def main():
    parser = argparse.ArgumentParser(description="Process data from an external file.")
    parser.add_argument('-i', '--input_5CT', help="Path to the input data file with a five-column table.", required=True)
    parser.add_argument('-j', '--input_8CT', help="Path to the input data file with an eight-column table.", required=True)
    parser.add_argument('-f', '--inputfasta', help="Path to the input fasta file.", required=True)
    parser.add_argument('-t', '--genome_type', help="Genome type, line or circle.", required=True)
    parser.add_argument('-o', '--output_file', help="Path to the output data file.")
    parser.add_argument('-auto', action='store_true', help="Enable auto mode to automatically process data.")
    args = parser.parse_args()

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # Read the data from the 5-column table input file
    data_5CT = pd.read_csv(args.input_5CT, sep='\t')

    # Create a new data table with four columns: 'fragment_id', 'length', 'type', 'paired_id'
    new_data = pd.DataFrame({'fragment_id': data_5CT['fragment_id'],
                         'length': (data_5CT['end'] - data_5CT['start']).apply(lambda x: x + 1 if x >= 0 else x - 1),
                         'type': data_5CT['type'],
                         'paired_id': data_5CT['paired_id']})

    # Read the data from the 8-column table input file
    data_8CT = pd.read_csv(args.input_8CT, sep='\t')
    data_8CT = drop_duplicates_in_8CT(data_8CT)
    
    # 删除含有空值的行，为查找id做准备。
    data_8CT_clean = data_8CT.dropna()

    filtered_pairs = [
    sorted(pair) for pair in data_8CT_clean[data_8CT_clean['direction'] == data_8CT_clean['paired_direction']]
    .dropna(subset=['fragment_id', 'paired_id'])
    [['fragment_id', 'paired_id']]
    .drop_duplicates()
    .values
    .tolist()]

    if not filtered_pairs:
        print("ATTENTION: No Direct Repeat pairs found. Please check the repeat type carefully.")
        sys.exit(1)

    if not args.auto and filtered_pairs:
        print("The following IDs are from paired Direct Repeats based on your repeat info, which may mediate genome recombination:")
        for i, pair in enumerate(filtered_pairs, start=1):
            print(f"{i}. {pair}")

        selected_pairs = []
        error_count = 0
        while filtered_pairs:
            user_input = input("Enter paired Direct Repeat sequence IDs (e.g., RU1a RU1b, case sensitive), 'a/all' to process all, or press 'Enter' to exit: ").strip()
            
            if user_input.lower() == 'a' or user_input.lower() == 'all':
                for pair in filtered_pairs:
                    print(f"Processing Direct Repeats: {pair}")
                    selected_pairs.append(pair)
                    find_indices_and_split(new_data, pair[0], pair[1], args.output_file)
                    update_inputfasta(args.input_5CT, pair[0], pair[1], args.inputfasta, args.output_file, args.genome_type)
                filtered_pairs.clear()
                print("All Direct Repeat pairs will be processed. Program will now continue...")
                break
            elif user_input == '':
                print("To process entered repeats...") if user_input else None
                print("\nNo IDs entered. Program was terminated by user.") if user_input else None
                break
            else:
                input_pair = sorted(user_input.split())
                if input_pair in filtered_pairs:
                    print(f"You selected the Direct Repeats: {input_pair}")
                    selected_pairs.append(input_pair)
                    filtered_pairs.remove(input_pair)
                    find_indices_and_split(new_data, input_pair[0], input_pair[1], args.output_file)
                    update_inputfasta(args.input_5CT, input_pair[0], input_pair[1], args.inputfasta, args.output_file, args.genome_type)
                    error_count = 0  # Reset error count after a successful selection
                else:
                    error_count += 1
                    if error_count >= 3:
                        print("ERROR: Too many incorrect attempts. Program will now exit.")
                        sys.exit(1)
                    print("Invalid or duplicate ID pair.")
                    
                if len(filtered_pairs) == 0:     # 当用户输入的重复序列对 等于 用户提供的信息中的序列对数量 则自动终止
                    print("All possible pairs have been added.")
                    break 

    if args.auto and filtered_pairs:
        print("The following IDs are from paired Direct Repeats based on your repeat info, which may mediate genome recombination:")
        for i, pair in enumerate(filtered_pairs, start=1):
            print(f"{i}. {pair}")
    
        for (id1, id2) in filtered_pairs:
            find_indices_and_split(new_data, id1, id2, args.output_file)
            update_inputfasta(args.input_5CT, id1, id2, args.inputfasta, args.output_file, args.genome_type)

##########################################################################################################################################
if __name__ == "__main__":
    main()

