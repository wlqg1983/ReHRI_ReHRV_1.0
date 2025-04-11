#!/usr/bin/env python

import argparse
import csv
import sys
import pandas as pd
import numpy as np
from Bio import SeqIO

#######################################################################################################################################################################
def split_and_modify_data(new_data, id1, id2,output_file):
    # Ensure new_data is a pandas DataFrame
    if not isinstance(new_data, pd.DataFrame):
        new_data = pd.DataFrame(new_data, columns=['fragment_id', 'length'])

    # Find the row numbers for id1 and id2
    id1_index = None
    id2_index = None
    for i, row in new_data.iterrows():
        if row['fragment_id'] == id1:
            id1_index = i
        elif row['fragment_id'] == id2:
            id2_index = i

    # Ensure both IDs were found
    if id1_index is not None and id2_index is not None:
        # Determine the start and end row indices for the sub-list
        start_row = min(id1_index, id2_index)
        end_row = max(id1_index, id2_index)

        # Create the sub-list between id1 and id2
        sub_list = new_data.iloc[start_row:end_row + 1].copy()

        # Vertically rotate the sub-list 180 degrees and multiply the length by -1
        sub_list['length'] = -sub_list['length']

        # Reverse the order of the sub-list
        sub_list = sub_list.iloc[::-1]

        # Split new_data into three sub-lists
        sub_list_1 = new_data.iloc[:start_row]
        sub_list_2 = new_data.iloc[end_row+1:]

        # Combine Remaining Sub-list 1, Modified Sub-list, and Remaining Sub-list 2
        combined_list = pd.concat([sub_list_1, sub_list, sub_list_2])

        # Define the new header
        new_header = ["fragment_id", "start", "end", "direction"]

        # Create a list to store the modified rows, starting with the header
        modified_combined_list = [new_header]

        # Initialize start and end positions
        start = 0
        end = 0

        # Iterate through the rows in combined_list, excluding the header
        for index, row in combined_list.iterrows():
            fragment_id, original_length = row['fragment_id'], row['length']
            original_length = int(original_length)

            # Calculate start and end based on the length
            new_start = end + 1
            new_end = end + abs(original_length)
            end = new_end

            # Determine the direction based on the start and end positions
            direction = 'plus' if original_length > 0 else 'minus'

            # Append the modified row to the new list
            if direction == 'plus':
                modified_combined_list.append([fragment_id, new_start, new_end, direction])
            elif direction == 'minus':
                modified_combined_list.append([fragment_id, new_end, new_start, direction])
            
        # Define the new header for the five-column list for map
        new_header_5CT = ["fragment_id", "start", "end", "type", "paired_id"]

        # Create a list to store the modified rows for the five-column list, starting with the header
        modified_combined_list_5CT = [new_header_5CT]

        # Iterate through the rows in modified_combined_list, excluding the header
        for row in modified_combined_list[1:]:
            fragment_id, start, end, direction = row
            # Determine the type based on the fragment_id
            type = "inter_region" if fragment_id.startswith("ctg") else "-"
            # Append the modified row to the new list for the five-column list
            modified_combined_list_5CT.append([fragment_id, start, end, type, "-"])
        
        # Define the output file name with the specified format
        output_file = f"{output_file}_IR_{id1}_{id2}_5CT.tsv"

        # Save the modified_combined_list_5CT to the output file
        with open(output_file, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')  # Set the delimiter to tab
            writer.writerows(modified_combined_list_5CT)

        return modified_combined_list
            
################################################################################
def update_data_8CT(data_8CT, modified_combined_list, output_file, id1, id2):
    # Create a dictionary from modified_combined_list for easy lookup
    modified_dict = {row[0]: row[1:] for row in modified_combined_list[1:]}

    # Ensure data_8CT is a pandas DataFrame
    if not isinstance(data_8CT, pd.DataFrame):
        data_8CT = pd.DataFrame(data_8CT, columns=['fragment_id', 'start', 'end', 'direction', 'paired_id', 'paired_start', 'paired_end', 'paired_direction'])

    # 在5CT的基础上进行的重组变换，然后根据变换后的5CT，Update data_8CT。
    for index, row in data_8CT.iterrows():
        if row['fragment_id'] in modified_dict:
            data_8CT.loc[index, ['start', 'end', 'direction']] = modified_dict[row['fragment_id']]
        elif row['paired_id'] in modified_dict:
            data_8CT.loc[index, ['paired_start', 'paired_end', 'paired_direction']] = modified_dict[row['paired_id']]

    # Define the output file name
    output_file = f"{output_file}_IR_{id1}_{id2}_map.tsv"

    # Save the updated data_8CT to the output file
    data_8CT.to_csv(output_file, sep='\t', index=False)

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
    
#######################################################################################################################################################################
def update_inputfasta(data_5CT, id1, id2, inputfasta, output_file, genome_type):  

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
    with open(data_5CT, 'r') as file:  
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
        # repeat1_rev, repeat2_rev = repeat2_rev, repeat1_rev
        
    #检查重复单元是否重合
    if repeat1_end >= repeat2_start:
        print(f"Disregard two overlapping repeats marked with {id1, id2}!!")
        return  # Exit the function early

    # 提取 rep1.fasta 与 rep2.fasta 之间的序列并取其反向序列为 revseq1.fasta
    # rep1.fasta 与 rep2.fasta 也一起取了反向重复序列
    seq_between = genome.seq[repeat1_start-1:repeat2_end]
    revseq_sequence = seq_between.reverse_complement()

    # 替换原基因组内对应的序列为 revseq.fasta
    genome_rev = genome.seq[:repeat1_start-1] + revseq_sequence + genome.seq[repeat2_end:] 
    
    # 构造输出文件名  
    output_filename = f"{output_file}_IR_{id1}_{id2}_map.fasta"  
      
    # 将更新后的序列写入到新的fasta文件中  
    with open(output_filename, 'w') as fasta_file:  
        fasta_header = f">IR_{id1}_{repeat1_start}_{repeat1_end}_{id2}_{repeat2_start}_{repeat2_end}_{genome_type}_{len(genome.seq)}"  
        fasta_file.write(f"{fasta_header}\n")
        for i in range(0, len(genome_rev), 100):  
            fasta_file.write(str(genome_rev[i:i+100]) + "\n") 
    
#######################################################################################################################################################################
def main():
    parser = argparse.ArgumentParser(description="Process data from an external file.")
    parser.add_argument('-i', '--input_5CT', help="Path to the input data file with a five-column table. Obtained from the eight-column table.", required=True)
    parser.add_argument('-j', '--input_8CT', help="Path to the input data file with an eight-column table.", required=True)
    parser.add_argument('-f', '--inputfasta', help="Path to the input fasta file.", required=True)
    parser.add_argument('-t', '--genome_type', help="Genome type, linear or circle.", required=True)
    parser.add_argument('-o', '--output_file', help="Path to the output data file.", required=True)
    parser.add_argument('-auto', action='store_true', help="Enable auto mode to automatically process data.")
    args = parser.parse_args()

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # Read the data from the 5-column table input file
    data_5CT = pd.read_csv(args.input_5CT, sep='\t')

    # Create a new data table with 'fragment_id' and 'length'
    # For reverse repeats (end < start), length is negative
    new_data = pd.DataFrame({'fragment_id': data_5CT['fragment_id'],
                             'length': (data_5CT['end'] - data_5CT['start']).apply(lambda x: x + 1 if x >= 0 else x - 1)})

    # Read the data from the 8-column table input file
    data_8CT = pd.read_csv(args.input_8CT, sep='\t')
    data_8CT = data_8CT.drop_duplicates()    # 去除重复项，含空值的行在此处就被去重了
    # 去重复项, 防止用户提供的数据中含有重复项, 已经删除了含空值的行
    data_8CT = drop_duplicates_in_8CT(data_8CT)
    
    # 删除含有空值的行，为查找id做准备。
    data_8CT_clean = data_8CT.dropna()

    # Filter out pairs with different directions and remove pairs with nan values
    filtered_pairs = data_8CT_clean[data_8CT_clean['direction'] != data_8CT_clean['paired_direction']].dropna(subset=['fragment_id', 'paired_id'])[['fragment_id', 'paired_id']].drop_duplicates().values.tolist()

    if not filtered_pairs:
        print("ERROR: No Inverted Repeat found. Please carefully check the repeat type.")
        sys.exit(1)

    if filtered_pairs:
        print()
        print("The following IDs are from paired Inverted Repeats based on your repeat info, which may mediate genome recombination:")
        for i, pair in enumerate(filtered_pairs, start=1):
            print(f"{i}. {pair}")

    if not args.auto and filtered_pairs:
        selected_pairs = []
        error_count = 0
        while filtered_pairs:
            user_input = input("Enter a pair of ABOVE IDs (e.g., RP1a RP1c, case sensitive), 'a/all' to process all, or press 'Enter' Button to exit: ").strip()
            
            if user_input.lower() == 'a' or user_input.lower() == 'all':
                for pair in filtered_pairs:
                    selected_pairs.append(pair)
                    modified_combined_list = split_and_modify_data(new_data, pair[0], pair[1], args.output_file)
                    update_inputfasta(args.input_5CT, pair[0], pair[1], args.inputfasta, args.output_file, args.genome_type)
                    update_data_8CT(data_8CT, modified_combined_list, args.output_file, pair[0], pair[1])
                filtered_pairs.clear()
                print("All Inverted Repeat pairs will be processed.")
                break
            elif user_input == '':
                print("To process entered repeats...") if selected_pairs else None
                print("\nNo IDs entered. Program was terminated by user.") if not selected_pairs else None
                break
            else:
                input_pair = sorted(user_input.split())
                if input_pair in filtered_pairs:
                    print(f"You selected the Inverted Repeats: {input_pair}")
                    selected_pairs.append(input_pair)
                    filtered_pairs.remove(input_pair)
                    modified_combined_list = split_and_modify_data(new_data, input_pair[0], input_pair[1], args.output_file)
                    update_inputfasta(args.input_5CT, input_pair[0], input_pair[1], args.inputfasta, args.output_file, args.genome_type)
                    update_data_8CT(data_8CT, modified_combined_list, args.output_file, input_pair[0], input_pair[1])
                    error_count = 0  # Reset error count after a successful selection
                else:
                    error_count += 1
                    if error_count >= 3:
                        print("ERROR: Too many incorrect attempts. Program will now exit.")
                        sys.exit(1)
                    print("Invalid or duplicate ID pair. Please check again!")
                    
                if len(filtered_pairs) == 0:     # 当用户输入的重复序列对 等于 用户提供的信息中的序列对数量 则自动终止
                    print("All possible pairs have been added.")
                    break 

    if args.auto and filtered_pairs:
        for (id1, id2) in filtered_pairs:
            modified_combined_list = split_and_modify_data(new_data, id1, id2, args.output_file)
            update_inputfasta(args.input_5CT, id1, id2, args.inputfasta, args.output_file, args.genome_type)
            # 升级8CT，用于其他的重组图片的绘制。
            update_data_8CT(data_8CT, modified_combined_list, args.output_file, id1, id2)

################################################################################################################################################################
if __name__ == "__main__":
    main()

