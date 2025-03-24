#!/usr/bin/env python

import argparse, csv
import sys, os, re, time
import pandas as pd
import numpy as np
from Bio import SeqIO 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord

##########################################################################################################################################
def complementary_strand_5CT_and_fasta(data):
    # 将5CT由plus转换为minus
    df = pd.DataFrame(data)  
  
    # 步骤1: 将length列中的每个值乘以-1  
    df['length'] = df['length'] * -1  
  
    # 步骤2: 将direction列中的'minus'替换为'plus'，'plus'替换为'minus'  
    df['direction'] = df['direction'].replace({'minus': 'plus_temp', 'plus': 'minus'})  # 先将minus替换为临时值  
    df['direction'] = df['direction'].replace('plus_temp', 'plus')  # 再将临时值替换为plus，以避免直接替换导致的冲突  
  
    # 步骤3: 将整个数据框上下翻转180度  
    df = df.iloc[::-1]  # 使用降序索引来翻转数据框 
    
    return df

##########################################################################################################################################
def save_complement_sequences(fasta_file, output_file_path):  
    # 将5CT由plus转换为minus时，对应的fasta序列要取互补链
    # Check if the input file is in FASTA format
    file_format = check_file_format_efficient(fasta_file)
    if file_format != 'FASTA':
        print(f"ERROR: The infile {fasta_file} is not in FASTA format.")
        sys.exit(1)
        
    # Read the sequences from the FASTA file
    sequences = []  
    with open(fasta_file, 'r') as fasta_file:  
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
        print(f"ERROR: Only one sequence should be in: {fasta_file}.")
        sys.exit(1) 
    elif len(sequences) == 0:  
        print(f"ERROR: No sequence found in: {fasta_file}.")
        sys.exit(1)  

    # Write the reverse complement of the sequence to the output file
    with open(output_file_path, "w") as output_file:
        for seq in sequences:
            record = SeqRecord(Seq(seq), id="sequence", description="")
            reverse_complement = record.seq.reverse_complement()
            output_file.write(f">sequence\n{reverse_complement}\n")
            
    return output_file_path

##########################################################################################################################################
def split_data_by_direct_repeat_line(df, id1):
    # unique_ids_with_direction_1：从输入的8CT转换来的5CT，ctg和重复单元用长度表示了，此步是根据id1拆分这个5CT，可以直接绘图
    # Find the index of the row with the direct repeat ID
    direct_repeat_index = df.index[df['fragment_id'] == id1].tolist()

    # If the direct repeat ID is found, split the DataFrame
    if len(direct_repeat_index) == 1:
        # The split index is the index of the direct repeat ID
        split_index = direct_repeat_index[0]
        
        # Split the DataFrame into two parts
        first_subset = df.iloc[:split_index+1]
        second_subset = df.iloc[split_index+1:]

        return first_subset, second_subset
    else:
        print(f"ERROR: Direct repeat ID {id1} is not found in the {df}.")
        sys.exit(1)

##########################################################################################################################################
def split_fasta_by_direct_repeat_line(chrom1_filter, id1, long_genome_seq):
    # 将fasta序列根据id的起点和终点切分
    # chrom1_filter的标头：['fragment_id', 'start', 'end', 'direction', 'paired_id']
    # chrom1_filter：从输入的8CT转换来的5CT，ctg和重复单元用start和end表示的，是unique_ids_with_direction_1之前的数据
    # 将chrom1_filter转换成DataFrame，如果它还不是  
    df_filter = pd.DataFrame(chrom1_filter)  
      
    # 查找id1对应的记录  
    target_fragment = df_filter[df_filter['fragment_id'] == id1]  
      
    if target_fragment.empty:  
        raise ValueError(f"No fragment found with id: {id1}")  
      
    # 假定只有一个匹配项，获取start和end  
    start = target_fragment['start'].iloc[0]  
    end = target_fragment['end'].iloc[0]  

    #先检查是不是fasta格式
    file_format = check_file_format_efficient(long_genome_seq)
    if file_format != 'FASTA':
        print(f"The infile {long_genome_seq} is not in FASTA format.")
        sys.exit(1)
    
    sequences = []  
    with open(long_genome_seq, 'r') as fasta_file:  
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
        print(f"Only one sequence should be in: {long_genome_seq}.")
        sys.exit(1) 
    elif len(sequences) == 0:  
        print(f"Only one sequence should be in: {long_genome_seq}.")
        sys.exit(1)  
    
    genome_seq = sequences[0]
    
    # 根据start和end拆分原始序列  
    seq_before_end = genome_seq[:end+1]  # 包含end位置的碱基  
    seq_after_end = genome_seq[end+1:]  
    
    # 计算seq_before_end的长度  
    seq_before_end_length = len(seq_before_end)
    seq_after_end_length = len(seq_before_end)

    return seq_before_end, seq_after_end, seq_before_end_length, seq_after_end_length, start, end

##########################################################################################################################################
def convert_length_to_start_end(table):
    # Define the header for the new 4-column table
    new_header = ["fragment_id", "start", "end", "direction"]
    converted_table = [new_header]
    end = 0

    for index, row in table.iterrows():
        fragment_id, original_length, direction = row['fragment_id'], row['length'], row['direction']
        
        if original_length != 'length':  # Skip the header row
            original_length = int(original_length)
            
            # 计算起点和终点
            new_start = end + 1
            new_end = end + abs(original_length)

            # Adjust the start and end positions based on the direction
            if direction == 'minus':
                new_start, new_end = new_end, new_start

            converted_table.append([fragment_id, new_start, new_end, direction])
            end = max(new_end, new_start)  # Update the end position for the next iteration

    return converted_table

##########################################################################################################################################
def add_ctg(processed_data, genome_length, chr):
    
    final_data = []
    prev_end = 0
    ctg_counter = 1
    
    first_start = 0
    last_end = 0

    for data in processed_data: #[1:]:  # Skip the header row
        fragment_id, start, end, direction = data[:4]

        adjusted_start = min(int(start), int(end))
        adjusted_end = max(int(start), int(end))
        
        # 查找绘图数据中的起点和终点
        if ctg_counter == 1:    # 先赋初始值
            first_start = adjusted_start
            last_end = adjusted_end
        if adjusted_start < first_start:    # 基于初始值，再寻找起点和终点
            first_start = adjusted_start
        if adjusted_end > last_end: 
            last_end = adjusted_end

        inter_region_start = prev_end + 1
        inter_region_end = adjusted_start - 1

        if inter_region_start < inter_region_end-1:
            ctg_name = f"ctg{chr}.{ctg_counter:d}"
            final_data.append((ctg_name, inter_region_start, inter_region_end, 'inter_region', "-"))
            ctg_counter += 1

        #final_data.append((fragment_id, adjusted_start, adjusted_end, direction))  # Keep the original direction
        final_data.append((fragment_id, start, end, direction))  # Keep the original direction
        prev_end = adjusted_end

    if last_end > genome_length:
        print(f"ERROR: The position of last repeat exceeds the 'chrom{chr}.fasta' length {genome_length}. Please confirm again.")
        sys.exit(1)

    if prev_end < genome_length:
        final_data.append((f"ctg{chr}.{ctg_counter:d}", prev_end + 1, genome_length, "inter_region", "-"))

    return final_data

##########################################################################################################################################
def filter_ctg_rows(eight_col_df):
    # 过滤掉以"ctg"开头的行
    return eight_col_df[~eight_col_df['fragment_id'].str.startswith('ctg')]

##########################################################################################################################################
def check_fasta_sequence_length(fasta_file_path):  
    #先检查是不是fasta格式
    file_format = check_file_format_efficient(fasta_file_path)
    if file_format != 'FASTA':
        print(f"ERROR: The infile {fasta_file_path} is not in FASTA format.")
        sys.exit(1)
        
    sequences = []  
    with open(fasta_file_path, 'r') as fasta_file:  
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
        print(f"ERROR: Only one sequence should be in: {fasta_file}.")
        sys.exit(1) 
    elif len(sequences) == 0:  
        print(f"ERROR: Only one sequence should be in: {fasta_file}.")
        sys.exit(1)  
    else:  
        # 返回序列的长度  
        return len(sequences[0]) 

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

################################################################################
def save_5CT_for_mapping(combined_new_data1, label, args):
    # 将combined_new_data保存为5列表，用于绘图
    processed_data_with_paired_id = []
    
    for row in combined_new_data1:
        if row[0] == 'fragment_id':
            # Change the header row
            row = ['fragment_id', 'start', 'end', 'type', 'paired_id']
        else:
            # Add '-' to the other rows for 'paired_id'
            row.append('-')
        processed_data_with_paired_id.append(row)
    
    file_name = f"{args.output_file}_DR_{label}_2to2_5CT.tsv"    # 中间结果，含有project_id
    # Save the data to a TSV file
    with open(file_name, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(processed_data_with_paired_id)

################################################################################
def form_8CT_from_5CT(combined_new_data1, chrom1_file, args):
    ####  形成新的8CT
    #对new_data_1_calibrate1, new_data_1_calibrate2等拆开的5CT分别形成8CT
    #根据5CT的第一列id，在8CT中找到配对的id，然后找到的id再在5CT中寻找，能找到就配对，否则就表明没有配对id

    # 分离列名和数据  
    column_names = combined_new_data1[0]
    data_rows = combined_new_data1[1:]  
  
    data_rows = [row[:4] for row in data_rows]  # 只保留每行的前四个元素  
    combined_new_data1 = pd.DataFrame(data_rows, columns=column_names[:4])  # 只取前4个列名，因为数据只有4列  

    # 创建一个空的DataFrame来存储8CT结果，并预设列名  
    eight_column_table = pd.DataFrame(columns=['fragment_id', 'start', 'end', 'direction', 'paired_id', 'paired_start', 'paired_end', 'paired_direction'])  

    # 遍历 combined_new_data1 中的每个ID  
    for index, row in combined_new_data1.iterrows():  
        current_id = row['fragment_id']  
        paired_id = None  
        paired_data = [np.nan, np.nan, np.nan, np.nan]  # 初始化配对ID的数据为空  

        # 在 chrom1_file 中查找当前ID的配对ID  
        for uid_index, uid_row in chrom1_file.iterrows():  
            if current_id == uid_row['fragment_id']:  
                paired_id = uid_row['paired_id']  
                break  
            elif current_id == uid_row['paired_id']:  
                paired_id = uid_row['fragment_id']  
                break  

        # 提取当前ID的前四列数据  
        current_data = row[0:4].tolist()  # 假设前四列是按顺序排列的'fragment_id', 'start', 'end', 'direction'  

        # 如果找到了配对ID，并且在 combined_new_data1 中也存在，则提取配对ID的前四列数据  
        if paired_id is not None and paired_id in combined_new_data1['fragment_id'].values:  
            paired_row = combined_new_data1[combined_new_data1['fragment_id'] == paired_id].iloc[0]  
            paired_data = paired_row[0:4].tolist()  # 假设前四列是按顺序排列的  

        # 将当前ID和配对ID的数据合并到eight_column_table中  
        new_row = {
            'fragment_id': current_data[0],
            'start': current_data[1],
            'end': current_data[2],
            'direction': current_data[3],
            'paired_id': paired_data[0] if pd.notnull(paired_data[0]) else None,
            'paired_start': paired_data[1] if pd.notnull(paired_data[1]) else None,
            'paired_end': paired_data[2] if pd.notnull(paired_data[2]) else None,
            'paired_direction': paired_data[3] if pd.notnull(paired_data[3]) else None
            }
            
        eight_column_table = pd.concat([eight_column_table, pd.DataFrame([new_row])], ignore_index=True)

    # 注意：return 语句现在位于循环外部，确保所有行都被处理后才返回结果  
    return eight_column_table

##########################################################################################################################################
def split_data_and_fasta_by_paired_direct_repeats(paired_ids, unique_ids_with_direction_1, unique_ids_with_direction_2, chrom1_file, chrom1_filter, chrom2_file, chrom2_filter, long_genome, short_genome, args):

    # 根据paired_ids，形成重组基因组对应的5CT，5CT用于绘制基因组图谱
    # 根据paired_ids，再形成8CT
    selected_pairs = []
    error_count = 0
    auto_processed = False
    
    if len(paired_ids) == 0:
        print("ATTENTION: No Direct Repeat pairs found. Please check the repeat type carefully.")
        exit(1)
    
    while True:
        if args.auto and not auto_processed:
            # Automatically add all possible pairs and proceed
            selected_pairs.extend(paired_ids)
            print("Automatically processing all possible pairs...")
            time.sleep(2)
            auto_processed = True
            break
        else:
    
            user_input = input("Enter paired Direct Repeat IDs (e.g., RP1a RP1c, case sensitive), 'a/all' to process all, or press 'Enter' to exit:").strip()
    
            if user_input == '':
                if not selected_pairs:
                    print("\nNo IDs entered. Program was terminated by user.")
                    sys.exit(1)
                else:
                    print("Processing entered IDs...\n")
                    break
            elif user_input.lower() in ['a', 'all']:
                # Add all possible pairs to selected_pairs
                selected_pairs.extend(paired_ids)
                print("All possible pairs have been added.")
                break
            else:
                input_ids = sorted(user_input.split())
                if len(input_ids) == 2:
                    match1 = re.match(r'([A-Za-z]+\d+)', input_ids[0])
                    match2 = re.match(r'([A-Za-z]+\d+)', input_ids[1])
    
                    if match1 and match2:
                        prefix1 = match1.group(1)
                        prefix2 = match2.group(1)
    
                        if prefix1 == prefix2 and input_ids[0] != input_ids[1]:
                            direction1 = unique_ids_with_direction_1.loc[unique_ids_with_direction_1['fragment_id'] == input_ids[0], 'direction'].values
                            direction2 = unique_ids_with_direction_2.loc[unique_ids_with_direction_2['fragment_id'] == input_ids[1], 'direction'].values
                            if not direction1.size and not direction2.size:
                                direction1 = unique_ids_with_direction_1.loc[unique_ids_with_direction_1['fragment_id'] == input_ids[1], 'direction'].values
                                direction2 = unique_ids_with_direction_2.loc[unique_ids_with_direction_2['fragment_id'] == input_ids[0], 'direction'].values
                                input_ids[0],input_ids[1] = input_ids[1],input_ids[0]       # 调整位置，以使后面在数据集中能找到id
                                
                            if direction1.size > 0 and direction2.size > 0 and direction1[0] == direction2[0]:
                                selected_pairs.append(input_ids)
                                if input_ids in paired_ids:
                                    paired_ids.remove(input_ids)
                                else:
                                    print("Invalid or duplicate IDs.")
                                    error_count += 1
                                if len(paired_ids) == 0:     # 当用户输入的重复序列对 等于 用户提供的信息中的序列对数量 则自动终止
                                   print("All possible pairs have been added.")
                                   break 
                            else:
                                print("Invalid input. The two IDsshould be from different chromosomes.")
                                error_count += 1
                        else:
                            lprint("Invalid input. Please enter two IDs with the same prefix and different suffixes, or press 'Enter' to skip.")
                            error_count += 1
                    else:
                        print("Invalid input. Please enter two IDs with the same prefix and different suffixes, or press 'Enter' to skip.")
                        error_count += 1
                else:
                    print("Invalid input. Please enter two IDs with the same prefix and different suffixes, or press 'Enter' to skip.")
                    error_count += 1
    
            if error_count >= 3:
                print("ERROR: Too many incorrect attempts. Program will now exit.")
                sys.exit(1)
    
################################################################################
    # 在 while 循环之外处理 selected_pairs
    for pair in selected_pairs:
        id1, id2 = pair
        new_data_1_calibrate1, new_data_1_calibrate2 = split_data_by_direct_repeat_line(unique_ids_with_direction_1, id1)    
        # unique_ids_with_direction_1：从输入的8CT转换来的3CT，ctg和重复单元用长度表示了，此步是根据id1拆分这个3CT，可以直接绘图
        new_data_2_calibrate1, new_data_2_calibrate2 = split_data_by_direct_repeat_line(unique_ids_with_direction_2, id2)
        
        seq_before_end1, seq_after_end1, seq_before_end_length1, seq_after_end_length1, start1, end1 = split_fasta_by_direct_repeat_line(chrom1_filter, id1, long_genome)
        # chrom1_filter：从输入的8CT转换来的5CT，ctg和重复单元用start和end表示的，是unique_ids_with_direction_1之前的数据
        seq_before_end2, seq_after_end2, seq_before_end_length2, seq_after_end_length2, start2, end2 = split_fasta_by_direct_repeat_line(chrom2_filter, id2, short_genome)
        # 将切分的fasta序列，以交叉互换的方式，形成两条新的fasta
        new_fasta1 = seq_before_end1 + seq_after_end2
        new_fasta2 = seq_before_end2 + seq_after_end1
        
        filename1 = f"{args.output_file}_DR_{id1}_{id2}_c1_2to2.fasta" 
        with open(filename1, 'w') as fasta_file:                      # 根据从两个8CT形成的一个5CT，产生5CT对应的fasta  f"{id1}_{id2}_c1",
            fasta_file.write(f">DR_{id1}_{start1}_{end1}_Linear_{seq_before_end_length1 + seq_after_end_length2}_2to2\n")  
            # Split the sequence into lines of 100 characters for readability  
            for i in range(0, len(new_fasta1), 100):  
                fasta_file.write(str(new_fasta1[i : i+100]) + "\n") 
                
        filename2 = f"{args.output_file}_DR_{id1}_{id2}_c2_2to2.fasta"
        with open(filename2, 'w') as fasta_file:                      # 根据从两个8CT形成的一个5CT，产生5CT对应的fasta
            fasta_file.write(f">DR_{id2}_{start2}_{end2}_Linear_{seq_before_end_length2 + seq_after_end_length1}_2to2\n")  
            # Split the sequence into lines of 100 characters for readability  
            for i in range(0, len(new_fasta2), 100):  
                fasta_file.write(str(new_fasta2[i : i+100]) + "\n") 
        
        # 交叉互换后，组合形成新的5CT，用于绘图
        combined_new_data1 = pd.concat([new_data_1_calibrate1, new_data_2_calibrate2], ignore_index=True)
        combined_new_data1 = convert_length_to_start_end(combined_new_data1)
        combined_new_data2 = pd.concat([new_data_2_calibrate1, new_data_1_calibrate2], ignore_index=True)
        combined_new_data2 = convert_length_to_start_end(combined_new_data2)
        
        # 将combined_new_data保存为5列表，用于绘图
        save_5CT_for_mapping(combined_new_data1, f"{id1}_{id2}_c1", args)
        save_5CT_for_mapping(combined_new_data2, f"{id1}_{id2}_c2", args)
        
        # 从拆开的5CT出发，产生8CT，即绘制的图谱对应的8CT。
        eight_column_table1 =  form_8CT_from_5CT(combined_new_data1, chrom1_file, args)
        eight_column_table2 =  form_8CT_from_5CT(combined_new_data2, chrom2_file, args)
        
        # 去重产生的8CT，并保存
        drop_duplicates_and_save_8CT(eight_column_table1, f"{id1}_{id2}_c1", args)
        drop_duplicates_and_save_8CT(eight_column_table2, f"{id1}_{id2}_c2", args)
        
################################################################################
def drop_duplicates_and_save_8CT(merged_file, label, args):

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
    # 对含有空值的行进行去重
    #nan_rows = drop_duplicates_in_nan_rows(nan_rows)
    #nan_rows = process_unpaired_columns(nan_rows)
    
    # 最终的DataFrame，包含含有空值的行和去重后的不含空值的行  
    final_merged_file = pd.concat([nan_rows, unique_non_nan_rows], ignore_index=True)
    merged_file = final_merged_file.iloc[:, :8]
    
    # 步骤2: 保存merged_file为 TSV 文件 
    file_name = f"{args.output_file}_DR_{label}_2to2_map.tsv" 
    # 首先将 DataFrame 转换为列表的列表  
    data_to_write = merged_file.values.tolist()  
    with open(file_name, 'w', newline='') as file:  
        writer = csv.writer(file, delimiter='\t')
        # 首先写入表头
        writer.writerow(merged_file.columns)
        # 再写入数据
        writer.writerows(merged_file.values.tolist())

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
##########################################################################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process data from external files.")
    parser.add_argument('-i', '--chrom1_file', required=True, help="Path to the chrom1 data file (8CT).")
    parser.add_argument('-j', '--chrom2_file', required=True, help="Path to the chrom2 data file (8CT).")
    parser.add_argument('-o', '--output_file', required=True, help="Prefix of the output data file.")
    parser.add_argument('-l', '--long_genome', required=True, help="Organelle genome1 sequence.")
    parser.add_argument('-s', '--short_genome', required=True, help="Organelle genome2 sequence.")
    parser.add_argument('-auto', action='store_true', help="Enable auto mode to automatically process data.")
    parser.add_argument('-g', '--log', required=True, help="If consider complementary chains.")

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)
        
    print("The following operation is to recombinate two linear chromosomes mediated by Direct Repeat.")
    print("Make sure the two chromosomes come from the same DNA strand, e.g., all from plus or minus strand!")
    time.sleep(2)
    
    chrom1_len = check_fasta_sequence_length(args.long_genome)
    chrom2_len = check_fasta_sequence_length(args.short_genome)  
    
    comp_ch_2to2_log = args.log

################################################################################
    # Read the 8-column data
    if os.path.isfile(f"{args.chrom1_file}"):
        chrom1_file = pd.read_csv(args.chrom1_file, sep='\t', header=0)
        chrom1_file = drop_duplicates_in_8CT(chrom1_file)
    else:
        print(f"ERROR: {args.chrom1_file} non-exist.")
        sys.exit(1)

    if os.path.isfile(f"{args.chrom2_file}"):
        chrom2_file = pd.read_csv(args.chrom2_file, sep='\t', header=0)
        chrom2_file = drop_duplicates_in_8CT(chrom2_file)
    else:
        print(f"ERROR: {args.chrom2_file} non-exist.")
        sys.exit(1)

    # Create a new DataFrame with the first four columns
    chrom1_filter_top = chrom1_file.iloc[:, :4].copy()
    # Create another DataFrame with the last four columns, renaming them to match the first four
    chrom1_filter_bottom = chrom1_file.iloc[:, 4:].copy()  # Original last four columns  
    chrom1_filter_bottom.columns = chrom1_filter_top.columns  # Rename columns to match chrom1_filter_top
    
    # Concatenate the two DataFrames vertically
    chrom1_filter = pd.concat([chrom1_filter_top, chrom1_filter_bottom], ignore_index=True)
    # Drop any rows with missing values (if needed)
    chrom1_filter.dropna(subset=['start', 'end', 'direction'], inplace=True)
    chrom1_filter.drop_duplicates(inplace=True)
    # 在转换为列表之前进行排序  
    # 添加一个辅助列，用于存储第二列和第三列中的较小值  
    chrom1_filter['start'] = chrom1_filter['start'].astype(int)  # 将 start 列转换为整数  
    chrom1_filter['end'] = chrom1_filter['end'].astype(int)  # 将 end 列转换为整数  
    chrom1_filter['sort_key'] = chrom1_filter.iloc[:, 1].where(chrom1_filter.iloc[:, 1] < chrom1_filter.iloc[:, 2], chrom1_filter.iloc[:, 2])  
    # 根据新添加的辅助列进行排序  
    chrom1_filter = chrom1_filter.sort_values(by='sort_key')  
    # 如果不需要辅助列，可以删除它  
    chrom1_filter = chrom1_filter.drop(columns=['sort_key'])  

    chrom1_filter = chrom1_filter.values.tolist()
    chrom1_filter = add_ctg(chrom1_filter, chrom1_len, 1)

    # Create a new DataFrame with the first four columns
    chrom2_filter_top = chrom2_file.iloc[:, :4].copy()
    # Create another DataFrame with the last four columns, renaming them to match the first four
    chrom2_filter_bottom = chrom2_file.iloc[:, 4:].copy()  # Original last four columns  
    chrom2_filter_bottom.columns = chrom2_filter_top.columns
    
    # Concatenate the two DataFrames vertically
    chrom2_filter = pd.concat([chrom2_filter_top, chrom2_filter_bottom], ignore_index=True)
    
    # Drop any rows with missing values (if needed)
    chrom2_filter.dropna(subset=['start', 'end', 'direction'], inplace=True)
    chrom2_filter.drop_duplicates(inplace=True)
    # 在转换为列表之前进行排序  
    # 添加一个辅助列，用于存储第二列和第三列中的较小值  
    chrom2_filter['start'] = chrom2_filter['start'].astype(int)  # 将 start 列转换为整数  
    chrom2_filter['end'] = chrom2_filter['end'].astype(int)  # 将 end 列转换为整数  
    chrom2_filter['sort_key'] = chrom2_filter.iloc[:, 1].where(chrom2_filter.iloc[:, 1] < chrom2_filter.iloc[:, 2], chrom2_filter.iloc[:, 2])  
    
    # 根据新添加的辅助列进行排序  
    chrom2_filter = chrom2_filter.sort_values(by='sort_key')  
    # 如果不需要辅助列，可以删除它   
    chrom2_filter = chrom2_filter.drop(columns=['sort_key'])  

    chrom2_filter = chrom2_filter.values.tolist()
    chrom2_filter = add_ctg(chrom2_filter, chrom2_len, 2)
    
    chrom1_filter = pd.DataFrame(chrom1_filter, columns=['fragment_id', 'start', 'end', 'direction', 'paired_id'])    # list转换为dataframe
    chrom2_filter = pd.DataFrame(chrom2_filter, columns=['fragment_id', 'start', 'end', 'direction', 'paired_id'])

    # Create a new data table with two columns: 'fragment_id', 'length', 'direction', 'paired_id'
    new_data1 = pd.DataFrame({'fragment_id': chrom1_filter['fragment_id'],
                              'length': np.where(chrom1_filter['end'] >= chrom1_filter['start'],
                                                 chrom1_filter['end'] - chrom1_filter['start'] + 1,
                                                 chrom1_filter['end'] - chrom1_filter['start'] - 1),
                              'direction': chrom1_filter['direction']
                             })
    new_data2 = pd.DataFrame({'fragment_id': chrom2_filter['fragment_id'],
                              'length': np.where(chrom2_filter['end'] >= chrom2_filter['start'],
                                                 chrom2_filter['end'] - chrom2_filter['start'] + 1,
                                                 chrom2_filter['end'] - chrom2_filter['start'] - 1),
                              'direction': chrom2_filter['direction']
                             })
                             
    # 找出在 new_data2 的 fragment_id 列中也出现在 new_data1 的 fragment_id 列中的值  
    duplicates = new_data2['fragment_id'].isin(new_data1['fragment_id'])  
    # 检查是否存在重复值  
    if duplicates.any():  
        print("ERROR: There are duplicate fragment_ids in two datasets, please confirm.")  
        sys.exit(1)  # 如果存在重复值，则停止程序运行  

    # 根据 new_data1 的 'fragment_id' , 'length'和 'direction' 列，去除重复的行
    unique_ids_with_direction_1 = new_data1[['fragment_id', 'length', 'direction']].drop_duplicates()

    # 根据 new_data2 的 'fragment_id' 和 'direction' 列，去除重复的行
    unique_ids_with_direction_2 = new_data2[['fragment_id', 'length', 'direction']].drop_duplicates()

    # 定义一个函数来提取 ID 的前缀和数字部分
    def extract_prefix_and_number(id):
        match = re.match(r'([A-Za-z]+)(\d+)', id)
        return match.group(1) + match.group(2) if match else None
    
    # 创建一个空列表来存储可以配对的 DR_ID
    paired_ids = []
    # 创建一个空列表，存储将一条染色体互补翻转(complementary flipped strand)后可以配对的 DR_ID
    paired_ids_CS = [] 
    
    # 遍历 chromosome1 中的每个 ID
    for index1, row1 in unique_ids_with_direction_1.iterrows():
        if row1['fragment_id'].startswith('ctg'):
            continue
        prefix_and_number1 = extract_prefix_and_number(row1['fragment_id'])
        direction1 = row1['direction']
    
        # 遍历 chromosome2 中的每个 ID
        for index2, row2 in unique_ids_with_direction_2.iterrows():
            if row2['fragment_id'].startswith('ctg'):
                continue
            prefix_and_number2 = extract_prefix_and_number(row2['fragment_id'])
            direction2 = row2['direction']
    
            # 检查前缀和数字部分是否相同，且方向相同
            if prefix_and_number1 == prefix_and_number2 and direction1 == direction2:
                # 如果可以配对，将它们添加到列表中
                paired_ids.append([row1['fragment_id'], row2['fragment_id']])
                
            # 检查前缀和数字部分是否相同，且方向相反
            if prefix_and_number1 == prefix_and_number2 and direction1 != direction2:
                # 如果可以配对，将它们添加到列表中
                paired_ids_CS.append([row1['fragment_id'], row2['fragment_id']])

################################################################################
    # 打印正向重复序列配对的 ID
    if paired_ids:
        print("The following IDs are from paired Direct Repeats based on your repeat info:")
        for i, pair in enumerate(paired_ids, start=1):
            print(f"{i}. {pair}")
        # 根据paired_ids，形成重组基因组对应的5CT，5CT用于绘制基因组图谱，# unique_ids_with_direction_1, unique_ids_with_direction_2, 包含['fragment_id', 'length', 'direction'] 三列内容
        split_data_and_fasta_by_paired_direct_repeats(paired_ids, unique_ids_with_direction_1, unique_ids_with_direction_2, chrom1_file, chrom1_filter, chrom2_file, chrom2_filter, args.long_genome, args.short_genome, args)
    
################################################################################
    # 打印反向重复序列配对的 ID，但是将其中一条染色体取互补链，反向重复序列对就可以变为正向重复序列对
    if paired_ids_CS and comp_ch_2to2_log:    # == True:
        print("The following IDs are from paired Inverted Repeats, which can become paired Direct Repeats when one of the provided chromosomes takes the complementary strand.")
        print(f"Take the complementary strand of the second chromosome {args.chrom2_file} and draw the possible genome rearrangement maps.")
        for i, pair in enumerate(paired_ids_CS, start=1):
            print(f"{i}. {pair}")
            
        # 对于paired_ids_CS，需要对一条染色体的5CT进行翻转处理，仅对unique_ids_with_direction_2进行处理
        unique_ids_with_direction_2 = complementary_strand_5CT_and_fasta(unique_ids_with_direction_2)
        short_genome = save_complement_sequences(args.short_genome, 'output_temp_file.fasta') 
        split_data_and_fasta_by_paired_direct_repeats(paired_ids_CS, unique_ids_with_direction_1, unique_ids_with_direction_2, chrom1_file, chrom1_filter, chrom2_file, chrom2_filter, args.long_genome, short_genome, args)
        os.remove('output_temp_file.fasta')    # 删除临时的fasta文件

