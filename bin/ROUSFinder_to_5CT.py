#!/usr/bin/env python

import pandas as pd
import argparse, os, sys, time
from itertools import combinations
from concurrent.futures import ThreadPoolExecutor

import io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, write_through=True)

########################################################################################################################################################################
def convert_columns_to_int(grouped):
    '''对第二、第三、第四列进行整数转换'''
    for group in grouped:
        for row in group:
            # 确保 row 是列表而不是嵌套列表或其他类型
            if isinstance(row, list):
                try:
                    # 对第二、第三、第四列进行整数转换，假设这些列对应的是索引 1、2、3
                    row[1] = int(row[1])  # 转换 'length'
                    row[2] = int(row[2])  # 转换 'start'
                    row[3] = int(row[3])  # 转换 'end'
                except ValueError as ve:
                    print(f"ValueError encountered: {ve}. Row: {row}")
            else:
                print(f"Unexpected row format: {row}")
    
    return grouped

################################################################################
def sort_group_and_row(merged_groups, column_names):
    '''对各组以及每一组内的元素进行排序'''
    # Process each merged group
    processed_merged_groups = []
    for merged_group in merged_groups:
        # Assign the same 'fragment_id' to all rows in the merged group
        fragment_id_to_assign = merged_group[0][column_names.index('fragment_id')]
        merged_group = [row[:column_names.index('fragment_id')] + [fragment_id_to_assign] + row[column_names.index('fragment_id')+1:] for row in merged_group]

        # Step 1: Collect order_info information
        order_info = []
        for row in merged_group:
            is_swap_needed = row[2] > row[3]
            #is_swap_needed = int(row[2]) > int(row[3])
            if is_swap_needed:
                # Swap "start" and "end" values
                row[2], row[3] = row[3], row[2]
            order_info.append((is_swap_needed, row))

        # Step 2: Sort based on "start" and then "end" if "start" values are equal
        sorted_order = sorted(order_info, key=lambda x: (x[1][2], x[1][3]))

        # Step 3: Restore the original order
        merged_group_sorted = []
        for is_swap_needed, row in sorted_order:
            if is_swap_needed:
                # Swap "start" and "end" values
                row[2], row[3] = row[3], row[2]
            merged_group_sorted.append(row)    # Now, merged_group_sorted contains the rows in the desired order
        processed_merged_groups.append(merged_group_sorted)

################################################################################
    final_df_processed = []
    # Step 1: Extract min_value_first_row for each sublist's first row
    min_values_first_row = [min(group[0][2], group[0][3]) for group in processed_merged_groups]

    # Step 2: Sort min_values_first_row and get the corresponding indices
    sorted_indices = sorted(range(len(min_values_first_row)), key=lambda k: min_values_first_row[k])

    # Step 3: Sort processed_group based on the sorted indices
    for index in sorted_indices:
        processed_group_sorted = sorted(processed_merged_groups[index], key=lambda x: (min_values_first_row[index], x[2], x[3]))
        final_df_processed.append(processed_group_sorted)  # Now, final_df_processed contains the processed and sorted data from all the merged groups

################################################################################
    # Process each merged group
    processed_merged_groups_again = []
    for merged_group_again in final_df_processed:
        # Assign the same 'fragment_id' to all rows in the merged group
        fragment_id_to_assign = merged_group_again[0][column_names.index('fragment_id')]
        merged_group_again = [row[:column_names.index('fragment_id')] + [fragment_id_to_assign] + row[column_names.index('fragment_id')+1:] for row in merged_group_again]

        # Step 1: Collect order_info information
        order_info = []
        for row in merged_group_again:
            is_swap_needed = row[2] > row[3]
            #is_swap_needed = int(row[2]) > int(row[3])
            if is_swap_needed:
                # Swap "start" and "end" values
                row[2], row[3] = row[3], row[2]
            order_info.append((is_swap_needed, row))

        # Step 2: Sort based on "start" and then "end" if "start" values are equal
        sorted_order = sorted(order_info, key=lambda x: (x[1][2], x[1][3]))

        # Step 3: Restore the original order
        merged_group_sorted = []
        for is_swap_needed, row in sorted_order:
            if is_swap_needed:
                # Swap "start" and "end" values
                row[2], row[3] = row[3], row[2]
            merged_group_sorted.append(row)    # Now, merged_group_sorted contains the rows in the desired order
        processed_merged_groups_again.append(merged_group_sorted)

################################################################################
    # Convert to list of dictionaries and create DataFrame
    final_df_processed = [dict(zip(column_names, row)) for group in processed_merged_groups_again for row in group]
    final_processed = pd.DataFrame(final_df_processed)

    # Remove square brackets from 'start' and 'end' columns
    final_processed['start'] = final_processed['start'].apply(lambda x: tuple(x[0]) if isinstance(x, (list, tuple)) and x else x)
    final_processed['end'] = final_processed['end'].apply(lambda x: tuple(x[0]) if isinstance(x, (list, tuple)) and x else x)
    
    return final_processed

########################################################################################################################################################################
def preprocess_input_file(infile):
    """ Preprocess the calibrated input file from the ROUSFinder2.0 script. """
    # Combine the current working directory with the input file's relative path

    # Read the entire TSV file as a DataFrame, skipping the first row
    df = pd.read_csv(infile, delimiter='\t')
    
    # Check for NaN or missing values
    if df.isnull().values.any():
        print("The input file contains NaN or missing values.")
        print(df[df.isnull().any(axis=1)])  # sys.stderr.write rows with NaN
        sys.exit(1)
        
    """Check for unique occurrence of fragment_id."""
    unique_fragment_ids = df['fragment_id'].value_counts()
    to_delete = []    # 存储要删除的fragment_id
    for fragment_id, occurrences in unique_fragment_ids.items():
        if occurrences == 1:
            print(f"Repeat sequence unit {fragment_id} appears only once.\n"
                  f"Press 'Ctrl + C' to stop the program now. Or it will be deleted.")
            to_delete.append(fragment_id)    # 如果fragment_id只出现一次，则加入到删除列表中
            try:
                # 等待10秒以给用户取消的机会
                print(f"Or, the {fragment_id} will be delete automatically.")
                sys.stdout.flush()    # 确保消息立即打印
                time.sleep(10)
            except KeyboardInterrupt:
                print("Program terminated by the user.")
                sys.exit(1)
                
    # 删除那些只出现一次的fragment_id的行
    if to_delete:
        # 找到需要删除的行
        rows_to_delete = df[df['fragment_id'].isin(to_delete)]
        # 打印出需要删除的行
        print(f"warning: The 'fragment_id' in the following rows are unpaired.\nThe rows were deleted:\n{rows_to_delete}")
        time.sleep(5)
        # 删除这些行
        df = df[~df['fragment_id'].isin(to_delete)]
                
    # 对组间和组内排序
    column_names = ['fragment_id', 'length', 'start', 'end', 'direction']
    # Group by 'fragment_id' and convert to a nested list
    grouped = [group.reset_index(drop=True).values.tolist() for _, group in df.groupby('fragment_id')]
    grouped = convert_columns_to_int(grouped)               # 对中间三列进行整数转换，便于后续操作
    sorted_group_and_row = sort_group_and_row(grouped, column_names)

    return sorted_group_and_row

########################################################################################################################################################################
def change_id(df):
    # 数据处理
    current_id = None
    counter = 0
    label = "RP"
    new_id = {}

    for index, row in df.iterrows():
        fragment_id = row['fragment_id']
        if fragment_id != current_id:
            current_id = fragment_id
            counter += 1
            sub_counter = 97  # ASCII value of 'a'
        else:
            sub_counter += 1

        label_suffix = ''
        temp_sub_counter = sub_counter
        while temp_sub_counter > 122:
            label_suffix += chr((temp_sub_counter - 97) % 26 + 97)
            temp_sub_counter = 97 + (temp_sub_counter - 97) // 26 - 1

        new_label = label + str(counter) + chr(temp_sub_counter) + label_suffix
        df.at[index, 'fragment_id'] = new_label

    return df

#######################################################################################################################################################################
def check_tsv_format(input_file):
    try:
        # 尝试读取文件，不限制行数，这样即使文件行数少于5行也能正确处理
        df = pd.read_csv(input_file, delimiter='\t')
        # 检查列数是否为5
        if df.shape[1] != 5:
            return False
        return True
    except Exception as e:
        # 如果读取过程中出现任何错误，打印错误并假定文件格式不正确
        print(f"Incorrect file format: {e}")
        return False

#######################################################################################################################################################################
def process_pair(pair, df):
    # 找到每个fragment_id对应的行数据
    fragment_id1, fragment_id2 = pair
    row1 = df[df['fragment_id'] == fragment_id1].iloc[0]
    row2 = df[df['fragment_id'] == fragment_id2].iloc[0]
    processed_data = []

    # 对每一对fragment_id进行处理，然后分别存入processed_data
    for row, paired_id in [(row1, fragment_id2), (row2, fragment_id1)]:
        id = row['fragment_id']
        start, end, direction = row['start'], row['end'], row['direction']
        original_order = start < end
        if not original_order:
            start, end = end, start
        # 将每一行的信息，连同其配对的fragment_id，存入processed_data
        processed_data.append((id, start, end, direction, original_order, paired_id))

    return processed_data
    
#######################################################################################################################################################################
# 使用多线程并行处理配对
def parallel_process_pairs(all_pairs, df, threads):
    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = executor.map(lambda pair: process_pair(pair, df), all_pairs)
    
    # 合并所有线程的结果
    processed_data = [item for sublist in results for item in sublist]
    return processed_data

########################################################################################################################################################################
def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='Process genomic data from a TSV file.')
    parser.add_argument('-i', '--input_file', type=str, required=True, help='Path to the input TSV file containing genomic data.')
    parser.add_argument('-l', '--genome_length', type=int, required=True, help='Total genome length as an integer.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output file prefix for the generated output file.')
    parser.add_argument('-t', '--threads', type=int, help='Threads for speeding up.')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)
    
    # 检查文件格式是否为5CT
    is_valid_format = check_tsv_format(args.input_file)
    if not is_valid_format:
        print(f"Please check if the file format of {args.input_file} is TSV and should be 5 columns.")
        sys.exit(1)

    # 读取TSV文件
    df = preprocess_input_file(args.input_file)
    # 更改fragment_id的前缀为RP
    df = change_id(df)

################################################################################
    #形成配对
    # 首先，提取 fragment_id 的前缀用于分组
    df['prefix'] = df['fragment_id'].str.extract(r'(RP\d+)')

    # 创建一个空列表来存储所有可能的配对
    all_pairs = []

    # 按前缀分组
    for prefix, group in df.groupby('prefix'):
        # 获取当前组的所有 fragment_id
        fragment_ids = group['fragment_id'].tolist()

        # 检查分组大小，确保每个分组至少有两个元素
        if len(fragment_ids) < 2:
            raise ValueError(f"Group {prefix} does not have enough items for pairing. Each group must have at least 2 items.")

        # 生成当前分组的所有可能的两两配对组合，并将它们添加到 all_pairs 列表中
        all_pairs.extend(list(combinations(fragment_ids, 2)))
    
    # 采用多进程，对每一配对的重复单元进行处理
    processed_data = parallel_process_pairs(all_pairs, df, args.threads)

################################################################################
    # 标记重复序列类型并插入间区
    final_data = []
    processed_data.sort(key=lambda x: x[1])  # 按起点排序
    prev_end = 0
    ctg_counter = 1
    
    # 初始化 directions 字典, 用于确定重复序列的类型
    directions = {}
    # 填充 directions 字典
    for index, row in df.iterrows():
        fragment_id = row['fragment_id']
        direction = row['direction']
        directions[fragment_id] = direction

    for data in processed_data:
        fragment_id, start, end, direction, original_order, paired_id = data

        # 确保start总是小于等于end
        adjusted_start = min(start, end)
        adjusted_end = max(start, end)

        # 如果是原始顺序（即非反向重复序列），则使用adjusted_start，否则使用adjusted_end + 1
        inter_region_start = prev_end + 1
        inter_region_end = adjusted_start - 1 if original_order else adjusted_end

        if inter_region_start <= inter_region_end:
            ctg_name = f"ctg{ctg_counter:02d}"
            final_data.append((ctg_name, inter_region_start, inter_region_end, 'inter_region', '-'))
            ctg_counter += 1

        # 重复序列类型
        repeat_type = 'inverted_repeat' if directions[fragment_id] != directions[paired_id] else 'direct_repeat'  #根据direction是否相同来判断重复序列类型

        # 恢复原始顺序
        if not original_order:
            start, end = end, start

        final_data.append((fragment_id, start, end, repeat_type, paired_id))
        prev_end = adjusted_end

################################################################################
    def print_error_and_exit(message):  
        print(f"{message}\n")  
        sys.exit(1)  

    if prev_end > args.genome_length:  
        print_error_and_exit("Repeat sequences are located beyond the length of the mitogenome.")

    # 当 prev_end 正好等于基因组长度时，没有剩余的间区需要添加
    if prev_end < args.genome_length:
        # 添加最后的间区
        final_data.append(("ctg{:02d}".format(ctg_counter), prev_end + 1, args.genome_length, "inter_region", "-"))
    
    # 创建DataFrame并直接保存到文件
    output_df = pd.DataFrame(final_data, columns=['fragment_id', 'start', 'end', 'type', 'paired_id'])
    output_df.to_csv(f"{args.output_file}_rep_5CT.tsv", sep='\t', index=False)
    
########################################################################################################################################################################
if __name__ == "__main__":
    main()


