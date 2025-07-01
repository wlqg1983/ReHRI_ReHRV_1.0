#!/usr/bin/env python

import pandas as pd
import argparse
import sys, shutil
import os
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

#######################################################################################################################################################################
def split_and_rename_columns(df):
    split_columns = df.iloc[:, 0].str.split('_', expand=True)
    df = pd.concat([split_columns, df.iloc[:, 1:]], axis=1)
    column_names = ['split_1', 'split_2', 'split_3', 'split_4', 'direction', 'trim_len', 'seq_length', 'spanning_read_number', 'spanning_read_remark']
    df.columns = column_names[:df.shape[1]]  # Ensure the column names match up to the actual number of columns
    return df

#######################################################################################################################################################################
# Generate new columns based on the specified rules. 考虑正负链
def generate_new_columns(row):
    # Generate new columns based on specified rules, considering the chain direction.
    new_col1, new_col3, new_col4, new_col5 = row['split_1'], row['split_3'], row['split_4'], row['direction']   # Initialize with existing values as default
    # Define mappings for new values based on conditions
    if row['split_2'] == 'LD':
        new_col2 = 'LR'
    elif row['split_2'] == 'UR':
        new_col2 = 'UD'
    elif row['split_2'] == 'LU':
        new_col2 = 'UD' 
        new_col5 = 'minus' if row['direction'] == 'plus' else 'plus'
    elif row['split_2'] == 'RD':
        new_col2 = 'LR'
        new_col5 = 'minus' if row['direction'] == 'plus' else 'plus'
    return pd.Series([new_col1, new_col2, row['split_3'], row['split_4'], new_col5], index=['new_split_1', 'new_split_2', 'new_split_3', 'new_split_4', 'new_direction'])

#######################################################################################################################################################################
# 定义匹配函数
def match_row(row, lookup_dict):
    new_row_values = tuple(row[['new_split_1', 'new_split_2', 'new_split_3', 'new_split_4', 'new_direction']].values)
    return (row.name, lookup_dict.get(new_row_values)) if new_row_values in lookup_dict else None

# 使用ThreadPoolExecutor进行多线程处理
def find_matches(new_columns, df_main_split, threads):
    # 创建一个字典用于快速查找
    lookup_dict = {tuple(row[['split_1', 'split_2', 'split_3', 'split_4', 'direction']].values): idx for idx, row in df_main_split.iterrows()}
    matches = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_row = {executor.submit(match_row, row, lookup_dict): row for index, row in new_columns.iterrows()}
        for future in as_completed(future_to_row):
            result = future.result()
            if result is not None and result[1] is not None:
                matches.append(result)
    return matches

#######################################################################################################################################################################
# 并行删除文件夹的函数
def delete_folder(folder_path):
    try:
        if os.path.isdir(folder_path):
            shutil.rmtree(folder_path)
    except Exception as e:
        print(f"Error processing folder: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check correspondence of trimmed sequence.")
    parser.add_argument("-mp", "--map_path", help="Path of mapping results in the mainconfiguration.")
    parser.add_argument("-sr", "--subconfig_reads", required=True, help="The aggregated result file of subconfiguration.")
    parser.add_argument("-mr", "--mainconfig_reads", required=True, help="The aggregated result file of mainconfiguration.")
    parser.add_argument("-o", "--output", required=True, help="Path of output result file for the final combined data.")
    parser.add_argument("-rn", "--read_number", default=1, type=int, help="Spanning read number of repeat supporting mitogenome recombination.")
    parser.add_argument("-ml", "--meaningless_results", default='K', help="Remove the mapping results of those trimmed sequences with no repeat-spanning RELDs.")
    parser.add_argument("-t", '--threads', type=int, required=True, help='Threads for speeding up.')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # 优化数据加载
    try:
        df_main = pd.read_csv(args.mainconfig_reads, sep='\t', header=0)
        df_sub = pd.read_csv(args.subconfig_reads, sep='\t', header=0)
    except FileNotFoundError:
        print("ERROR: One or both of the input files were not found.")
        sys.exit(1)

    # Process columns and split
    df_main_split = split_and_rename_columns(df_main)
    df_sub_split = split_and_rename_columns(df_sub)

    columns = ['sub_trimmed_seq', 'sub_trimmed_seq_length(bp)', 'sub_spanning_read_number', 'spanning_read_scfg',
               'main_trimmed_seq', 'main_trimmed_seq_length(bp)', 'main_spanning_read_number', 'spanning_read_mcfg']
    output_file_path = f"{args.output}/recomb-supporting_paired_repeat_correspondence.tsv"

    # Generate new columns and add to df_sub_split
    new_columns = df_sub_split.apply(generate_new_columns, axis=1)

    # Step 1: Match rows in new_columns with rows in df_main_split based on the first five columns
    # 调用多线程匹配函数
    matches = find_matches(new_columns, df_main_split, args.threads)

    # Step 2: Read corresponding rows from df_sub and df_main based on the matched row indices
    matched_rows = []
    for sub_idx, main_idx in matches:
        sub_row = df_sub.iloc[sub_idx]
        main_row = df_main.iloc[main_idx]
        matched_rows.append(pd.concat([sub_row, main_row], axis=0))

    # Step 3: Concatenate the corresponding rows horizontally
    # Ensure the column names are as specified
    matched_df = pd.DataFrame(matched_rows, columns=columns)

    # Step 4: Save to file
    matched_df.to_csv(output_file_path, sep='\t', index=False)

    # Step 5: Convert meaningless_results to uppercase for case-insensitive comparison
    meaningless_results_flag = args.meaningless_results.upper()

    # After finding matches and before saving to the output file
    remaining_indices = set(df_main.index) - {main_idx for _, main_idx in matches}
    folders_to_delete = []
    for idx in remaining_indices:
        prefix = df_main.iloc[idx, 0]
        try:
            for folder in os.listdir(args.map_path):
                if folder.startswith(prefix):
                    folder_path = os.path.join(args.map_path, folder)
                    if meaningless_results_flag == 'D' and os.path.isdir(folder_path):
                        folders_to_delete.append(folder_path)
        except Exception as e:
            print(f"Error processing folder: {e}")

    # 使用多进程并行删除文件夹
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        executor.map(delete_folder, folders_to_delete)

    # 汇总的结果为空的时候，如何处理对应关系。
    if df_sub_split.empty or df_main_split.empty:
        pd.DataFrame(columns=columns).to_csv(output_file_path, sep='\t', index=False)


