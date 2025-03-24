#!/usr/bin/env python

import pandas as pd
import sys
import re
import argparse
import warnings
import numpy as np

# Filter out the specific FutureWarning
warnings.filterwarnings('ignore', category=FutureWarning)

#####################################################################################################################################################################
def process_repeat_type(combined_df):
    # 准备5列的DataFrame
    unique_fragment_ids = combined_df['fragment_id'].dropna().unique()  # Exclude NaN values
    
    # 定义一个正则表达式模式来匹配以字母开头，字母+数字的形式
    pattern = re.compile(r'^[A-Za-z]+[0-9]{1,}[A-Za-z]{1,}$')

    # 检查 unique_fragment_ids 的格式
    for fragment_id in unique_fragment_ids:
        if not pattern.match(fragment_id):
            print(f"ERROR: The format of Repeat name '{fragment_id}' is invalid.")
            print(f"ATTNTION: It must start with letters and be followed by one or more numbers (e.g., 'A1', 'AB2').")
            print("ATTNTION: The same Repeat unit can be distinguished by different letter suffixes, which is optional.")
            sys.exit(1)
    
    # 检查 unique_fragment_ids 的数量,超过30，则停止程序。
    if len(unique_fragment_ids) > 30:
        print()
        print("ERROR: Repeat is more than 30. Too many duplicate units to display, currently accepting 30 Repeats.\n")
        sys.exit(1)
    
    five_col_df = pd.DataFrame(columns=["fragment_id", "start", "end", "type", "paired_id"])

    for fragment_id in unique_fragment_ids:
        fragment_rows = combined_df[combined_df['fragment_id'] == fragment_id]

        direct_repeat_ids = []
        inverted_repeat_ids = []

        for _, row in fragment_rows.iterrows():
            if pd.notna(row['paired_id']):  # Check if 'paired_id' is not NaN
                if row['direction'] == row['paired_direction']:
                    direct_repeat_ids.append(row['paired_id'])
                else:
                    inverted_repeat_ids.append(row['paired_id'])

        if len(direct_repeat_ids) == 0 and len(inverted_repeat_ids) == 0:
            type = "-"
            paired_ids = "-"
        elif len(fragment_rows) == len(direct_repeat_ids):
            type = "direct_repeat"
            paired_ids = ', '.join(direct_repeat_ids)
        elif len(fragment_rows) == len(inverted_repeat_ids):
            type = "inverted_repeat"
            paired_ids = ', '.join(inverted_repeat_ids)
        else:
            type = "two_types"
            paired_ids = ', '.join(direct_repeat_ids) + " | " + ', '.join(inverted_repeat_ids)

        first_row = fragment_rows.iloc[0]
        new_row = pd.DataFrame({
            "fragment_id": [fragment_id],
            "start": [first_row['start']],
            "end": [first_row['end']],
            "type": [type],
            "paired_id": [paired_ids]
        })

        if not new_row.isnull().values.all():  # Check if new_row is not all-NA
            five_col_df = pd.concat([five_col_df, new_row], ignore_index=True)

    return five_col_df

#####################################################################################################################################################################
def create_ctg_final_tsv(processed_data, args):
    """Create the final data structure and save it to a file."""
    processed_data.sort(key=lambda x: x[1])  # Sort by start position
    prev_end = 0
    ctg_counter = 1
    final_data = []

    for data in processed_data:
        fragment_id, start, end, repeat_type, paired_id = data

        # Ensure start is always less than or equal to end
        adjusted_start = min(start, end)
        adjusted_end = max(start, end)

        # Calculate the intergenic region before the current fragment
        inter_region_start = prev_end + 1
        inter_region_end = adjusted_start - 1 if start <= end else end - 1

        # Add the intergenic region to final_data if it's valid
        if inter_region_start <= inter_region_end:
            ctg_name = f"ctg{ctg_counter:02d}"
            final_data.append((ctg_name, inter_region_start, inter_region_end, 'inter_region', '-'))
            ctg_counter += 1

        # Add the current fragment to final_data
        final_data.append((fragment_id, start, end, repeat_type, paired_id))
        prev_end = adjusted_end

    # Add the last intergenic region if needed
    if prev_end < args.genome_length:
        ctg_name = f"ctg{ctg_counter:02d}"
        final_data.append((ctg_name, prev_end + 1, args.genome_length, 'inter_region', '-'))

    return final_data

#####################################################################################################################################################################
def create_combined_dataframe(data):
    # 创建原始DataFrame，提取前八列
    df = pd.DataFrame(data).iloc[:, :8]
    df.columns = ["fragment_id", "start", "end", "direction", "paired_id", "paired_start", "paired_end", "paired_direction"]

    # 创建翻转后的DataFrame，对调前四列和后四列的位置
    swapped_df = df.copy()
    swapped_df.columns = ["paired_id", "paired_start", "paired_end", "paired_direction", "fragment_id", "start", "end", "direction"]

    # 合并原始DataFrame和翻转后的DataFrame
    combined_df = pd.concat([df, swapped_df], ignore_index=True)

    return combined_df

#####################################################################################################################################################################
def main():
    parser = argparse.ArgumentParser(description="Process paired information and generate 5-column table.")
    parser.add_argument('-i', '--input', type=str, required=True, help='Input aggregate results with 17CT.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file prefix.')
    parser.add_argument('-l', '--genome_length', type=int, required=True, help='Mitogenome length.')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    file_path = args.input
    output_prefix = args.output

    try:
        data = pd.read_csv(file_path, sep='\t')

        data = create_combined_dataframe(data)

        processed_data = process_repeat_type(data)
        processed_data_list = processed_data.values.tolist()  # Convert DataFrame to list of tuples
        final_data = create_ctg_final_tsv(processed_data_list, args)

        output_file = f"{output_prefix}"
        pd.DataFrame(final_data, columns=['fragment_id', 'start', 'end', 'type', 'paired_id']).to_csv(output_file, sep='\t', index=False)

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)  # 终止程序

#####################################################################################################################################################################
if __name__ == "__main__":
    main()
