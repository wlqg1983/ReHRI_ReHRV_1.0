#!/usr/bin/env python

import pandas as pd
import argparse, sys
from concurrent.futures import ThreadPoolExecutor

#######################################################################################################################################################################
def check_tsv_format(input_file):
    column_names = ['fragment_id', 'length', 'start', 'end', 'direction']
    try:
        # 尝试读取文件，不限制行数，这样即使文件行数少于5行也能正确处理
        df = pd.read_csv(input_file, delimiter='\t', header=None, names=column_names, skiprows=1)
        # 检查列数是否为5
        if df.shape[1] != 5:
            return False
        return True
    except Exception as e:
        # 如果读取过程中出现任何错误，打印错误并假定文件格式不正确
        print(f"Incorrect file format: {e}")
        return False

#######################################################################################################################################################################
def merge_group(group, merged_groups, flexibility, column_names):
    merge_condition_met = False
    for merged_group in merged_groups:
        for i in range(len(group)):
            for j in range(len(merged_group)):
                if (
                    abs(group[i][2] - merged_group[j][column_names.index('start')]) <= flexibility and
                    abs(group[i][3] - merged_group[j][column_names.index('end')]) <= flexibility
                ):
                    merge_condition_met = True
                    break
            if merge_condition_met:
                break
        if merge_condition_met:
            merged_group.extend(group)
            break
    if not merge_condition_met:
        merged_groups.append(group)
    return merged_groups
    
#######################################################################################################################################################################
def main(input_file, flexibility, output_prefix):
    # Define column names
    column_names = ['fragment_id', 'length', 'start', 'end', 'direction']
    
    is_valid_format = check_tsv_format(input_file)
    if not is_valid_format:
        print(f"Please check if the file format of {input_file} is TSV and should be 5 columns.")
        sys.exit(1)

    # Read data, skipping the first row, and specify column names
    df = pd.read_csv(input_file, delimiter='\t', header=None, names=column_names, skiprows=1)
    df.drop_duplicates(subset=['start', 'end'], inplace=True)  # Remove duplicates based on 'start', 'end'

    # Group by 'fragment_id' and convert to a nested list
    grouped = [group.reset_index(drop=True).values.tolist() for _, group in df.groupby('fragment_id')]

    # Merge groups
    merged_groups = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for group in grouped:
            futures.append(executor.submit(merge_group, group, merged_groups, flexibility, column_names))
        for future in futures:
            merged_groups = future.result()
            
    sort_group_and_row(merged_groups, column_names, output_prefix)

###################################################################################################################################################################
def sort_group_and_row(merged_groups, column_names, output_prefix):
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

    # Remove duplicates in the processed DataFrame
    final_processed.drop_duplicates(subset=['start', 'end'], inplace=True)

    # Check if a group has only one row and remove it
    final_processed = final_processed.groupby('fragment_id').filter(lambda x: len(x) > 1)

    # Save the processed calibrated results
    output_file_processed = f"{output_prefix}_rep_calibration_table.txt"
    final_processed.to_csv(output_file_processed, index=False, sep='\t')
    
    return final_processed
    
##############################################################################################################################################################
if __name__ == '__main__':
    # Argument parser
    parser = argparse.ArgumentParser(description='Calibrate ROUSFinder results.')
    parser.add_argument('-i', '--input', required=True, help='*_rep_table.txt of ROUSFinder results.')
    parser.add_argument('-fl', '--flexibility', type=int, default=5, help='Flexibility of merging overlapping repeats, default=5.')
    parser.add_argument('-o', '--output', required=True, help='Output file prefix for the calibrated results.')
    parser.add_argument("-t", '--threads', type=int, required=True, help='Threads for speeding up.')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # Run the main function
    main(args.input, args.flexibility, args.output)

