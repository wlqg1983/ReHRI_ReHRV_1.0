#!/usr/bin/env python

import argparse
import pandas as pd

#############################################################################################################################################################
def preprocess_sub_data(input_file_path, read_number):
    # 读取数据
    df = pd.read_csv(input_file_path, sep='\t', header=0)
    
    # 添加新列，基于条件设置值
    df['spanning_read_scfg'] = df.apply(lambda row: 'sufficient' if row.iloc[2] >= int(read_number) else 'insufficient', axis=1)
    
    # 将包含新列的DataFrame写回原始文件，原始文件中添加spanning_read_scfg列
    df.to_csv(input_file_path, sep='\t', index=False)

    # 读取spanning reads数量充足的行
    df_filtered = df[df['spanning_read_scfg'] == 'sufficient'] #.drop(columns=['spanning_read_scfg'])
    
    return df_filtered

#############################################################################################################################################################
def preprocess_main_data(input_file_path, read_number):
    # 读取数据
    df = pd.read_csv(input_file_path, sep='\t', header=0)
    
    # 添加新列，基于条件设置值
    df['spanning_read_mcfg'] = df.apply(lambda row: 'sufficient' if row.iloc[2] >= int(read_number) else 'insufficient', axis=1)
    
    # 将包含新列的DataFrame写回原始文件，原始文件中添加spanning_read_mcfg列
    df.to_csv(input_file_path, sep='\t', index=False)
    
    return df
    
######################################################################################################################################################################
def create_empty_output_file(filename, message):
    with open(filename, 'w', encoding='utf-8') as file:
        file.write(message)
        
######################################################################################################################################################################
def split_and_rename_columns(df): #将第一列按照"_"拆分后形成新的表，共9列，并重命名每一列的表头
    split_columns = df.iloc[:, 0].str.split('_', expand=True)
    df = pd.concat([split_columns, df.iloc[:, 1:]], axis=1)
    df.rename(columns={0: 'split_1', 1: 'split_2', 2: 'split_3', 3: 'split_4', 4: 'direction', 5: 'trim_len', 6: 'seq_length', 7: 'spanning_read_number', 8: 'spanning_read_remark'}, inplace=True)
    return df
    
######################################################################################################################################################################
def find_matching_element(pair_types, sorted_pair_types):
    # Check if pair_types is an element in sorted_pair_types
    return next((element for element in sorted_pair_types if pair_types == {element[0], element[1]}), None)
    
######################################################################################################################################################################
def create_paired_data_sub(df):
    df = split_and_rename_columns(df)
    paired_groups = []  # Save groups with paired fragment_id
    unpaired_groups = []  # Save groups without paired fragment_id

    pair_lookup = {'LU': 'RD', 'RD': 'LU', 'LD': 'UR', 'UR': 'LD'}
    pair_lookup_list = list(pair_lookup.items())

    # Split into two major groups: DR and IR
    for group_name in ['DR', 'IR']:
        group_df = df[df.iloc[:, 0] == group_name]

        # Ensure that there are at least 4 columns in the DataFrame
        if group_df.shape[1] >= 4:
            # Split into smaller groups based on the values in the third and fourth columns
            grouped = group_df.groupby([group_df.iloc[:, 2], group_df.iloc[:, 3]])

            for (split_3, split_4), sub_df in grouped:
                # Check if pairing is possible
                pair_types = set(sub_df.iloc[:, 1].unique())

                # Modify the check for matching elements using the helper function
                matching_element = find_matching_element(pair_types, pair_lookup_list)

                # Check if matching_element is not None before using it
                if matching_element is not None:
                    # Convert sorted_key to the desired format ('LU', 'RD')
                    formatted_sorted_key = ''.join(matching_element)
                    formatted_sorted_key = (formatted_sorted_key[:2], formatted_sorted_key[2:])

                    # Add the group to the corresponding list based on the pairing status
                    if matching_element == formatted_sorted_key:
                        paired_groups.append(sub_df.values.tolist())
                    else:
                        unpaired_groups.append(sub_df.values.tolist())
                else:
                    # Handle the case when matching_element is None
                    unpaired_groups.append(sub_df.values.tolist())

    return paired_groups, unpaired_groups

######################################################################################################################################################################
def create_paired_data_main(df):
    # 新的配对规则
    df = split_and_rename_columns(df)
    paired_groups = []  # Save groups with paired fragment_id
    unpaired_groups = []  # Save groups without paired fragment_id

    pair_lookup = {'LR': 'UD'}
    pair_lookup_list = list(pair_lookup.items())

    # Split into two major groups: DR and IR
    for group_name in ['DR', 'IR']:
        group_df = df[df.iloc[:, 0] == group_name]

        # Ensure that there are at least 4 columns in the DataFrame
        if group_df.shape[1] >= 4:
            # Split into smaller groups based on the values in the third and fourth columns
            grouped = group_df.groupby([group_df.iloc[:, 2], group_df.iloc[:, 3]])

            for (split_3, split_4), sub_df in grouped:
                # Check if pairing is possible
                pair_types = set(sub_df.iloc[:, 1].unique())

                # Modify the check for matching elements using the helper function
                matching_element = find_matching_element(pair_types, pair_lookup_list)

                # Check if matching_element is not None before using it
                if matching_element is not None:
                    # Convert sorted_key to the desired format ('LU', 'RD')
                    formatted_sorted_key = ''.join(matching_element)
                    formatted_sorted_key = (formatted_sorted_key[:2], formatted_sorted_key[2:])

                    # Add the group to the corresponding list based on the pairing status
                    if matching_element == formatted_sorted_key:
                        paired_groups.append(sub_df.values.tolist())
                    else:
                        unpaired_groups.append(sub_df.values.tolist())
                else:
                    # Handle the case when matching_element is None
                    unpaired_groups.append(sub_df.values.tolist())

    return paired_groups, unpaired_groups

######################################################################################################################################################################
def get_start_end(fragment_id, rep_5ct_file_path, rep_calibr_file_path):
    # Load rep_5ct_file and rep_calibr_file
    rep_5ct_df = pd.read_csv(rep_5ct_file_path, sep='\t', header=0)
    rep_calibr_df = pd.read_csv(rep_calibr_file_path, sep='\t', header=0)
    
    # Extract relevant information from rep_5ct_df based on fragment_id
    try:
        row_rep_5ct = rep_5ct_df[rep_5ct_df['fragment_id'] == fragment_id].iloc[0]
    except IndexError:
        print(f"No matching row found in rep_5ct_df for fragment_id: {fragment_id}")
        return None, None, None, None
        
    # Extract length, start, end, and direction from rep_calibr_df based on start and end from rep_5ct_df
    matching_rows = rep_calibr_df[(rep_calibr_df['start'] == row_rep_5ct['start']) & (rep_calibr_df['end'] == row_rep_5ct['end'])]
    if matching_rows.empty:
        #print(f"No matching rows found in rep_calibr_df for start={row_rep_5ct['start']}, end={row_rep_5ct['end']}")
        return None, None, None, None
        
    row_rep_calibr = matching_rows.iloc[0]
    
    return row_rep_calibr['length'], row_rep_calibr['start'], row_rep_calibr['end'], row_rep_calibr['direction']


######################################################################################################################################################################
def calculate_pair_ratio(sub_paired_groups, main_paired_groups, rep_5ct_file_path, rep_calibr_file_path): #, output_file):
    # 整合结果到一个列表
    result_list = []
    unique_sub_rows = set()  # Set to store unique sub_rows

    # 遍历 sub_paired_groups 和 main_paired_groups
    for sub_group in sub_paired_groups:
        for main_group in main_paired_groups:
            if sub_group[0][2] == main_group[0][2] and sub_group[0][3] == main_group[0][3]:
                ratio_result = calculate_ratio(sub_group, main_group)
                # Get unique key based on columns
                unique_key = (sub_group[0][0], sub_group[0][2], sub_group[0][3])
                # Check if the unique key is not in the set (not a duplicate)
                if unique_key not in unique_sub_rows:
                    unique_sub_rows.add(unique_key)  # Add the unique key to the set
                    # Get start and end for sub_row
                    sub_a_length, sub_a_start, sub_a_end, sub_a_direction = get_start_end(unique_key[1], rep_5ct_file_path, rep_calibr_file_path)
                    sub_b_length, sub_b_start, sub_b_end, sub_b_direction = get_start_end(unique_key[2], rep_5ct_file_path, rep_calibr_file_path)
                    if sub_a_length is None or sub_b_length is None or sub_a_start is None or sub_b_start is None:
                        continue

                    # Determine sequence type based on sub_row[0] sub_row
                    if sub_group[0][0] == "DR":
                        sequence_type = "direct_repeat" 
                    if sub_group[0][0] == "IR":
                        sequence_type = "inverted_repeat" 

                    # Append to result_list
                    if sequence_type == "direct_repeat":
                        result_list.append([
                            unique_key[1], sub_a_length, sub_a_start, sub_a_end, sub_a_direction, *ratio_result[:3], sequence_type,
                            unique_key[2], sub_b_length, sub_b_start, sub_b_end, sub_b_direction, *ratio_result[3:6]])
                    if sequence_type == "inverted_repeat":
                        result_list.append([
                            unique_key[1], sub_a_length, sub_a_start, sub_a_end, sub_a_direction, *ratio_result[6:9], sequence_type,
                            unique_key[2], sub_b_length, sub_b_start, sub_b_end, sub_b_direction, *ratio_result[9:12]])

    # Update the header
    header = ["fragment_id", "length", "start", "end", "direction", "plus_ratio(s/m)", "minus_ratio(s/m)", "combined_ratio", "type",
              "paired_id", "paired_length", "paired_start", "paired_end", "paired_direction", "paired_plus_ratio(s/m)", "paired_minus_ratio(s/m)", "paired_combined_ratio"]
    
    # Create a DataFrame
    result_df = pd.DataFrame(result_list, columns=header)

    return result_df

######################################################################################################################################################################
def calculate_unpair_ratio(sub_unpaired_groups, main_paired_groups, rep_5ct_file_path, rep_calibr_file_path):
    # 整合结果到一个列表
    i=1
    result_list = []
    unique_sub_rows = set()  # Set to store unique sub_rows

    # 遍历 sub_paired_groups 和 main_paired_groups
    for sub_group in sub_unpaired_groups:
        for main_group in main_paired_groups:
            if sub_group[0][2] == main_group[0][2] and sub_group[0][3] == main_group[0][3]:
                ratio_result = calculate_ratio(sub_group, main_group)
                # Get unique key based on columns
                unique_key = (sub_group[0][0], sub_group[0][2], sub_group[0][3])
                # Check if the unique key is not in the set (not a duplicate)
                if unique_key not in unique_sub_rows:
                    unique_sub_rows.add(unique_key)  # Add the unique key to the set
                    # Get start and end for sub_row
                    sub_a_length, sub_a_start, sub_a_end, sub_a_direction = get_start_end(unique_key[1], rep_5ct_file_path, rep_calibr_file_path)
                    sub_b_length, sub_b_start, sub_b_end, sub_b_direction = get_start_end(unique_key[2], rep_5ct_file_path, rep_calibr_file_path)

                    # Determine sequence type based on sub_row[0] sub_row
                    if sub_group[0][0] == "DR":
                        sequence_type = "direct_repeat" 
                    if sub_group[0][0] == "IR":
                        sequence_type = "inverted_repeat" 

                    # Append to result_list
                    if sequence_type == "direct_repeat" and sub_group[0][1] == "LD":
                        result_list.append([unique_key[1], sub_a_length, sub_a_start, sub_a_end, sub_a_direction, *ratio_result[:3], sequence_type, unique_key[2], sub_b_length, sub_b_start, sub_b_end, sub_b_direction, *ratio_result[3:6]])
                    if sequence_type == "direct_repeat" and sub_group[0][1] == "UR":
                        result_list.append([unique_key[2], sub_b_length, sub_b_start, sub_b_end, sub_b_direction, *ratio_result[3:6], sequence_type, unique_key[1], sub_a_length, sub_a_start, sub_a_end, sub_a_direction, *ratio_result[:3]])
                    if sequence_type == "inverted_repeat" and sub_group[0][1] == "LU":
                        result_list.append([unique_key[1], sub_a_length, sub_a_start, sub_a_end, sub_a_direction, *ratio_result[6:9], sequence_type, unique_key[2], sub_b_length, sub_b_start, sub_b_end, sub_b_direction, *ratio_result[9:12]])
                    if sequence_type == "inverted_repeat" and sub_group[0][1] == "RD":
                        result_list.append([unique_key[2], sub_b_length, sub_b_start, sub_b_end, sub_b_direction, *ratio_result[9:12], sequence_type, unique_key[1], sub_a_length, sub_a_start, sub_a_end, sub_a_direction, *ratio_result[6:9]])

    # Update the header
    header = ["fragment_id", "length", "start", "end", "direction", "plus_ratio(s/m)", "minus_ratio(s/m)", "combined_ratio", "type","paired_id", "paired_length", "paired_start", "paired_end", "paired_direction", "paired_plus_ratio(s/m)", "paired_minus_ratio(s/m)", "paired_combined_ratio"]
    
    # Create a DataFrame
    result_df = pd.DataFrame(result_list, columns=header)
    
    return result_df
    
######################################################################################################################################################################
def calculate_ratio(sub_group, main_group):
    # 初始化各个比例的变量
    ratio_DR_LD_C_LR_C_numerator, ratio_DR_LD_C_LR_C_denominator = 0, 0
    ratio_DR_LD_LR_numerator, ratio_DR_LD_LR_denominator = 0, 0
    ratio_DR_UR_C_UD_C_numerator, ratio_DR_UR_C_UD_C_denominator = 0, 0
    ratio_DR_UR_UD_numerator, ratio_DR_UR_UD_denominator = 0, 0
    ratio_IR_LU_UD_C_numerator, ratio_IR_LU_UD_C_denominator = 0, 0
    ratio_IR_LU_C_UD_numerator, ratio_IR_LU_C_UD_denominator = 0, 0
    ratio_IR_RD_C_LR_numerator, ratio_IR_RD_C_LR_denominator = 0, 0
    ratio_IR_RD_LR_C_numerator, ratio_IR_RD_LR_C_denominator = 0, 0
    ratio_DR_LD_LR = ratio_DR_LD_C_LR_C = ratio_DR_LD = ratio_DR_UR_UD = "N/A"
    ratio_DR_UR_C_UD_C = ratio_DR_UR = ratio_IR_LU_UD_C = ratio_IR_LU_C_UD = "N/A"
    ratio_IR_LU = ratio_IR_RD_LR_C = ratio_IR_RD_C_LR = ratio_IR_RD = "N/A"

    for sub_row in sub_group:
        if sub_row[1] == 'LD' and sub_row[4] == 'minus':
            ratio_DR_LD_C_LR_C_numerator = sub_row[7]
        if sub_row[1] == 'LD' and sub_row[4] == 'plus':
            ratio_DR_LD_LR_numerator = sub_row[7]
        if sub_row[1] == 'UR' and sub_row[4] == 'minus':
            ratio_DR_UR_C_UD_C_numerator = sub_row[7]
        if sub_row[1] == 'UR' and sub_row[4] == 'plus':
            ratio_DR_UR_UD_numerator = sub_row[7]
        if sub_row[1] == 'LU' and sub_row[4] == 'plus':
            ratio_IR_LU_UD_C_numerator = sub_row[7]
        if sub_row[1] == 'LU' and sub_row[4] == 'minus':
            ratio_IR_LU_C_UD_numerator = sub_row[7]
        if sub_row[1] == 'RD' and sub_row[4] == 'minus': 
            ratio_IR_RD_C_LR_numerator = sub_row[7]
        if sub_row[1] == 'RD' and sub_row[4] == 'plus':
            ratio_IR_RD_LR_C_numerator = sub_row[7]
            
    for main_row in main_group:
        if main_row[1] == 'LR' and main_row[4] == 'minus':
            ratio_DR_LD_C_LR_C_denominator = main_row[7]
        if main_row[1] == 'LR' and main_row[4] == 'plus':
            ratio_DR_LD_LR_denominator = main_row[7]
        if main_row[1] == 'UD' and main_row[4] == 'minus':
            ratio_DR_UR_C_UD_C_denominator = main_row[7]
        if main_row[1] == 'UD' and main_row[4] == 'plus':
            ratio_DR_UR_UD_denominator = main_row[7]
        if main_row[1] == 'UD' and main_row[4] == 'minus':
            ratio_IR_LU_UD_C_denominator = main_row[7]
        if main_row[1] == 'UD' and main_row[4] == 'plus':
            ratio_IR_LU_C_UD_denominator = main_row[7]
        if main_row[1] == 'LR' and main_row[4] == 'plus':
            ratio_IR_RD_C_LR_denominator = main_row[7]
        if main_row[1] == 'LR' and main_row[4] == 'minus':
            ratio_IR_RD_LR_C_denominator = main_row[7]
                
    # 计算综合性的 ratio
    if (ratio_DR_LD_C_LR_C_denominator + ratio_DR_LD_LR_denominator) > 0:
        ratio_DR_LD = (ratio_DR_LD_C_LR_C_numerator + ratio_DR_LD_LR_numerator) / (ratio_DR_LD_C_LR_C_denominator + ratio_DR_LD_LR_denominator)
        ratio_DR_LD = round(ratio_DR_LD, 6)
    else:
        ratio_DR_LD = 'N/A'
    ratio_DR_LD_LR = f"{ratio_DR_LD_LR_numerator}/{ratio_DR_LD_LR_denominator}"
    ratio_DR_LD_C_LR_C = f"{ratio_DR_LD_C_LR_C_numerator}/{ratio_DR_LD_C_LR_C_denominator}"

    if (ratio_DR_UR_C_UD_C_denominator + ratio_DR_UR_UD_denominator) > 0:
        ratio_DR_UR = (ratio_DR_UR_C_UD_C_numerator + ratio_DR_UR_UD_numerator) / (ratio_DR_UR_C_UD_C_denominator + ratio_DR_UR_UD_denominator)
        ratio_DR_UR = round(ratio_DR_UR, 6)
    else:
       ratio_DR_UR = 'N/A'
    ratio_DR_UR_UD = f"{ratio_DR_UR_UD_numerator}/{ratio_DR_UR_UD_denominator}"
    ratio_DR_UR_C_UD_C = f"{ratio_DR_UR_C_UD_C_numerator}/{ratio_DR_UR_C_UD_C_denominator}"

    if (ratio_IR_LU_UD_C_denominator + ratio_IR_LU_C_UD_denominator) > 0:
        ratio_IR_LU = (ratio_IR_LU_UD_C_numerator + ratio_IR_LU_C_UD_numerator) / (ratio_IR_LU_UD_C_denominator + ratio_IR_LU_C_UD_denominator)
        ratio_IR_LU = round(ratio_IR_LU, 6)
    else:
        ratio_IR_LU = 'N/A'
    ratio_IR_LU_UD_C = f"{ratio_IR_LU_UD_C_numerator}/{ratio_IR_LU_UD_C_denominator}"
    ratio_IR_LU_C_UD = f"{ratio_IR_LU_C_UD_numerator}/{ratio_IR_LU_C_UD_denominator}"


    if (ratio_IR_RD_C_LR_denominator + ratio_IR_RD_LR_C_denominator) > 0:
        ratio_IR_RD = (ratio_IR_RD_C_LR_numerator + ratio_IR_RD_LR_C_numerator) / (ratio_IR_RD_C_LR_denominator + ratio_IR_RD_LR_C_denominator)
        ratio_IR_RD = round(ratio_IR_RD, 6)
    else:
        ratio_IR_RD = 'N/A'
    ratio_IR_RD_LR_C = f"{ratio_IR_RD_LR_C_numerator}/{ratio_IR_RD_LR_C_denominator}"
    ratio_IR_RD_C_LR = f"{ratio_IR_RD_C_LR_numerator}/{ratio_IR_RD_C_LR_denominator}"
    
    return [ratio_DR_LD_LR, ratio_DR_LD_C_LR_C, ratio_DR_LD, ratio_DR_UR_UD, ratio_DR_UR_C_UD_C, ratio_DR_UR, ratio_IR_LU_UD_C, ratio_IR_LU_C_UD, ratio_IR_LU, ratio_IR_RD_LR_C, ratio_IR_RD_C_LR, ratio_IR_RD] 

######################################################################################################################################################################
def extract_number(ratio_str):
    try:
        parts = ratio_str.split('/')
        sub = int(parts[0])  # The first part of the ratio
        main = int(parts[1])  # The second part of the ratio
        return sub, main
    except (ValueError, IndexError, AttributeError):
        # Handle cases where ratio_str is not a string or doesn't follow the 'main/sub' format
        return 0, 0

######################################################################################################################################################################
def safe_ratio_convert(value):
    # Define a helper function to safely handle and convert ratio strings or 'N/A'
    if pd.isna(value) or value == 'N/A':  # Check for NaN or 'N/A'
        return 0, 0
    elif isinstance(value, str) and '/' in value:  # If it's a ratio string, extract both numbers
        return extract_number(value)
    try:
        # If value can be converted directly to float, treat it as the main value with no sub value
        return float(value), 0
    except (ValueError, TypeError):
        return 0, 0

######################################################################################################################################################################
def check_sufficient_main(row, read_number, rev_compl):
    try:
        # Safely get values from row, handling missing columns and converting ratios
        sub_plus, main_plus = safe_ratio_convert(row.get('plus_ratio(s/m)', 0))
        sub_minus, main_minus = safe_ratio_convert(row.get('minus_ratio(s/m)', 0))
        sub_plus_pair, main_plus_pair = safe_ratio_convert(row.get('paired_plus_ratio(s/m)', 0))
        sub_minus_pair, main_minus_pair = safe_ratio_convert(row.get('paired_minus_ratio(s/m)', 0))

        # Adjust the condition based on the value of rev_compl
        if rev_compl:
            # Use 'and' for all conditions if rev_compl is True
            condition_met = (main_plus >= read_number and main_minus >= read_number) and (main_plus_pair >= read_number and main_minus_pair >= read_number)
        else:
            # Use 'or' for any of the conditions if rev_compl is False
            condition_met = (main_plus >= read_number or main_minus >= read_number) and (main_plus_pair >= read_number or main_minus_pair >= read_number)
            
        return 'sufficient' if condition_met else 'insufficient'
        
    except Exception as e:
        loggig.error(f"Error processing row with exception: {e}")
        # sys.stderr.write only the problematic part for debugging
        problematic_part = {k: row.get(k, 'N/A') for k in ['plus_ratio(s/m)', 'minus_ratio(s/m)', 'paired_plus_ratio(s/m)', 'paired_minus_ratio(s/m)']}
        sys.stderr.info(problematic_part)
        return 'insufficient'
        
######################################################################################################################################################################
def check_sufficient_main_unpaired(row, read_number, rev_compl):
    try:
        # Safely get values from row, handling missing columns and converting ratios
        sub_plus, main_plus = safe_ratio_convert(row.get('plus_ratio(s/m)', 0))
        sub_minus, main_minus = safe_ratio_convert(row.get('minus_ratio(s/m)', 0))
        sub_plus_pair, main_plus_pair = safe_ratio_convert(row.get('paired_plus_ratio(s/m)', 0))
        sub_minus_pair, main_minus_pair = safe_ratio_convert(row.get('paired_minus_ratio(s/m)', 0))

        # Adjust the condition based on the value of rev_compl
        if rev_compl:
            # Use 'and' for all conditions if rev_compl is True
            condition_met = (main_plus >= read_number and main_minus >= read_number) and (main_plus_pair >= read_number and main_minus_pair >= read_number)
        else:
            # Use 'or' for any of the conditions if rev_compl is False
            condition_met = (main_plus >= read_number or main_minus >= read_number) and (main_plus_pair >= read_number or main_minus_pair >= read_number)
            
        return 'sufficient' if condition_met else 'insufficient'
        
    except Exception as e:
        print(f"Error processing row with exception: {e}")
        # sys.stderr.write only the problematic part for debugging
        problematic_part = {k: row.get(k, 'N/A') for k in ['plus_ratio(s/m)', 'minus_ratio(s/m)', 'paired_plus_ratio(s/m)', 'paired_minus_ratio(s/m)']}
        print(problematic_part)
        return 'insufficient'

######################################################################################################################################################################
def check_sufficient_sub(row, read_number):

    # Safely get values from row, handling missing columns and converting ratios
    sub_plus, main_plus = safe_ratio_convert(row.get('plus_ratio(s/m)', 0))
    sub_minus, main_minus = safe_ratio_convert(row.get('minus_ratio(s/m)', 0))
    sub_plus_pair, main_plus_pair  = safe_ratio_convert(row.get('paired_plus_ratio(s/m)', 0))
    sub_minus_pair, main_minus_pair = safe_ratio_convert(row.get('paired_minus_ratio(s/m)', 0))

    # Use 'and' for all conditions if rev_compl is True
    condition_met = (sub_plus >= read_number and sub_minus >= read_number) and (sub_plus_pair >= read_number and sub_minus_pair >= read_number)
        
    return condition_met
        
######################################################################################################################################################################
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

######################################################################################################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process TSV data files.")
    parser.add_argument("-sr", "--subconfig_reads", required=True, help="The aggregated result file (repeats with spanning reads) of subconfiguration.")
    parser.add_argument("-mr", "--mainconfig_reads", required=True, help="The aggregated result file (repeats with spanning reads) of mainconfiguration.")
    parser.add_argument("-rf", "--repfile", required=True, help="A 5-column TSV file (*_rep_5CT.tsv) obtained from ROUSFinder' results.")
    parser.add_argument("-rc", "--rep_calibration", required=True, help="A calibrated table (*_rep_calibration_table.txt) obtained from ROUSFinder' results.")
    parser.add_argument("-o", "--output", required=True, help="Output file name for the final combined data.")
    parser.add_argument("-rn", "--read_number", default=1, type=int, help="Spanning read number of repeat supporting mitogenome recombination.")
    parser.add_argument("-cr", "--reverse_complement", required=True, type=str, help="If to consider complementary chains (Yes/No).")
    
    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # Check the number of specific arguments
    required_args = ['subconfig_reads', 'mainconfig_reads', 'repfile', 'rep_calibration', 'output','read_number','reverse_complement']
    if not all(getattr(args, arg, None) for arg in required_args):
        print("ERROR: Incorrect number of external parameters. Please provide the correct parameters.")
        # Add additional information or usage instructions if needed
        exit(1)

    # Assign the file paths directly without opening
    df_sub_path = args.subconfig_reads
    df_main_path = args.mainconfig_reads
    rep_5ct_file_path = args.repfile
    rep_calibr_file_path = args.rep_calibration
        
    read_number = args.read_number
    output_file = args.output
    
    rev_compl = parse_extramaincon_compl_chain(args.reverse_complement)

################################################################################
    #df_sub的数据是根据read_number经过了删减的数据
    df_sub = preprocess_sub_data(df_sub_path, read_number) #预处理subconfiguration数据，根据read_number 标注每一行数据，并写入读取的文件中。
    df_main = preprocess_main_data(df_main_path, read_number) #预处理mainconfiguration数据，根据read_number 标注每一行数据，并写入读取的文件中。

    paired_data_sub, unpaired_data_sub = create_paired_data_sub(df_sub) #将第一列按照"_"拆分后形成新的表，共9列，并重命名每一列的表头，获得配对与不配对的行
    paired_data_main, unpaired_data_main = create_paired_data_main(df_main) #主构型没有删减数据，所以repeat应全部是配对的
    
    # 计算ratio的时候，已经剔除了重组构型中某些repeat中的spanning read数量小于read_number的情况.
    result_paired_ratio = calculate_pair_ratio(paired_data_sub, paired_data_main+unpaired_data_main, rep_5ct_file_path, rep_calibr_file_path)
    # 计算subconfig中仅一个重复单元有spanning reads的repeat的重组比率
    result_unpaired_ratio = calculate_unpair_ratio(unpaired_data_sub, paired_data_main+unpaired_data_main, rep_5ct_file_path, rep_calibr_file_path)

################################################################################
    # In result_paired_ratio and result_unpaired_ratio, add a new column based on conditions, for marking the number of spanning reads in the main configuration
    result_paired_ratio['spanning_read_mcfg'] = result_paired_ratio.apply(check_sufficient_main, args=(read_number, rev_compl), axis=1)
    result_unpaired_ratio['spanning_read_mcfg'] = result_unpaired_ratio.apply(check_sufficient_main_unpaired, args=(read_number, rev_compl), axis=1)

    # Save result_paired_ratio, including the header
    if rev_compl:
        sufficient_rows = result_paired_ratio[result_paired_ratio.apply(check_sufficient_sub, args=(read_number,), axis=1)]
        sufficient_output_path = f"{output_file}/paired_repeats_recomb-supporting_ratio.tsv"
        if not sufficient_rows.empty:
            sufficient_rows.to_csv(sufficient_output_path, sep='\t', index=False)
        else:
            result_paired_ratio.to_csv(sufficient_output_path, sep='\t', index=False)
        
        insufficient_rows = result_paired_ratio[~result_paired_ratio.apply(check_sufficient_sub, args=(read_number,), axis=1)]
        if not insufficient_rows.empty:  # Corrected line
            insufficient_output_path = f"{output_file}/one_chain_without_sufficient_spanning_reads.tsv"
            insufficient_rows.to_csv(insufficient_output_path, sep='\t', index=False)
    else:
        output_file_path = f"{output_file}/paired_repeats_recomb-supporting_ratio.tsv"  # Modify to your desired file save path
        result_paired_ratio.to_csv(output_file_path, sep='\t', index=False)  # index=False means don't save row index, ensuring header is saved
        
    # 保存 result_unpaired_ratio，包含表头
    if not result_unpaired_ratio.empty:
        output_file_path = f"{output_file}/one_repeat_unit_without_spanning_reads.tsv"  # 修改为你要保存文件的路径
        result_unpaired_ratio.to_csv(output_file_path, sep='\t', index=False)  # index=False 表示不保存行索引，确保表头被保存
    
################################################################################
    columns_to_select = ["fragment_id", "start", "end", "direction", "paired_id", "paired_start", "paired_end", "paired_direction", "spanning_read_mcfg"]

    # 从 result_paired_ratio 中选择指定列
    selected_df = result_paired_ratio[columns_to_select]
    # 指定保存文件的文件名
    output_file_name = f"{output_file}/paired_repeats_for_mapping.tsv"
    # 保存选定的数据帧到文件
    selected_df.to_csv(output_file_name, sep='\t', index=False)


