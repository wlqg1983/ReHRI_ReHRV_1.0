#!/usr/bin/env python

import pandas as pd
import argparse, os
import warnings

#######################################################################################################################################################################
# Filter out MatplotlibDeprecationWarning
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Filter out FutureWarning for DataFrame concatenation
warnings.filterwarnings("ignore", category=FutureWarning)

#######################################################################################################################################################################
def convert_5ct_to_8ct(five_col_df):
    # 初始化8列格式的DataFrame
    eight_col_df = pd.DataFrame(columns=["fragment_id", "start", "end", "direction", "paired_id", "paired_start", "paired_end", "paired_direction"])

    # 创建一个字典，用于存储fragment_id与其start和end值的映射
    fragment_id_mapping = {}
    for _, row in five_col_df.iterrows():
        direction = 'plus' if row['start'] < row['end'] else 'minus'
        fragment_id_mapping[row['fragment_id']] = (row['start'], row['end'], direction)
    

    # 遍历5列格式的DataFrame
    for _, row in five_col_df.iterrows():
        fragment_id, start, end, repeat_type, paired_ids = row
        direction = 'plus' if row['start'] < row['end'] else 'minus'

        if repeat_type == "direct_repeat":
            for paired_id in paired_ids.split(", "):
                paired_start, paired_end, paired_direction = fragment_id_mapping.get(paired_id, (None, None, None))
                # In direct_repeat, the paired_direction is the same as the paired_id's original direction
                new_row = pd.DataFrame({
                    "fragment_id": [fragment_id],
                    "start": [start],
                    "end": [end],
                    "direction": [direction],
                    "paired_id": [paired_id],
                    "paired_start": [paired_start],
                    "paired_end": [paired_end],
                    "paired_direction": [direction]  # Same direction as paired_id's original direction
                })
                eight_col_df = pd.concat([eight_col_df, new_row], ignore_index=True)

        if repeat_type == "inverted_repeat":
            for paired_id in paired_ids.split(", "):
                paired_start, paired_end, paired_direction = fragment_id_mapping.get(paired_id, (None, None, None))
                # In inverted_repeat, the paired_direction is opposite to the paired_id's original direction
                new_paired_direction = "plus" if direction == "minus" else "minus"
                new_row = pd.DataFrame({
                    "fragment_id": [fragment_id],
                    "start": [start],
                    "end": [end],
                    "direction": [direction],
                    "paired_id": [paired_id],
                    "paired_start": [paired_start],
                    "paired_end": [paired_end],
                    "paired_direction": [new_paired_direction]  # Opposite direction to paired_id's original direction
                })
                eight_col_df = pd.concat([eight_col_df, new_row], ignore_index=True)

        if repeat_type == "two_types":
            direct_ids, inverted_ids = paired_ids.split(" | ")
            for paired_id in direct_ids.split(", "):
                paired_start, paired_end, _ = fragment_id_mapping.get(paired_id, (None, None, None))
                # The direction for direct_ids should be the same as the fragment_id's direction
                new_row = pd.DataFrame({
                    "fragment_id": [fragment_id],
                    "start": [start],
                    "end": [end],
                    "direction": [direction],  # Direction of the fragment_id
                    "paired_id": [paired_id],
                    "paired_start": [paired_start],
                    "paired_end": [paired_end],
                    "paired_direction": [direction]  # Same direction as the fragment_id
                })
                eight_col_df = pd.concat([eight_col_df, new_row], ignore_index=True)
            
            for paired_id in inverted_ids.split(", "):
                paired_start, paired_end, _ = fragment_id_mapping.get(paired_id, (None, None, None))
                # The direction for inverted_ids should be opposite to the fragment_id's direction
                inverted_direction = "minus" if direction == "plus" else "plus"
                new_row = pd.DataFrame({
                    "fragment_id": [fragment_id],
                    "start": [start],
                    "end": [end],
                    "direction": [direction],  # Direction of the fragment_id
                    "paired_id": [paired_id],
                    "paired_start": [paired_start],
                    "paired_end": [paired_end],
                    "paired_direction": [inverted_direction]  # Opposite direction to the fragment_id
                })
                eight_col_df = pd.concat([eight_col_df, new_row], ignore_index=True)
        
        if repeat_type == "-":
            # 对于非重复类型的行，只添加 fragment_id, start, end，并根据 start 和 end 的值设置 direction
            direction = "plus" if start < end else "minus"
            new_row = pd.DataFrame({
                "fragment_id": [fragment_id],
                "start": [start],
                "end": [end],
                "direction": [direction],  # 根据 start 和 end 的值设置 direction
                "paired_id": [None],  # 设置 paired_id 为 None
                "paired_start": [None],
                "paired_end": [None],
                "paired_direction": [None]
            })
            eight_col_df = pd.concat([eight_col_df, new_row], ignore_index=True)

    return eight_col_df

################################################################################
def filter_ctg_rows(eight_col_df):
    # 过滤掉以"ctg"开头的行
    return eight_col_df[~eight_col_df['fragment_id'].str.startswith('ctg')]

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
def main():
    parser = argparse.ArgumentParser(description="Convert 5-column format to 8-column format.")
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file path of 5-column format.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file path for 8-column format.')
    args = parser.parse_args()

    # 读取5列格式的数据
    five_col_df = pd.read_csv(args.input, sep='\t')

    # 转换为8列格式
    eight_col_df = convert_5ct_to_8ct(five_col_df)

    # 过滤掉以"ctg"开头的行
    filtered_df = filter_ctg_rows(eight_col_df)
    
    # 删除重复的行
    #unique_df = remove_duplicates(filtered_df)
    unique_df = drop_duplicates_in_8CT(filtered_df)
    # 检查paired_start和paired_end列是否为空或NaN，并相应地修改paired_id和paired_direction列  
    mask = unique_df[['paired_start', 'paired_end']].isnull().any(axis=1) 
    unique_df.loc[mask, ['paired_id', 'paired_direction']] = ''

    # 构建输出文件路径
    input_basename = os.path.basename(args.input)  # 获取输入文件的基本名称
    output_filename = os.path.splitext(input_basename)[0] + '_map.tsv'  # 构建输出文件名
    output_filepath = os.path.join(args.output, output_filename)  # 构建完整的输出文件路径

    # 输出到文件
    unique_df.to_csv(output_filepath, sep='\t', index=False)

#######################################################################################################################################################################
if __name__ == "__main__":
    main()
