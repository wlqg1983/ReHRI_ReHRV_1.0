#!/usr/bin/env python

import pandas as pd
import os
import argparse

################################################################################################################################################################
def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="Split TSV file by chromosome")
    parser.add_argument('-i', '--input', required=True, help="Path to the input TSV file")
    parser.add_argument('-o', '--output', default='.', help="Path to the output directory (default is the current directory)")
    
    args = parser.parse_args()

    # 读取 TSV 文件，保持所有列为字符串类型
    df = pd.read_csv(args.input, sep='\t', dtype=str).fillna('')

    # 获取所有唯一的染色体名称
    chromosomes = set(df['chromosome']).union(set(df['paired_chromosome']))

    # 逐个染色体进行处理
    for chrom in chromosomes:
        # 选取当前染色体相关的行
        chrom_df = df[(df['chromosome'] == chrom) | (df['paired_chromosome'] == chrom)].copy()

        # 复制必要的列
        new_df = chrom_df[['fragment_id', 'start', 'end', 'direction', 
                           'paired_id', 'paired_start', 'paired_end', 'paired_direction']].copy()

        # 清除不属于当前染色体的部分
        new_df[['fragment_id', 'start', 'end', 'direction']] = new_df[['fragment_id', 'start', 'end', 'direction']].mask(chrom_df['chromosome'] != chrom, '')
        new_df[['paired_id', 'paired_start', 'paired_end', 'paired_direction']] = new_df[['paired_id', 'paired_start', 'paired_end', 'paired_direction']].mask(chrom_df['paired_chromosome'] != chrom, '')

        # 创建 dedup_key（忽略空值，仅比较非空列）
        new_df['dedup_key'] = new_df.apply(lambda row: tuple(row[row != '']), axis=1)

        # 去重
        new_df = new_df.drop_duplicates(subset=['dedup_key']).drop(columns=['dedup_key'])

        # 创建输出目录
        os.makedirs(args.output, exist_ok=True)

        # 保存结果
        output_path = os.path.join(args.output, f'{chrom}.tsv')
        new_df.to_csv(output_path, sep='\t', index=False)

    print(f"Files successfully split and saved to directory: {args.output}")

################################################################################################################################################################
if __name__ == "__main__":
    main()

