import argparse
import pandas as pd
from Bio import SeqIO
import sys, time, re

################################################################################################################################################################
def clean_chromosome_number(chromosome):
    # 去除所有不可见字符（包括空格、制表符等）
    cleaned_chromosome = re.sub(r'\s+', '', chromosome)
    return cleaned_chromosome

################################################################################
def extract_chromosome_number(chromosome):
    cleaned_chromosome = clean_chromosome_number(chromosome)
    # 检查编号是否符合chr后跟数字的格式
    match = re.search(r'^chr\d+$', cleaned_chromosome)
    if not match:
        raise ValueError("The chromosome number format should be chr1, chr2,..., chr10,..., chr100, etc.")
        exit(1)
    # 提取编号中的数字部分，用于排序
    return int(re.search(r'\d+', cleaned_chromosome).group())

################################################################################
def validate_manual_calibrated_list(manual_calibrated_list, count):
    # 检查列数
    num_columns = manual_calibrated_list.shape[1]

    if num_columns != 6:
        raise ValueError("ERROR: The format of the input repetitive sequence information is incorrect.")
        raise ValueError("The correct header should be: fragment_id, length, start, end, direction, chromosome. (TSV format)\n")
        exit(1)

    # 去重检查染色体编号, # 检查并排序染色体编号
    chromosomes = manual_calibrated_list["chromosome"].unique()
    # 对染色体编号进行清理并检查格式
    cleaned_chromosomes = [clean_chromosome_number(chrom) for chrom in chromosomes]
    # 检查所有编号是否符合格式
    for chromosome in cleaned_chromosomes:
        extract_chromosome_number(chromosome)  # 如果不符合格式，会抛出错误
    
    # 按照清理后的编号进行数字大小排序
    sorted_chromosomes = sorted(cleaned_chromosomes, key=extract_chromosome_number)

    if len(sorted_chromosomes) != count:
        raise ValueError(f"The number of chromosomes you provided is inconsistent with the the number of chromosome numbers in file {manual_calibrated_list}.")
        exit(1)
    else:
        return sorted_chromosomes

################################################################################
def calculate_chromosome_lengths(fasta_file, count):
    seq_lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_lengths.append(len(record.seq))
    if count == 1:
        return [seq_lengths[0]]  # 仅返回一条染色体的长度
    elif count > 1:
        cumulative_lengths = [sum(seq_lengths[:i+1]) for i in range(len(seq_lengths))]
        return cumulative_lengths

################################################################################
def validate_chromosome_headers(fasta_file, chromosomes):
    fasta_headers = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
    if sorted(fasta_headers) != sorted(chromosomes):
        print(f"Warning: The chromosomes in the {fasta_file} will be processed in the given order, treated as chr1, chr2, and so on.")
        time.sleep(5)  # 等待用户终止程序

################################################################################
def adjust_start_end(manual_calibrated_list, cumulative_lengths):
    # Adjust 'start' and 'end' based on cumulative_lengths using apply
    manual_calibrated_list["start"] = manual_calibrated_list.apply(lambda row: row["start"] + cumulative_lengths[int(row["chromosome"][3:]) - 2] if int(row["chromosome"][3:]) > 1 else row["start"], axis=1)
    manual_calibrated_list["end"] = manual_calibrated_list.apply(lambda row: row["end"] + cumulative_lengths[int(row["chromosome"][3:]) - 2] if int(row["chromosome"][3:]) > 1 else row["end"], axis=1)
    return manual_calibrated_list

################################################################################
def check_start_end_max_length(manual_calibrated_list, sorted_chromosomes, cumulative_lengths):
    for i, chromosome in enumerate(sorted_chromosomes):
        # 提取manual_calibrated_list中对应于该染色体的行
        chrom_rows = manual_calibrated_list[manual_calibrated_list["chromosome"] == chromosome]

        if chrom_rows.empty:
            continue

        # 提取start和end的最大值
        max_start = chrom_rows["start"].max()
        max_end = chrom_rows["end"].max()

        # 获取该染色体的长度
        chromosome_length = cumulative_lengths[i] if i == 0 else cumulative_lengths[i] - cumulative_lengths[i - 1]

        # 检查start和end最大值是否超出染色体长度
        if max_start > chromosome_length or max_end > chromosome_length:
            raise ValueError(f"The repetitive sequence position of Chromosome {chromosome} exceeds the length of this chromosome.\nMaximum start: {max_start}, maximum end: {max_end}, chromosome length: {chromosome_length}.\n")

################################################################################
def check_direction(df):
    # 检查 direction 是 plus 时，end > start；direction 是 minus 时，end < start
    for index, row in df.iterrows():
        direction = row['direction']
        start = row['start']
        end = row['end']
        if direction == 'plus' and end <= start:
            raise ValueError(f"Error in row {index + 1}: When direction is 'plus', 'end' should be > 'start'.")
            raise ValueError(row)
            exit(1)
        elif direction == 'minus' and end >= start:
            raise ValueError(f"Error in row {index + 1}: When direction is 'minus', 'end' should be < 'start'.")
            raise ValueError(row)
            exit(1)
            
    # 检查 fragment_id 相同的行，direction 不能全部为 minus
    #grouped = df.groupby('fragment_id')
    #for fragment_id, group in grouped:
    #    directions = group['direction'].tolist()
    #    if all(d == 'minus' for d in directions):
    #        raise ValueError(f"Error: For fragment_id '{fragment_id}', please use 'plus' position information to represent the repeats.")
    #        raise ValueError(group)
    #        exit(1)


######################################################################################################################################################################
################################################################################
def main():
    parser = argparse.ArgumentParser(description="Detect the data format of input files and process chromosome sequences")
    parser.add_argument("-f", "--fasta_file", required=True, help="FASTA file storing chromosome sequences")
    parser.add_argument("-l", "--list_file", required=True, help="manual_calibrated_list file (TSV format)")
    parser.add_argument("-c", "--count", type=int, required=True, help="Chromosome count (count >= 1)")
    
    args = parser.parse_args()
    
    if not isinstance(args.count, int) or args.count < 1:
        parser.error("The number of chromosomes must be a positive integer greater than or equal to 1.")
        exit(1)

    # 读取manual_calibrated_list文件
    manual_calibrated_list = pd.read_csv(args.list_file, sep=r'\s+|,', engine='python')
    
    # 调用 check_direction 函数进行检查
    check_direction(manual_calibrated_list)

    # 计算染色体长度及连加矩阵
    cumulative_lengths = calculate_chromosome_lengths(args.fasta_file, args.count)

    if args.count > 1: 
        # 验证染色体编号
        validate_chromosome_headers(args.fasta_file, manual_calibrated_list['chromosome'].unique())
        # 检查列数是否和染色体的条数相符合
        sorted_chromosomes = validate_manual_calibrated_list(manual_calibrated_list, args.count)
         
        # 使用 pd.Categorical() 将 manual_calibrated_list 中的 chromosome 列设为分类变量，并指定其排序顺序为 sorted_chromosomes。
        manual_calibrated_list["chromosome"] = pd.Categorical(manual_calibrated_list["chromosome"], categories=sorted_chromosomes, ordered=True)
        # 使用 sort_values() 按照 chromosome 列的顺序进行排序，这样 manual_calibrated_list 中的染色体顺序将与 sorted_chromosomes 一致。
        manual_calibrated_list = manual_calibrated_list.sort_values(by="chromosome")
        
        # 检查start和end最大值是否超过染色体长度
        check_start_end_max_length(manual_calibrated_list, sorted_chromosomes, cumulative_lengths)
        
        # 调整start和end
        adjusted_list = adjust_start_end(manual_calibrated_list, cumulative_lengths)
        # 取前5列
        adjusted_list = adjusted_list.iloc[:, :5]
        adjusted_list.to_csv("adjusted_manual_calibrated_list.tsv", sep='\t', index=False)

    if args.count == 1: 
        # 检查列数
        num_columns = manual_calibrated_list.shape[1]
        if num_columns != 5:
            print("The format of the input repetitive sequence information is incorrect.")
            print("The correct header should be: fragment_id, length, start, end, direction. (tsv format)\n")
            exit(1)
        # 取前5列
        adjusted_list = manual_calibrated_list.iloc[:, :5]
        adjusted_list.to_csv("adjusted_manual_calibrated_list.tsv", sep='\t', index=False)       # 保存调整后的列表

if __name__ == "__main__":
    main()
