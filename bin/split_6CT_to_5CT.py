import pandas as pd
import sys
import os

def split_tsv_by_chromosome(input_file):
    # 尝试读取文件并检查分隔符是否为 TSV 格式
    try:
        # 先读取前几行，检查是否为 TSV 格式
        with open(input_file, 'r') as file:
            first_line = file.readline()
            if '\t' not in first_line:
                print("Error: The file does not appear to be in TSV format.")
                sys.exit(1)
        
        # 确认是TSV格式后，读取整个文件
        data = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading the file: {e}")
        sys.exit(1)

    # 检查是否为6列
    if data.shape[1] != 6:
        print("Error: The file does not have exactly 6 columns.")
        sys.exit(1)
        
    # 检查文件是否只有一行（表头）
    if len(data) == 0:
        print("Error: The repeat sequence information file may be empty!")
        sys.exit(1)
        
    # 获取原始文件的表头
    original_headers = data.columns[:5]

    # 按第六列 'chromosome' 进行分组
    grouped = data.groupby('chromosome')
    
    # 获取当前目录作为输出目录
    output_dir = os.getcwd()

    # 遍历每个组，生成新表并保存为文件
    for chromosome, group in grouped:
        # 只保留前五列
        new_table = group.iloc[:, :5]
        new_table.columns = original_headers  # 保留原始表头

        # 保存新表为以 'chromosome' 命名的文件，格式为 TSV
        filename = os.path.join(output_dir, f'{chromosome}_5CT.tsv')
        new_table.to_csv(filename, sep='\t', index=False)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python split_tsv_by_chromosome.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    split_tsv_by_chromosome(input_file)
