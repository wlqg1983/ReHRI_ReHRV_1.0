#!/usr/bin/env python

import argparse, os, sys
import configparser
import subprocess
from Bio import SeqIO
import shutil, glob
import time, re, io
import logging, signal
from itertools import accumulate  
import pandas as pd 
from colorama import Fore, Style, init  
import numpy as np
from multiprocessing import Pool, cpu_count

#######################################################################################################################################################################
def read_ini_file(ini_file_path):
    config = configparser.ConfigParser()
    config.optionxform = str
    ini_content = []
    with open(ini_file_path, 'r') as file:
        for line in file:
            # 跳过日志行，只保留INI配置行
            if line.strip().startswith('[') or '=' in line:
                ini_content.append(line)
    config.read_string(''.join(ini_content))
    return config

#######################################################################################################################################################################
def extract_ini_from_log(log_file, project_id):
    ini_content = []
    is_ini_section = False

    with open(log_file, 'r') as file:
        for line in file:
            if re.match(r'\[general\]', line):  # Detect the start of the INI section
                is_ini_section = True
            if is_ini_section:
                ini_content.append(line)
                if 'project_id' in line:
                    project_id = line.split('=')[1].strip()
                if 'spanning_read_number' in line:  # Detect the end of the INI section
                    break

    ini_content_str = ''.join(ini_content)
    
    if project_id:
        output_filename = f'old_log_{project_id}.ini'
        with open(output_filename, 'w') as output_file:
            output_file.write(ini_content_str)
            
#######################################################################################################################################################################
def compare_ini_contents(old_ini_file, new_ini_file, keys_to_compare):
    old_config = configparser.ConfigParser()
    old_config.optionxform = str
    new_config = configparser.ConfigParser()
    new_config.optionxform = str

    old_config.read(old_ini_file)
    new_config.read(new_ini_file)

    has_differences = False  # 默认没有差异
    
    for section in keys_to_compare:
        if not old_config.has_section(section):
            has_differences = True
            print(f"Section '{section}' not found in old INI.")
            continue
        
        for key in keys_to_compare[section]:
            old_value = old_config.get(section, key, fallback="Not found")
            new_value = new_config.get(section, key, fallback="Not found")

            if old_value != new_value:
                has_differences = True
                print("The modified parameters are as follows:")
                print(f"Section: {section}, Key: {key}, Old Value: {old_value}, New Value: {new_value}")

    return has_differences

#######################################################################################################################################################################
def extract_trimmed_length_from_filename(filename):
    # 修改正则表达式以匹配新的文件名格式
    match = re.search(r'_(\d+)\.fasta$', filename)
    if match:
        length = int(match.group(1))
        if length == 0:
            # 对 0 进行特殊处理
            return False  # 或者返回其他你希望的值
        return length
    else:
        return False

#######################################################################################################################################################################
def run_command(command, output_file=None, append=False, add_header=False, header=None):
    try:
        # 记录即将执行的命令
        logging.info(f"Executing command: {' '.join(command)}")

        # 修改subprocess.run调用，捕获stdout和stderr
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)

        # 记录标准输出
        if result.stdout:
            logging.info(f"Command stdout: {result.stdout}")

        # 记录标准错误（如果有的话）
        if result.stderr:
            # 将特定警告信息记录为WARNING，其余信息保持为ERROR
            if "No reads span the repeat unit" in result.stderr:
                logging.warning(f"Command stderr: {result.stderr}")
            elif "Disregard two overlapping repeats marked with" in result.stderr:
                logging.warning(f"Command stderr: {result.stderr}")
            else:
                logging.error(f"Command stderr: {result.stderr}")

        # 处理输出文件
        if output_file:
            mode = 'a' if append else 'w'
            with open(output_file, mode) as f:
                if add_header and header and (not os.path.isfile(output_file) or os.path.getsize(output_file) == 0):
                    f.write(header + '\n')
                f.write(result.stdout)
        return result.stdout
    except subprocess.CalledProcessError as e:
        # 使用logging记录异常
        logging.error(f"Error executing command: {' '.join(command)}")
        logging.error(f"Error message: {e.stderr}")
        sys.exit(1)
    except Exception as e:
        # 捕获其他异常
        logging.error(f"Unexpected error: {e}")
        sys.exit(1)
        
#######################################################################################################################################################################
def parse_redundant_intermediate_results(value):
    if isinstance(value, str):
        normalized_value = value.strip().upper()
        if normalized_value == "D":
            return True
        elif normalized_value == "K":
            return False
        else:
            # 记录错误日志
            logging.error(f"Invalid input for redundant_intermediate_results: {value}. Must be 'D' or 'K'.")
            raise ValueError("Valid input: must be 'D' or 'K'.")
    else:
        # 记录错误日志
        logging.error(f"Type error in redundant_intermediate_results: Input type is {type(value)}, but must be a string.")
        raise TypeError("Input must be a string.")

#######################################################################################################################################################################
def validate_positive_integer(value, param_name):
    # Validate that the parameters are greater than zero and integers
    if not isinstance(value, int) or value <= 0:
        error_message = f"{param_name} must be a positive integer (>=1). Given value: {value}"
        # Log the error message
        logging.error(error_message)
        # Optionally, you can still raise the error to stop execution or signal that an error condition has occurred
        raise ValueError(error_message)

#######################################################################################################################################################################
def is_valid_project_id(project_id):
    """ 检查 project_id 是否有效： - 只能包含字母、数字、下划线  - 必须以字母开头 """
    # 定义正则表达式
    pattern = r'^[A-Za-z][A-Za-z0-9_]*$'
    # 使用 re.match 检查 project_id 是否匹配正则表达式
    if re.match(pattern, project_id):
        return True
    else:
        return False

#######################################################################################################################################################################
def parse_filter_again(value):
    if isinstance(value, str):
        normalized_value = value.strip().upper()
        if normalized_value == "YES" or normalized_value == "Y":
            return True
        elif normalized_value == "NO" or normalized_value == "N":
            return False

#######################################################################################################################################################################
def setup_logging(mode, project_id, log_file_name=None, include_file_handler=True):
    # 检查mode参数是否正确
    if mode not in ['w', 'a', 'w+']:
        raise ValueError("mode参数必须是'w'、'a'或'w+'之一")

    handlers = []
    if include_file_handler:
        if log_file_name is None:
            log_file_name = f"custom_log_{project_id}.txt"
        try:
            handlers.append(logging.FileHandler(log_file_name, mode=mode))
            print(f"Log file will be created: {log_file_name}")  # 添加调试信息
        except Exception as e:
            print(f"Failed to create log file: {e}")
    handlers.append(logging.StreamHandler(sys.stdout))

    # 检查是否已经配置过basicConfig
    if not logging.root.handlers:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=handlers)
    return handlers

#######################################################################################################################################################################
def setup_refilter_logging(mode, refilter_id, log_file_name=None, include_file_handler=True):
    handlers = [logging.StreamHandler(sys.stdout)]
    if include_file_handler:
        if log_file_name is None:
            log_file_name = f"refilter_log_{refilter_id}.txt"
        handlers.append(logging.FileHandler(log_file_name, mode=mode))
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', handlers=handlers)
    return handlers  # 返回处理器列表以便后续操作
    
#######################################################################################################################################################################
def write_ini_to_log(ini_file_path):
    config = configparser.ConfigParser()
    config.optionxform = str
    config.read(ini_file_path)

    ini_content = ""
    for section in config.sections():
        ini_content += f"[{section}]\n"
        for key, value in config.items(section):
            ini_content += f"{key} = {value}\n"
        ini_content += "\n"

    # 记录日志时包含 ini_file_path
    logging.info(f"INI file path: {ini_file_path}\nINI file content:\n{ini_content}")

#######################################################################################################################################################################
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

#######################################################################################################################################################################
def check_fasta_sequence(fasta_file_path):  
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
        return len(sequences) 
    elif len(sequences) == 0:  
        return 0 
    else:  
        return 1 

#######################################################################################################################################################################
def process_fastq_chunk(chunk):
    """
    处理FASTQ文件的一个块，并返回转换后的FASTA内容
    :param chunk: FASTQ文件的一个块（包含多条记录）
    :return: 转换后的FASTA内容
    """
    fasta_lines = []
    for i in range(0, len(chunk), 4):
        header = chunk[i].strip()
        sequence = chunk[i + 1].strip()
        plus_line = chunk[i + 2].strip()
        quality = chunk[i + 3].strip()

        # 检查格式
        if not header.startswith('@'):
            raise ValueError(f"FASTQ文件格式错误：第{i + 1}行不是以'@'开头")
        if not plus_line.startswith('+'):
            raise ValueError(f"FASTQ文件格式错误：第{i + 3}行不是以'+'开头")
        if len(sequence) != len(quality):
            raise ValueError(f"FASTQ文件格式错误：第{i + 2}行的序列长度与第{i + 4}行的质量信息长度不一致")

        # 转换为FASTA格式
        fasta_lines.append(f">{header[1:]}\n{sequence}\n")

    return "\n".join(fasta_lines)

###############################################################################
def fastq_to_fasta_multiprocess(fastq_file, fasta_file, num_processes=None):
    """
    使用多进程将FASTQ文件转换为FASTA文件
    :param fastq_file: 输入的FASTQ文件路径
    :param fasta_file: 输出的FASTA文件路径
    :param num_processes: 使用的进程数，默认为CPU核心数
    """
    if num_processes is None:
        num_processes = cpu_count()  # 默认使用所有CPU核心

    # 读取FASTQ文件并分块
    with open(fastq_file, 'r') as fq:
        lines = fq.readlines()

    # 检查文件是否完整
    if len(lines) % 4 != 0:
        raise ValueError("FASTQ文件格式错误：文件不完整，记录数不是4的倍数")

    # 每个进程处理的记录数
    chunk_size = len(lines) // (4 * num_processes)
    if chunk_size == 0:
        chunk_size = 1

    # 创建进程池
    with Pool(processes=num_processes) as pool:
        # 将FASTQ文件分块并分配给进程池
        chunks = [lines[i:i + 4 * chunk_size] for i in range(0, len(lines), 4 * chunk_size)]
        results = pool.map(process_fastq_chunk, chunks)

    # 将所有结果合并并写入FASTA文件
    with open(fasta_file, 'w') as fa:
        fa.write("\n".join(results))

    print(f"Conversion completed, FASTA file has been saved to {fasta_file}")
    return fasta_file  # 返回生成的FASTA文件路径

#######################################################################################################################################################################
def concatenate_fasta(input_fasta, output_fasta):
    records = []
    headers = []
    lengths = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        records.append(record)
        headers.append(record.description)
        lengths.append(len(record.seq))

    concatenated_records = [records[0]]
    for i in range(1, len(records)):
        concatenated_records[0].seq += records[i].seq
        concatenated_records[0].description += " " + records[i].description

    with open(output_fasta, "w") as output_handle:
        SeqIO.write(concatenated_records, output_handle, "fasta")

    return output_fasta, headers, lengths

################################################################################
def parse_repeat_lengths(repeat_length_str, sequence_lengths):
    # Parse the repeat_length configuration item and return a list of integers. 
    # sequence_lengths is the maximum length of the generated sequence, used to limit the upper range for open-ended lengths (like 30:). 
    repeat_length_str = re.sub(r'\s+', '', repeat_length_str)
    lengths = []
    min_value = 5
    parts = repeat_length_str.split(',')

    if is_valid_repeat_length_str(repeat_length_str):
        for part in parts:
            part = part.strip()
            if ':' in part:
                start, end = part.split(':')
                
                try:
                    start = int(start) if start else min_value
                    validate_positive_integer(start, 'start')
                except ValueError as e:
                    logging.error(f"Validation error for start: {e}")
                    sys.exit(1)
                
                try:
                    end = int(end) if end else sequence_lengths
                    validate_positive_integer(end, 'end')
                except ValueError as e:
                    logging.error(f"Validation error for end: {e}")
                    sys.exit(1)
            
                if start <= end:
                    lengths.extend(range(start, min(end + 1, sequence_lengths + 1)))
                else:
                    start, end = end, start
                    lengths.extend(range(start, min(end + 1, sequence_lengths + 1)))
            else:
                try:
                    length = int(part)
                    validate_positive_integer(length, 'length')
                except ValueError as e:
                    logging.error(f"Validation error for length: {e}")
                    sys.exit(1)

                if length <= sequence_lengths:
                    lengths.append(length)
        return lengths
    else:
        logging.error(f"The 'repeat_length' provided in the [ROUSFinder] section was invalid!\n")
        sys.exit(1)

################################################################################
def is_valid_repeat_length_str(repeat_length_str):
    # Validates repeat length string format: "number", "number:number", "number:", or ":number"
    pattern = r'^(\d+:\d*|\d*:\d+|\d+)(,\s*(\d+:\d*|\d*:\d+|\d+))*$'
    return bool(re.match(pattern, repeat_length_str))

################################################################################
# 重新归档产生的结果时，处理伪基因组内 start 和 end 的真实位点 
def find_chromosome(start, end, cumulative_lengths):
    # 将start和end转换为整数
    try:
        start = int(start)
        end = int(end)
    except ValueError:
        raise ValueError("The 'start' and 'end' must be integers or values that can be converted to integers.")

    for i, length in enumerate(cumulative_lengths):
        # 确保比较时数据类型一致
        if start <= length and end <= length:
            if i == 0:
                return f'chr{i+1}', start, end
            else:
                # 计算调整后的起始和结束位置
                adjusted_start = start - cumulative_lengths[i-1]
                adjusted_end = end - cumulative_lengths[i-1]
                return f'chr{i+1}', adjusted_start, adjusted_end

################################################################################
def is_file_empty_or_invisible_chars(file_path):  
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:  
        content = file.read()  
        # 检查文件是否为空  
        if not content:  
            return True  
        # 检查是否所有字符都是不可见的  
        for char in content:  
            if ord(char) >= 32:  # 假设ASCII码32以下的字符为不可见字符  
                return False  
        return True  

################################################################################
def get_ordinal_suffix(number):  
    # Helper function to get the ordinal suffix for a given number.  

    if number % 10 == 1 and number % 100 != 11:  
        suffix = 'st'  
    elif number % 10 == 2 and number % 100 != 12:  
        suffix = 'nd'  
    elif number % 10 == 3 and number % 100 != 13:  
        suffix = 'rd'  
    else:  
        suffix = 'th'  
    return suffix 

################################################################################
#对输入的fasta的header进行处理，以免意外报错
def replace_spaces_in_headers(fasta_file):
    """
    Replace spaces with underscores in the headers of a FASTA file.
    
    Parameters:
    fasta_file (str): Path to the input FASTA file.
    """
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    # 替换以 '>' 开头的行中的空格为下划线
    modified_lines = [line.replace(' ', '_') if line.startswith('>') else line for line in lines]

    with open(fasta_file, 'w') as f:
        f.writelines(modified_lines)

def modify_headers_in_place(fasta_file):
    """
    Modify the headers of each sequence in the FASTA file by removing all 'chr{i}_' prefixes.
    The 'i' in 'chr{i}_' is the sequence's order in the file.
    
    Parameters:
    fasta_file (str): Path to the input FASTA file.
    """
    records = []
    sequence_index = 1  # 用于记录序列的顺序
    
    # 读取并修改所有序列的header
    for record in SeqIO.parse(fasta_file, "fasta"):
        # 检查header是否包含'chr{i}_'格式的前缀，并删除它
        while record.id.startswith('chr') and '_' in record.id:
            prefix = 'chr' + str(sequence_index) + '_'
            if record.id.startswith(prefix):
                new_id = record.id[len(prefix):]
            else:
                new_id = record.id.split('_', 1)[1]  # 只取第一个'_'后面的部分
            record.id = new_id
            sequence_index += 1  # 更新序列的顺序
        
        # 更新record的description
        record.description = ''
        
        records.append(record)
    
    # 将修改后的序列写回原文件
    with open(fasta_file, 'w') as out_handle:
        SeqIO.write(records, out_handle, "fasta")
        
#######################################################################################################################################################################
#######################################################################################################################################################################
def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="MiRI: A tool to check spanning reads for supporting subconfig of your organelle genome.")
    parser.add_argument("-c", dest="config", help="Path to external configuration file.", required=True)
    parser.add_argument("-redo", help="Delete all previous results and start calculation anew.", action="store_true")
    parser.add_argument("-resume", action="store_true", help="Resume from a previous project.")
    parser.add_argument("-v", "--version", action="version", version="MiRI 1.0", help="Show the version number and exit.")

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # Load external configuration file
    config = read_ini_file(args.config)

    project_id = config['general']['project_id']
    
    # 检验 project_id 是否有效。无效时，退出程序。
    if is_valid_project_id(project_id):
        print(f"The 'project_id' {project_id} is valid.\n")
    else:
        print("Invalid 'project_id'. The 'project_id' must consist of letters, numbers, and underscores, and must start with a letter.")
        sys.exit(1)
        
    mode = config['general'].get('mode', 'A').upper().strip() or 'A'
    refilter_mode = config['refilter_params'].get('refilter_mode', 'no').upper().strip() or 'NO'
    
    if refilter_mode not in ['YES','Y','NO','N']:
        print("Invalid paramater 'refilter_mode' in the [refilter_params] section. It must be one of 'yes', 'y', 'no' and 'n'.")
        sys.exit(1)

    refilter_mode = parse_filter_again(refilter_mode)

    # 检查mode值是否有效
    if mode not in ['A', 'R', 'C', 'N']:
        print(f"Invalid 'mode' value '{mode}' in the [general]. It must be one of 'N', 'A', 'R' or 'C'.")
        sys.exit(1)

    ####### 对日志文件进行管理
    # 删除 log.txt，防止后面 判断与前一次运行结果冲突的逻辑 出现问题
    # log.txt 在进行二次数据结果过滤的时候要用到，所以
    if not args.resume: 
        os.remove(f"custom_log_{project_id}.txt") if os.path.isfile(f"custom_log_{project_id}.txt") and not refilter_mode else None

    # 配置日志系统，确保在 resume 模式下也能记录日志
    log_file_name = f"custom_log_{project_id}.txt"
    mode_for_log = 'a' if args.resume else 'w'  # 如果是 resume 模式，则以追加模式打开日志文件
    handlers = setup_logging(mode=mode_for_log, project_id=project_id, log_file_name=log_file_name, include_file_handler=True)

    def signal_handler(signal, frame):
        logging.info('Program has been terminated by user.\n')
        for handler in handlers:
            handler.flush()
            handler.close()
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)
    
    # write *.ini content to "log.txt".
    write_ini_to_log(args.config)
    
################################################################################
    ########################################################
    ####主程序开始

    if (mode == 'N' or mode == 'NO') and not refilter_mode:
        logging.info("No work to be executed. Please reset 'mode' values!\n")
        sys.exit(0)

    if mode in ['A', 'R', 'C'] and not refilter_mode:
        # 检查 args.redo 和 args.resume 是否同时被提供
        if args.redo and args.resume:
            raise ValueError("Parameters '-redo' and '-resume' cannot be provided at the same time.")
            sys.exit(1)

        # General parameters
        genome_type = config['general'].get('genome_type', 'N').upper().strip() or 'N'     # genome_type 当为空值的时候，会被赋值 N

        if genome_type not in ['C', 'L', 'N']:
            logging.error("Invalid paramater 'genome_type' in [general]. It must be 'C' or 'L' when one sequence is present, or left unconfigured when two or more sequences are present.")
            sys.exit(1)
        
        # 将文件夹内的f"{project_id}"文件及文件夹转移至当前工作目录下
        if args.resume and os.path.isdir(f"{project_id}"):
            # Use glob to list all files and directories in the project_id directory  
            for file_path_resume in glob.glob(f"{project_id}/*"):  
                # Extract the filename from the full path  
                file_name_resume = os.path.basename(file_path_resume)  
                # Construct the new path in the current working directory  
                new_path_resume = os.path.join(os.getcwd(), file_name_resume)  

                if os.path.exists(new_path_resume):
                    # 如果是文件，使用 os.remove 删除
                    if os.path.isfile(new_path_resume):
                        os.remove(new_path_resume)
                    # 如果是文件夹，使用 shutil.rmtree 删除
                    elif os.path.isdir(new_path_resume):
                        shutil.rmtree(new_path_resume)
                # Move the file or directory  
                shutil.move(file_path_resume, new_path_resume) 
                
                
        #### 断点续接，需要日志文件和处理过程记录文件
        ## 检查处理过程记录文件是否存在
        if args.resume and not os.path.isfile(f"record_for_resume_{project_id}.txt"):      # 用于断点续接的文件不存在，终止程序
            logging.error(f"The '-resume' can't be used now. Please use '-redo' instead.")
            print(f"You can check if the PROJECT labeled by '{project_id}' has already been completed. If so, reset a new 'project_id'.\n")
            sys.exit(1)
        ## 检查日志文件是否存在
        if args.resume and not os.path.isfile(f"custom_log_{project_id}.txt"): 
            logging.error(f"Previous '.ini' file {args.config} has not been found.")
            logging.info("You can use '-redo' instead of '-resume', or set a new 'project_id'.\n")
            sys.exit(1)

        # 在断点继续计算时，先检查日志文件有没有发生变化
        if args.resume:
            # Define the keys to compare
            keys_to_compare = {
                'general': ['project_id', 'mode', 'inputfasta', 'genome_type', 'complementary_chain', 'redundant_intermediate_results'],
                'ROUSFinder': ['repeat_length', 'reward', 'penalty'],
                'manually_calibrated_repeat_info': ['calibrated_repeat_file'],
                'mainconfiguration': ['flanked_sequence_length'],
                'subconfiguration': ['flanked_sequence_length'],
                'sequencing_depth': ['alignment_software', 'evalue', 'threads', 'NGS_single_end', 'NGS_pair_ends', 'TGS', 'TGS_type', 'filter_reads'],
                'check_spanning_reads': ['spanning_read_number'],
                }

            # Read old INI content from log.txt
            if os.path.isfile(f"custom_log_{project_id}.txt"):
                extract_ini_from_log(f"custom_log_{project_id}.txt", project_id)
                # 前后两个配置的内容进行比较
                has_difference = compare_ini_contents(f"old_log_{project_id}.ini", args.config, keys_to_compare)
                if has_difference:
                    logging.error(f"Previous '.ini' file {args.config} has been changed.")
                    print("You can use '-redo' instead of '-resume'.\n")
                    sys.exit(1)     # 配置信息出现变化，则终止程序
            
            if os.path.isfile(f"old_log_{project_id}.ini"):      # 删除用于检查配置文件有没有发生变化的临时配置文件
                os.remove(f"old_log_{project_id}.ini")
            
        logging.info("Current work is to SEARCH for paired repeats that support genomic recombination based on ROUSFinder results.") if mode == 'A' else None
        logging.info("Current work is to SEARCH for paired repeats that support genomic recombination based on user provided repeat info.") if mode == 'C' else None
        time.sleep(3)

################################################################################
        inputfasta = config['general']['inputfasta']
        output_dir_prefix = config['general'].get('output_directory_prefix', 'final_repeat-spanning_results').strip() or 'final_repeat-spanning_results'
        redundant_intermediate_results = config['general'].get('redundant_intermediate_results', 'D').upper().strip() or 'D' 
        complementary_chain = config['general'].get('complementary_chain', 'Yes').upper().strip() or 'YES'

        # Check mandatory parameters in the [general] section  resume
        try:
            project_id = config['general']['project_id']
            inputfasta = config['general']['inputfasta']
        except KeyError as e:
            logging.exception(f"Missing mandatory configuration parameter(s): {e}")
            exit(1)
                 
        if complementary_chain not in ['YES','Y','NO','N']:
            logging.error(f"Error: Invalid 'complementary_chain' '{complementary_chain}' parameter in the [general] section. It must be one of 'Yes', 'Y', 'No' or 'N'.")
            sys.exit(1)
            
        if redundant_intermediate_results not in ['D','K']:
            logging.error(f"Error: Invalid 'redundant_intermediate_results' '{redundant_intermediate_results}' parameter in the [general] section. It must be one of 'D' or 'K'.")
            sys.exit(1)

        # 删除pseudo-genome，防止后面 判断与前一次运行结果冲突的逻辑出现问题
        os.remove(f"cat_inputfasta_{project_id}.fasta") if os.path.isfile(f"cat_inputfasta_{project_id}.fasta") else None        
        
        # 输入的fasta序列，其序列的header中可能存在空格，为防止出现意外报错，事先用"_"替代空格
        replace_spaces_in_headers(inputfasta)    # Replace spaces with underscores in the headers of a FASTA file.
        logging.info(f"ATTENTION: The spaces in the header of {inputfasta} are replaced with '_'.")
        time.sleep(2)

        # Calculate the lengths of sequences in the FASTA file
        count = check_fasta_sequence(inputfasta)    # 计数fasta序列的条数
        # 检查 inputfasta 是否为空
        if not os.path.isfile(inputfasta):
            logging.error("When 'mode' in the [general] section is one of ['A', 'C', 'R'], 'inputfasta' in section [general] cannot be empty.")
            sys.exit(1)
            
        if count == 0:
            logging.error(f"No sequences in the {inputfasta} file.")  
            sys.exit(1)                
        if count == 1:
            logging.info(f"There are {count} sequences in the input fasta file {inputfasta}.")
            if genome_type not in ['C', 'L']:      # 设置基因组的类型
                logging.error("Invalid paramater 'genome_type' in [general]. It must be 'C' or 'L' when one sequence is present.")
                sys.exit(1)
                
        if count > 1:
            cat_inputfasta_path, headers, lengths = concatenate_fasta(inputfasta, f"cat_inputfasta_{project_id}.fasta")
            logging.info(f"There are {count} sequences in the input fasta file {inputfasta}. They are concatenated into a pseudo-genome for the further process.")
            # 使用itertools.accumulate进行  各条染色体长度的累积求和  
            cumulative_length = list(accumulate(lengths)) 
            #logging.info("The chromosomes you provided will be labeled as chr1, chr2, (case sensitive), and so on, based on their order in the FASTA file you provide.")
            time.sleep(5)
            if genome_type == 'N':
                genome_type = 'L'
            else:
                logging.error("Invalid paramater 'genome_type' in [general]. It should left unconfigured when two or more sequences are present.")
                sys.exit(1)
            
        sequence_lengths = []
        if  count == 1:
            for record in SeqIO.parse(inputfasta, "fasta"):
                sequence_lengths = len(record.seq)    # 一条序列的长度
            logging.info(f"The length of the provided genome is: {sequence_lengths} bp.")

        if  count > 1:
            for record in SeqIO.parse(cat_inputfasta_path, "fasta"):
                sequence_lengths =+ len(record.seq)    # 多条序列的长度和
            logging.info(f"The total length of the {count} chromsomes is: {sequence_lengths} bp.")
            
            os.remove(cat_inputfasta_path) if os.path.isfile(cat_inputfasta_path) else None        # 删除pseudo-genome，防止后面 判断与前一次运行结果冲突的逻辑出现问题
            
            # 为输入的fasta序列更改名字
            dir_path = os.path.dirname(os.path.realpath(__file__))
            subprocess.run(["python", os.path.join(dir_path, "modify_fasta_header.py"), inputfasta], check=True) 

################################################################################
        # [ROUSFinder] Section
        rous_output_file_prefix = config['ROUSFinder'].get('output_file_prefix')
        repeat_length_str = config['ROUSFinder'].get('repeat_length', '50:') or '50:'
        rous_reward = config['ROUSFinder'].getint('reward', 1)
        rous_penalty = config['ROUSFinder'].getint('penalty', 20)

        # Check if rous_output_file_prefix is empty, if so, use the prefix of inputfasta
        if not rous_output_file_prefix:
            # 使用 os.path.basename 去掉路径，只保留文件名  
            basename_input = os.path.basename(inputfasta)  
            rous_output_file_prefix =  '.'.join(basename_input.split('.')[:-1]) 

        repeat_lengths = parse_repeat_lengths(repeat_length_str, sequence_lengths)    # 返回一个包含各种重复序列长度的列表
        rous_repeat_minlength = min(repeat_lengths)    # 后面在查找重复序列的时候，ROUSFinder 仅能接受一个数值，作为最小重复序列的长度，此处提取 repeat_lengths 的最小值进行搜索

################################################################################
        # [manually_calibrated_repeat_info] Section
        calibrate_flexibility = config['manually_calibrated_repeat_info'].getint('flexibility', 0)    # 预留接口
        manually_calibrate = config['manually_calibrated_repeat_info'].get('calibrated_repeat_file', '').strip()     # 读取手工修改的重复序列信息文件
        
        #检查重复序列单元的长度不小于5bp
        if rous_repeat_minlength < 5:
            logging.error("The repeat length must be NOT less than 5 bp.")
            sys.exit(1)
                
        # Check if manually_calibrate is empty when mode is 'C'
        if mode == 'C' and not os.path.isfile(manually_calibrate):
            logging.error("The 'calibrated_repeat_file' parameter in the '[manually_calibrated_repeat_info]' section must be provided when 'mode == C' in the [general] section.")
            sys.exit(1)
            
        if mode == 'C' and not os.path.isfile(inputfasta):
            logging.error(f"The 'inputfasta' parameter in the '[general]' section and the 'calibrated_repeat_file' parameter in the '[manually_calibrated_repeat_info]' must be provided at the same time when 'mode == C' in the [general] section.")
            sys.exit(1)

        try:
            validate_positive_integer(rous_repeat_minlength, 'rous_repeat_minlength')
        except ValueError as e:
            logging.error(f"Validation error for rous_repeat_minlength: {e}")
            sys.exit(1)
        try:
            validate_positive_integer(rous_reward, 'rous_reward')
        except ValueError as e:
            logging.error(f"Validation error for rous_reward: {e}")
            sys.exit(1)
        try:
            validate_positive_integer(rous_penalty, 'rous_penalty')
        except ValueError as e:
            logging.error(f"Validation error for rous_penalty: {e}")
            sys.exit(1)

################################################################################
        # [sequencing_depth] Section  
        seqdepth_alignment_software = config['sequencing_depth'].get('alignment_software', 'minimap2').lower().strip() or 'minimap2'
        blast_evalue = config['sequencing_depth'].get('evalue', '1e-5').strip() or '1e-5'    
        seqdepth_single = config['sequencing_depth'].get('NGS_single_end', '').strip()
        seqdepth_pair = config['sequencing_depth'].get('NGS_pair_ends', '').strip()
        seqdepth_third = config['sequencing_depth'].get('TGS', '').strip()
        seqdepth_threads = config['sequencing_depth'].get('threads', '').strip()  
        seqdepth_type = config['sequencing_depth'].get('TGS_type', '').lower().strip()
        filter_reads = config['sequencing_depth'].get('filter_reads', 'Y').upper().strip() or 'Y'

        if seqdepth_alignment_software not in ['minimap2', 'bwa', 'blast']:
            logging.error("The mapping software must be selected 'blast', 'minimap2' or 'bwa'.")
            exit(1)
             
        if filter_reads not in ['YES','Y','NO','N']:
            logging.error(f"Error: Invalid 'filter_reads' '{filter_reads}' parameter in the [sequencing_depth] section. It must be one of 'Yes', 'Y', 'No' or 'N'.")
            sys.exit(1)

        # Initialize counter for provided data types
        provided_data_count = 0
        # Check sequencing read data and count provided types
        if seqdepth_single:
            provided_data_count += 1
            provided_data_type = 'NGS_single_end'
            provided_data = seqdepth_single
            if os.path.isfile(provided_data):
                file_format = check_file_format_efficient(provided_data)      # 检查输入文件的数据格式
                if file_format == 'FASTA':
                    logging.info(f"The infile {provided_data} is in FASTA format.")
                elif file_format == 'FASTQ':
                    logging.info(f"The infile {provided_data} is in FASTQ format.")
                else:
                        logging.error(f"The infile {provided_data} format is unknown, empty or invalid.")
                        sys.exit(1)
            else:
                logging.error(f"The infile {provided_data} does not exist.")
                sys.exit(1)

        if seqdepth_pair:
            provided_data_count += 1
            provided_data_type = 'NGS_pair_ends'
            provided_data_paths = seqdepth_pair.split()    # 拆分为两个文件路径 
            for provided_data in provided_data_paths:
                if os.path.isfile(provided_data):
                    file_format = check_file_format_efficient(provided_data)
                    if file_format == 'FASTA':
                        logging.info(f"The infile {provided_data} is in FASTA format.")
                    elif file_format == 'FASTQ':
                        logging.info(f"The infile {provided_data} is in FASTQ format.")
                    else:
                        logging.error(f"The infile {provided_data} format is unknown, empty or invalid.")
                        sys.exit(1)
                else:
                    logging.error(f"The infile {provided_data} does not exist.")
                    sys.exit(1)

        if seqdepth_third:
            provided_data_count += 1
            provided_data_type = 'TGS'
            provided_data = seqdepth_third
            # Check if third type is either 'ont' or 'pacbio'
            if seqdepth_type not in ['ont', 'pacbio']:
                logging.error("Third-generation sequencing data type is Required, which must be either 'ont' or 'pacbio'.")
                exit(1)
            if os.path.isfile(provided_data):
                file_format = check_file_format_efficient(provided_data)
                if file_format == 'FASTA':
                    logging.info(f"The infile {provided_data} is in FASTA format.")
                elif file_format == 'FASTQ':
                    logging.info(f"The infile {provided_data} is in FASTQ format.")
                elif file_format == 'Unknown': 
                        logging.error(f"The infile {provided_data} format is Unknown, empty or invalid.")
                        sys.exit(1)
            else:
                logging.error(f"The infile {provided_data} does not exist.")
                sys.exit(1)
                
        if not seqdepth_third:                # 当提供的不是三代数据的时候，但是又提供了三代数据的类型，则给出提示，将类型设置为 None
            if seqdepth_type:
                logging.warning(f"The third sequencing data type '{seqdepth_type}' in the [sequencing_depth] section is not required.")
                logging.warning("The 'TGS_type' will be reset to None defaultly.")
                time.sleep(3)
                seqdepth_type = ''
     
        # Check if exactly one data type is provided
        if provided_data_count != 1:
            logging.error("Exactly one of 'NGS_single_end', 'NGS_pair_ends', or 'TGS' must be provided in the [seqdepth] section.")
            exit(1)

        if seqdepth_threads.strip() == '':
            total_cores = os.cpu_count() or 1  # 处理 os.cpu_count() 可能返回 None 的情况
            seqdepth_threads = max(1, min(total_cores - 1, int(total_cores * 0.95)))
            logging.info(f"No thread parameter was provided, default thread has been used instead: {seqdepth_threads} threads.")
        else:
            # 如果不是空字符串，则尝试转换为整数
            try:
                seqdepth_threads = int(float(seqdepth_threads)) if int(float(seqdepth_threads)) <= os.cpu_count() else None       # 用于防止设置的进程数 > 计算机本身的进程数
            except ValueError:
                # 如果转换失败，则使用默认值，并记录错误信息
                total_cores = os.cpu_count() or 1
                seqdepth_threads = max(1, min(total_cores - 1, int(total_cores * 0.95)))
                logging.error(f"No thread parameter was provided, default parameter has been used instead: {seqdepth_threads} threads.")

################################################################################
        # [subconfiguration] Section
        extrsubcon_output_dir_prefix = config['subconfiguration'].get('output_directory_prefix', 'subconfig_repeat-spanned_results').strip() or 'subconfig_repeat-spanned_results'
        extrsubcon_trimmed_length = config['subconfiguration'].getint('flanked_sequence_length', 1000)

        # [mainconfiguration] Section
        extrmaincon_output_dir_prefix = config['mainconfiguration'].get('output_directory_prefix', 'mainconfig_repeat-spanned_results').strip() or 'mainconfig_repeat-spanned_results'
        extrmaincon_trimmed_length = config['mainconfiguration'].getint('flanked_sequence_length', 1000)

        # 从 [check_spanning_reads] 节获取参数
        check_spanning_length = config['check_spanning_reads'].getint('spanning_read_flanking_repeat_length', 1)
        read_number = config['check_spanning_reads'].getint('spanning_read_number', 1)

        try:      #检查亚构型中截取序列长度参数是否合法
            validate_positive_integer(extrsubcon_trimmed_length, 'flanked_sequence_length')
        except ValueError as e:
            logging.error(f"Invalid 'flanked_sequence_length' in [extrsubcon_trimmed_length] section.")

        try:      #检查主构型中截取序列长度参数是否合法
            validate_positive_integer(extrmaincon_trimmed_length, 'flanked_sequence_length')
        except ValueError as e:
            logging.error(f"Invalid 'flanked_sequence_length' in [extrmaincon_trimmed_length] section.")

        try:      #检查跨越重复序列的长度参数是否合法
            validate_positive_integer(check_spanning_length, 'check_spanning_length')
        except ValueError as e:
            logging.error(f"Invalid 'check_spanning_length' in [check_spanning_length] section.")
        
################################################################################
        if mode in ['A', 'C', 'R']:
            # Directory and file checks. 删除之前运行的痕迹
            chromosome_numbers = range(1, count + 1)  # Generate chromosome numbers from 1 to count
            check_files = [
                           f"{extrsubcon_output_dir_prefix}_{project_id}",
                           f"{extrmaincon_output_dir_prefix}_{project_id}",
                           f"{output_dir_prefix}_{project_id}",
                           f"record_for_resume_{project_id}.txt",
                           f"cat_inputfasta_{project_id}.fasta",
                           f"{project_id}",
                           f"repeat_positions_{project_id}.tsv",
                           f"repeat_sequences_{project_id}.fasta",
                           f"repeat_supp_subconfig_recomb_{project_id}.txt",
                           f"ROUSFinder_results_{project_id}",
                           #f"custom_log_{project_id}.txt",        # 前面对日志文件进行了管理，此处无需再次管理
                           f"cat_inputfasta_{project_id}.fasta",
                           ] + [
                           f"repeat_positions_{project_id}_chr{chr_num}.tsv" for chr_num in chromosome_numbers
                           ] + [
                           f"repeat_sequences_{project_id}_chr{chr_num}.fasta" for chr_num in chromosome_numbers
                           ]
            # Check if directories and files exist
            existing_paths = [path for path in check_files if os.path.exists(path)]
            # 单独处理带有通配符的文件名模式  
            wildcard_pattern = f"*_selected_{project_id}.fastq"    # 处理测序文件过滤后的结果文件
            # 使用glob模块检查匹配通配符模式的文件  
            wildcard_matches = glob.glob(wildcard_pattern)  
            # 将存在的通配符匹配文件添加到existing_paths列表中  
            existing_paths.extend(wildcard_matches)  

            if existing_paths:  # 如果列表非空，表示至少有一个文件或目录存在
                if args.redo:  
                    # 删除相关目录和文件  
                    for path in existing_paths:  
                        if os.path.isdir(path):  
                            logging.info(f"Delete previous folder: {path}.")  
                            shutil.rmtree(path)    # 删除目录及其内容  
                        elif os.path.isfile(path):  
                            logging.info(f"Delete previous file: {path}.")  
                            os.remove(path)  # 删除文件

                    logging.info("Previous results have been deleted. Starting fresh calculations.\n")
                    time.sleep(5)

                elif not args.resume and not args.redo:
                    # 提示用户决定
                    # Format the existing paths into a multi-line string
                    existing_paths_str = "\n".join(existing_paths)
                    logging.info(f"ATTENTION: Find previous folder(s)/file(s): {existing_paths_str}.")
                    print("You can eliminate this issue through one of the following methods:")
                    print("1. Choose '-resume' to continue from a previous run!")                # record_for_resume_{project_id}.txt' to resume from a previous run!")
                    print("2. Choose '-redo' to delete previous results and re-calculate!")
                    print(f"3. Reset the 'project_id' in the '{args.config}' file!\n")
                    sys.exit(1)

################################################################################
        # Get the directory of the current script
        dir_path = os.path.dirname(os.path.realpath(__file__))
        rous_output_files = glob.glob(f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}*")    # rous_output_file_prefix 输入文件前缀作为输出结果的文件前缀

######################################
        #此处为 resume 模式设置一个判断条件，即处于 resume 模式的时候，跳过某些运行程序，以节约时间
        record_file = f"record_for_resume_{project_id}.txt"          # 下面的程序中，record_file用于记录处理过的序列
        resume_logic = False
        # 检查三个条件：文件存在、文件非空、用户启用了resume模式
        if (os.path.exists(record_file) and os.path.getsize(record_file) > 0 and args.resume and os.path.exists(extrsubcon_output_dir_prefix + "_" + project_id) and os.path.exists(extrmaincon_output_dir_prefix + "_" + project_id) and os.path.exists("ROUSFinder_results_" + project_id)):
            resume_logic = True
######################################

        if count > 1:
            cat_inputfasta_path, headers, lengths = concatenate_fasta(inputfasta, f"cat_inputfasta_{project_id}.fasta")
        if mode == 'A':
            # 创建一个要检查的文件名列表，删除之前运行的痕迹
            check_files = [
                f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_5CT.tsv",
                f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt",
                f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt",
                f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep.fasta"
                ]
            # 筛选出存在的文件
            existing_rous_files = [file for file in rous_output_files if os.path.isfile(file) and file in check_files]
            if len(existing_rous_files) > 0 and not args.resume:
                logging.warning(f"Warning: Files {existing_rous_files} in ROUSFinder_results_{project_id} will be OVERWRITTEN.")
                print("Press 'Ctrl + C' to stop the program now if you do NOT wish to OVERWRITE the above files.")
                try:
                    # 等待10秒以给用户取消的机会
                    logging.info("Waiting for 10 seconds before continuing...\n")
                    sys.stdout.flush()  # 确保消息立即打印
                    time.sleep(10)
                except KeyboardInterrupt:
                    logging.info("Program terminated by the user.")
                    sys.exit(1)
            elif len(existing_rous_files) > 0 and args.redo:
                logging.info(f"Delete previous file: {existing_rous_files}")
                os.remove(existing_rous_files)

################################################################################
            # Execute ROUSFinder regardless of existing files
            if count == 1 and not resume_logic:
                run_command(["python", os.path.join(dir_path, "ROUSFinder2.0.py"), "-i", inputfasta, "-o", rous_output_file_prefix + "_" + project_id, "-m", str(rous_repeat_minlength), "-rew", str(rous_reward), "-pen", str(rous_penalty), "-id", project_id])
            if count > 1 and not resume_logic:
                run_command(["python", os.path.join(dir_path, "ROUSFinder2.0.py"), "-i", cat_inputfasta_path, "-o", rous_output_file_prefix + "_" + project_id, "-m", str(rous_repeat_minlength), "-rew", str(rous_reward), "-pen", str(rous_penalty), "-id", project_id])

        elif mode == 'R':
           # Check if the files exist and the number of files found
            check_files = [
                f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt",
                f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep.fasta"
                ]
            # 筛选出存在的文件
            existing_rous_files = [file for file in rous_output_files if os.path.isfile(file)]

            if len(existing_rous_files) > 0 and not args.resume:
                logging.warning(f"Files {existing_rous_files} in folder 'ROUSFinder_results_{project_id}' will be OVERWRITTEN.")
                print("Press 'Ctrl + C' to stop the program now if you do NOT wish to be OVERWRITTEN.")
                try:
                    # 等待10秒以给用户取消的机会
                    logging.info("Waiting for 10 seconds before continuing...")
                    sys.stdout.flush()  # 确保消息立即打印
                    time.sleep(10)
                except KeyboardInterrupt:
                    logging.exception("Program terminated by the user.")
                    sys.exit(1)
            elif len(existing_rous_files) > 0 and args.redo:
                logging.info(f"Delete previous file: {existing_rous_files}")
                os.remove(existing_rous_files)
                
            # 利用 ROUSFinder 生成重复序列信息的5列表 
            if count == 1 and not resume_logic:
                run_command(["python", os.path.join(dir_path, "ROUSFinder2.0.py"), "-i", inputfasta, "-o", rous_output_file_prefix + "_" + project_id, "-m", str(rous_repeat_minlength), "-rew", str(rous_reward), "-pen", str(rous_penalty), "-id", project_id])
            elif count > 1 and not resume_logic:
                run_command(["python", os.path.join(dir_path, "ROUSFinder2.0.py"), "-i", cat_inputfasta_path, "-o", rous_output_file_prefix + "_" + project_id, "-m", str(rous_repeat_minlength), "-rew", str(rous_reward), "-pen", str(rous_penalty), "-id", project_id])
                
                
            # 对 ROUSFinder 生成的复序列信息的5列表 进行重复单元长度的筛选
            ROUSFinder_results_file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt" 
            if os.path.isfile(ROUSFinder_results_file_path) and not resume_logic:  
                with open(ROUSFinder_results_file_path, 'r') as file:  
                    lines = file.readlines()  
                    line_count = len(lines)  
              
                    # 如果文件只有一行（表头），则报错并退出  
                    if line_count == 1:  
                        logging.warning("The ROUSFinder results file contains only a header row; no data rows were found.")  
                        sys.exit(1)
                        
                    #### 对 ROUSFinder 的结果根据 输入的 repeat 的长度信息进行删减 ####
                    #### 在 ROUSFinder_results_file_path 文件里直接修改
                    filtered_lines = []  
                    with open(ROUSFinder_results_file_path, 'r') as file:  
                        lines = file.readlines()  
                        # 保留表头  
                        filtered_lines.append(lines[0])  
                        for line in lines[1:]:  # 跳过表头  
                            parts = line.strip().split('\t')  # 假设字段之间是由空白字符分隔的  
                            # 检查length是否在repeat_lengths列表中  
                            if parts[1] in map(str, repeat_lengths):  # 将repeat_lengths中的元素转换为字符串类型进行比较  
                                filtered_lines.append(line)  
                    # 将筛选后的内容写回原文件  
                    with open(ROUSFinder_results_file_path, 'w') as file:  
                        file.writelines(filtered_lines)  
                        
                        
        # 处理 mode == "R"，染色体数量超过 1，追溯 start 和 end 的原始位置，增加 chromosome 编号信息
        if mode == "R" and count >1 and not resume_logic:
            if os.path.isfile(f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt"):
                file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt"
                
                default_headers = ['fragment_id', 'length', 'start','end','direction']
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path):   
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 

                # 读取文件，跳过第一行
                with open(file_path, 'r', newline='', encoding='utf-8') as f:
                    lines = f.readlines()
                # 将第一行替换为default_headers
                lines[0] = '\t'.join(default_headers) + '\n'
                # 将修改后的内容写回文件
                with open(file_path, 'w', newline='', encoding='utf-8') as f:
                    f.writelines(lines)

                df = pd.read_csv(file_path, sep='\t')  # 重新读取文件，带标题行
                df = df.dropna(how='all')    # 去掉空白行

                if len(df) == 0:
                    logging.warning(f"The DataFrame {file_path} is empty!!")
                else:
                    # 应用函数并分别更新列  
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = np.where(df['adjusted_start'] > 0, df['adjusted_start'], 1)
                    df['end'] = np.where(df['adjusted_end'] > 0, df['adjusted_end'], 1)
                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])    # 删除两列
                    # 将chromosome列插入到 direction 列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('direction') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  

                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]

                # 保存修改后的 DataFrame  
                file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt" 
                df.to_csv(file_path, sep='\t', index=False)              

                # 根据重复序列的名字，清除多余的重复序列单元。分三步完成。
                repeat_file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep.fasta"
                # 步骤1: 读取_rep_table.txt文件中的fragment_id列，并去除重复项 fasta_file_path
                df = pd.read_csv(file_path, sep='\t')
                unique_fragment_ids = df['fragment_id'].drop_duplicates()

                # 步骤2: 使用去重后的fragment_id作为header，从_rep.fasta文件中读取相应的序列
                with open(repeat_file_path, 'r') as fasta_file:
                    sequences = {}
                    current_id = None
                    for line in fasta_file:
                        if line.startswith('>'):  # FASTA header line
                            if current_id is not None:
                                sequences[current_id] = ''.join(sequence)
                            current_id = line[1:].strip()  # Remove the '>' and strip whitespace
                            sequence = []
                        else:
                            sequence.append(line.strip())  # Append the sequence lines

                   # Add the last sequence
                    if current_id is not None:
                        sequences[current_id] = ''.join(sequence)

                # 步骤3: 将读取的序列覆盖写入_rep.fasta文件，只包含去重后的fragment_id对应的序列
                with open(repeat_file_path, 'w') as fasta_file:
                    for fragment_id in unique_fragment_ids:
                        if fragment_id in sequences:
                            fasta_file.write(f">{fragment_id}\n")
                            fasta_file.write(f"{sequences[fragment_id]}\n")

################################################################################
        #### 开始处理前，增加一个判断，也就是有没有重复序列单元被找到，未找到 则报错退出
        #### 检查 ROUSFinder 的结果是不是为空
        if mode == 'A':
            ROUSFinder_results_file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt" 
            if os.path.isfile(ROUSFinder_results_file_path):  
                with open(ROUSFinder_results_file_path, 'r') as file:  
                    lines = file.readlines()  
                    line_count = len(lines)  
              
                    # 如果文件只有一行（表头），则报错并退出  
                    if line_count == 1:  
                        logging.warning("The ROUSFinder results file contains only a header row; no data rows were found.")  
                        sys.exit(1)
                        
                    #### 对 ROUSFinder 的结果根据 输入的 repeat 的长度信息进行删减 ####
                    #### 在 ROUSFinder_results_file_path 文件里直接修改
                    filtered_lines = []  
                    with open(ROUSFinder_results_file_path, 'r') as file:  
                        lines = file.readlines()  
                        # 保留表头  
                        filtered_lines.append(lines[0])  
                        for line in lines[1:]:  # 跳过表头  
                            parts = line.strip().split('\t')  # 假设字段之间是由空白字符分隔的  
                            # 检查length是否在repeat_lengths列表中  
                            if parts[1] in map(str, repeat_lengths):  # 将repeat_lengths中的元素转换为字符串类型进行比较  
                                filtered_lines.append(line)  
                    # 将筛选后的内容写回原文件  
                    with open(ROUSFinder_results_file_path, 'w') as file:  
                        file.writelines(filtered_lines)  
            else:  
                logging.error(f"The ROUSFinder result file {ROUSFinder_results_file_path} was not found.")  
                sys.exit(1)
        
################################################################################
        if mode == 'A' and not resume_logic:
            # 对ROUSFinder结果进行校准前，检查文件内是否存在重复序列的数据
            input_file = f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_table.txt"
            column_names = ['fragment_id', 'length', 'start', 'end', 'direction']
            df = pd.read_csv(input_file, delimiter='\t', header=None, names=column_names, skiprows=1)
            if len(df) == 0:
                logging.error(f"No repeat info was found in the file {input_file}.")
                print("Please reset the 'repeat_length' in the [ROUSFinder] section.")
                print("Or check if there are repeats in your provided genome.\n")
                sys.exit(1)
            # 对ROUSFinder结果进行校准
            run_command(["python", os.path.join(dir_path, "calibrate_ROUSFinder_results.py"), "-i", f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_table.txt", "-fl", str(calibrate_flexibility), "-o", f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id, "-t", str(seqdepth_threads)])
            # 将程序校准后的ROUSFinder结果转换为5CT
            run_command(["python", os.path.join(dir_path, "ROUSFinder_to_5CT.py"), "-i", f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_calibration_table.txt", "-l", str(sequence_lengths), "-o", f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id, "-t", str(seqdepth_threads)])
        
        if mode == 'C' and not resume_logic:
            logging.info("ATTENTION: Please ensure that the position information of repeats is consistent with the genome sequence!!")
            print()
            time.sleep(5)

            if not os.path.isdir(f"ROUSFinder_results_{project_id}/"):
                 # Create a new folder
                os.makedirs(f"ROUSFinder_results_{project_id}/")

            def validate_columns_in_manually_calibrate(manually_calibrate, count):
                """
                Validates the number of columns in the manually_calibrate file and checks the header.
                Parameters:
                manually_calibrate (str): Path to the manually_calibrate file.
                count (int): The count of sequences in the FASTA file.
                Raises:
                ValueError: If the number of columns does not match the expected count or header is incorrect.
                """
                try:
                    # Read the first line of the manually_calibrate file
                    with open(manually_calibrate, 'r') as file:
                        first_line = file.readline().strip()
        
                    # Split the first line into columns
                    columns = first_line.split('\t')
        
                    # Define the expected headers
                    if count == 1:
                        expected_header = "fragment_id\tlength\tstart\tend\tdirection\n"
                    elif count > 1:
                        expected_header = "fragment_id\tlength\tstart\tend\tdirection\tchromosome\n"
        
                    # Check the number of columns based on the count
                    if count == 1 and (len(columns) != 5 or first_line != expected_header.strip()):
                        raise ValueError(f"Only one chromosome was found in 'inputfasta' file. The 'calibrated_repeat_file' file must have 5 columns, but found 6 columns.\n")
                    elif count > 1 and (len(columns) != 6 or first_line != expected_header.strip()):
                        raise ValueError(f"More than two chromosomes were found. The manually_calibrate file must have 6 columns with header '{expected_header}', but found {len(columns)} columns with header '{first_line}'.")
        
                except FileNotFoundError:
                    print(f"The file {manually_calibrate} does not exist.\n")
                    exit(1)
                except ValueError as ve:
                    print(ve)
                    exit(1)
                except Exception as e:
                    print(f"An error occurred: {e}\n")
                    exit(1)
            
            # 检测 人工输入的列表与fasta序列是否一致，通过序列的条数进行检测，count=1,必须是5列，count>2，必须是6列
            validate_columns_in_manually_calibrate(manually_calibrate, count)

            # 对手工矫正的重复序列信息进行检测，使其符合一系列规则，产生的结果存放在 adjusted_manual_calibrated_list.tsv 中
            # count =1 或者 > 1，都产生这个文件
            run_command(["python", os.path.join(dir_path, "check_manual_calibrated_repeats.py"), "-c", str(count), "-l", manually_calibrate, "-f", inputfasta])    # 检查手工修正重复序列信息列表
            
            # 为手工修改的校正结果，放进ROUSFinder_results
            # 定义源文件路径和目标目录
            source_file = "adjusted_manual_calibrated_list.tsv"
            target_dir = f"ROUSFinder_results_{project_id}"

            # 确保目标目录存在，如果不存在则创建
            os.makedirs(target_dir, exist_ok=True)
            # 移动文件
            shutil.move(source_file, os.path.join(target_dir, os.path.basename(source_file)))
            
            manually_calibrate_adj = f"ROUSFinder_results_{project_id}/adjusted_manual_calibrated_list.tsv"     # 更改手工校准文件的名字

            # 提取输入文件的前缀
            dir_path_manual, file_name_manual = os.path.split(manually_calibrate_adj)       # 拆解为路径和文件名
            prefix_manual = os.path.splitext(file_name_manual)[0]    # 将文件名拆解为前缀和后缀
            
            # Read the file，更改表头
            with open(manually_calibrate_adj, 'r') as file:
                lines = file.readlines()
                if len(lines) <= 1:              ## 检查输入的repeat信息文件是否为空
                    raise ValueError("The file is empty, only contains the header, or is missing the header. Please provide a valid file with data.")
                    exit(1)
            # Replace the first line with the new header，替换原有表头，以免用户提供的有误
            #lines[0] = "fragment_id\tlength\tstart\tend\tdirection\n"
            
            manually_calibrate_reset = f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_calibration_table.txt"          ##更换表头之后，保存为calibration_table.txt
            # 将修改后的lines写入新的文件路径，达到更改表头的目的
            with open(manually_calibrate_reset, 'w') as file:
                file.writelines(lines)               # 完成由 manually_calibrate_adj 向 manually_calibrate_reset 的转变
            # 以上 manually_calibrate_adj 对应的文件 adjusted_manual_calibrated_list.tsv
            # 以上 manually_calibrate_reset 对应的文件 adjusted_manual_calibrated_list_rep_calibration_table.txt
            
            # 检查前一次运行的痕迹
            rous_output_files = glob.glob(f"ROUSFinder_results_{project_id}/{prefix_manual}*")
            # 创建一个要检查的文件名列表
            check_files = [
                f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_5CT.tsv",     # 此处是手工校准的结果，加上project_id
                f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_calibration_table.txt"        # 所以ROUSFinder_results中只应含有rep_5CT.tsv
                ]
            # 筛选出存在的文件
            existing_rous_files = [file for file in rous_output_files if os.path.isfile(file) and os.path.basename(file) in check_files]
            if len(existing_rous_files) > 0 and not args.resume:
                logging.warning(f"Files {existing_rous_files} in ROUSFinder_results_{project_id} will be OVERWRITTEN.")
                print("Press 'Ctrl + C' to STOP the program now. Or it will be OVERWRITTEN after 10s.")
                try:
                    # 等待10秒以给用户取消的机会
                    logging.info("Waiting for 10 seconds before continuing...\n")
                    sys.stdout.flush()  # 确保消息立即打印
                    time.sleep(10)
                except KeyboardInterrupt:
                    logging.exception("Program terminated by the user.")
                    sys.exit(1)
            elif len(existing_rous_files) > 0 and args.redo:
                logging.info(f"Delete previous file: {existing_rous_files}")
                os.remove(existing_rous_files)

            max_start = 0    # 检查输入的手工校准的重复序列文件与fasta序列文件是否一致。通过fasta序列的长度与 manually_calibrate 的位置是否匹配来进行检查
            max_end = 0
            for line in lines[1:]:
                start, end = map(int, line.strip().split('\t')[2:4])
                max_start = max(max_start, start)
                max_end = max(max_end, end)
            max_position = max(max_start, max_end)
            if max_position > sequence_lengths:
                logging.error(f"Some duplicates in the {manually_calibrate} are located outside the genome provided in the sequence {inputfasta}.")
                exit(1)

            if count >1:
                cat_inputfasta_path, headers, lengths = concatenate_fasta(inputfasta, f"cat_inputfasta_{project_id}.fasta")      # 再次创建pseudo-genome

            # 将手工校准后的ROUSFinder结果转换为5CT
            run_command(["python", os.path.join(dir_path, "ROUSFinder_to_5CT.py"), "-i", manually_calibrate_reset, "-l", str(sequence_lengths), "-o", f"ROUSFinder_results_{project_id}/" + prefix_manual + "_" + project_id, "-t", str(seqdepth_threads)])

################################################################################
        if mode == 'A' and not resume_logic:
            if count == 1:
                # subconfigure截取reference
                run_command(["python", os.path.join(dir_path, "extrsubcon.py"), "-r", inputfasta, "-rf", f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_5CT.tsv", "-tl", str(extrsubcon_trimmed_length), "-o", extrsubcon_output_dir_prefix + "_" + project_id, "-rc", complementary_chain, "-gt", genome_type])
                # mainconfigure截取reference
                run_command(["python", os.path.join(dir_path, "extrmaincon.py"), "-r", inputfasta, "-rf", f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_5CT.tsv", "-tl", str(extrmaincon_trimmed_length), "-o", extrmaincon_output_dir_prefix + "_" + project_id, "-rc", complementary_chain, "-gt", genome_type])
            if count > 1:  
                run_command(["python", os.path.join(dir_path, "extrsubcon.py"), "-r", cat_inputfasta_path, "-rf", f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_5CT.tsv", "-tl", str(extrsubcon_trimmed_length), "-o", extrsubcon_output_dir_prefix + "_" + project_id, "-rc", complementary_chain, "-gt", genome_type])
                run_command(["python", os.path.join(dir_path, "extrmaincon.py"), "-r", cat_inputfasta_path, "-rf", f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_5CT.tsv", "-tl", str(extrmaincon_trimmed_length), "-o", extrmaincon_output_dir_prefix + "_" + project_id, "-rc", complementary_chain, "-gt", genome_type])

        if mode == 'C' and not resume_logic:
            if count == 1:
                # subconfigure截取reference
                run_command(["python", os.path.join(dir_path, "extrsubcon.py"), "-r", inputfasta, "-rf", f"ROUSFinder_results_{project_id}/" + prefix_manual + "_" + project_id + "_rep_5CT.tsv", "-tl", str(extrsubcon_trimmed_length), "-o", extrsubcon_output_dir_prefix + "_" + project_id, "-rc", complementary_chain, "-gt", genome_type])
                # mainconfigure截取reference
                run_command(["python", os.path.join(dir_path, "extrmaincon.py"), "-r", inputfasta, "-rf", f"ROUSFinder_results_{project_id}/" + prefix_manual + "_" + project_id + "_rep_5CT.tsv", "-tl", str(extrmaincon_trimmed_length), "-o", extrmaincon_output_dir_prefix + "_" + project_id, "-rc", complementary_chain, "-gt", genome_type])
            if count > 1:
                run_command(["python", os.path.join(dir_path, "extrsubcon.py"), "-r", cat_inputfasta_path, "-rf", f"ROUSFinder_results_{project_id}/" + prefix_manual + "_" + project_id + "_rep_5CT.tsv", "-tl", str(extrsubcon_trimmed_length), "-o", extrsubcon_output_dir_prefix + "_" + project_id, "-rc", complementary_chain, "-gt", genome_type])
                run_command(["python", os.path.join(dir_path, "extrmaincon.py"), "-r", cat_inputfasta_path, "-rf", f"ROUSFinder_results_{project_id}/" + prefix_manual + "_" + project_id + "_rep_5CT.tsv", "-tl", str(extrmaincon_trimmed_length), "-o", extrmaincon_output_dir_prefix + "_" + project_id, "-rc", complementary_chain, "-gt", genome_type])

        # 用于输出的结果更好看
        if mode in ['A', 'C']:
            print()
            logging.info("@@@@@@@@@@@@ Start to detect genome recombination mediated by repeats. @@@@@@@@@@@@")
            logging.info("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            time.sleep(2)

        if mode in ['A', 'C']:         # 处理读取文件后的文件名，便于后续的数据处理
            # Determine the sequencing data type and its value
            #seq_data_flag = '-' + provided_data_type
            if provided_data_type == 'NGS_pair_ends':
                seq_data_value = provided_data_paths
            if provided_data_type == 'NGS_single_end' or provided_data_type == 'TGS':
                seq_data_value = provided_data
            # Ensure the final result directory exists
            final_output_dir = f"{output_dir_prefix}_{project_id}"
            # Check if the directory exists
            if not os.path.isdir(final_output_dir):
                # If it doesn't exist, create it
                os.makedirs(final_output_dir, exist_ok=True)

################################################################################
#### 提前处理测序数据：fastq_2_fasta，reads_filter
#### blast 只能接受fasta格式，minimap2 和 bwa 可以接受 fastq 和 fasta 两种格式，所以，不管是否进行过滤，均将序列文件转为 fasta 格式
#### 利用blast过滤reads 二代 三代 均要筛选 即使后面使用blast寻找跨越 repeat 的reads, 也同样进行reads的筛选
        if mode in ['A', 'C']:
            fastq2fasta = False       # 用于删除 fastq 转化来的 fasta
            if provided_data_type in ['NGS_single_end', 'TGS']:
                fastq_prefix = os.path.splitext(os.path.basename(seq_data_value))[0]           # 提取测序数据文件的前缀
                file_format = check_file_format_efficient(seq_data_value)                      # 检查输入文件的数据格式
                if filter_reads in ['Y','YES']:            # 开始过滤reads，提高后面比对的速度。
                    seq_data_value_filter = f"{fastq_prefix}_{project_id}_filtered.fasta"            # 用于 resume 模式
                    if not os.path.exists(seq_data_value_filter) or not os.path.exists(record_file): 
                        if file_format == 'FASTQ':         # 当输入文件为fastq格式，进行格式转换
                            print()
                            logging.info("Converting fastq to fasta ...") 
                            fasta_one = fastq_to_fasta_multiprocess(seq_data_value, fastq_prefix + ".fasta", num_processes=seqdepth_threads)
                            fastq2fasta = True
                        elif file_format == 'FASTA': 
                            fasta_one = seq_data_value                     # 为后面读取fasta设置路径
                        
                        fasta_prefix = os.path.splitext(os.path.basename(fasta_one))[0]          # 提取前缀, 此时的文件一定是fasta格式，与之前的一句重复
                        # 过滤reads
                        subprocess.run(["bash", os.path.join(dir_path, "blast.filter.sh"), inputfasta, fasta_one, f"{fasta_prefix}_{project_id}", str(seqdepth_threads), str(blast_evalue)], check=True)     
                        seq_data_value_filter = f"{fasta_prefix}_{project_id}_filtered.fasta"         # 重新定义用于比对的reads文件，为以上的"blast.filter.sh"产生的文件名
                
                if filter_reads in ['N','NO']:
                    fasta_one = os.path.join(os.getcwd(), fastq_prefix + ".fasta")
                    if not os.path.exists(fasta_one) or not os.path.exists(record_file):
                        if file_format == 'FASTQ':                         # 当输入文件为fastq格式，进行格式转换
                            logging.info("Converting fastq to fasta ...") 
                            fasta_one = fastq_to_fasta_multiprocess(seq_data_value, fastq_prefix + ".fasta", num_processes=seqdepth_threads)
                            fastq2fasta = True
                        elif file_format == 'FASTA':
                            fasta_one = seq_data_value   
                        
            if provided_data_type == 'NGS_pair_ends':
                fastq1, fastq2 = seq_data_value 
                
                fastq1_prefix = os.path.splitext(os.path.basename(fastq1))[0]       # 提取测序数据文件的前缀
                file_format = check_file_format_efficient(fastq1)      # 检查输入文件的数据格式
                if filter_reads in ['Y','YES']:                        # 开始过滤 reads，提高后面比对的速度。
                    fasta_pair1 = os.path.join(os.getcwd(), fastq1_prefix + ".fasta")
                    fasta1_prefix = os.path.splitext(os.path.basename(fasta_pair1))[0]          # 提取前缀
                    fasta_pair1_filter = f"{fasta1_prefix}_{project_id}_filtered.fasta"
                    if not os.path.exists(fasta_pair1_filter) or not os.path.exists(record_file):
                        if file_format == 'FASTQ':                                 # 当输入文件为 fastq 格式，进行格式转换
                            print()
                            logging.info("Converting fastq (5'-end) to fasta ...")
                            fasta_pair1 = fastq_to_fasta_multiprocess(fastq1, fastq1_prefix + ".fasta", num_processes=seqdepth_threads)
                            fastq2fasta = True
                        elif file_format == 'FASTA': 
                            fasta_pair1 = fastq1             # 为后面读取fasta设置路径

                        # 过滤reads
                        subprocess.run(["bash", os.path.join(dir_path, "blast.filter.sh"), inputfasta, fasta_pair1, f"{fasta1_prefix}_{project_id}", str(seqdepth_threads), str(blast_evalue)], check=True)     
                        fasta_pair1_filter = f"{fasta1_prefix}_{project_id}_filtered.fasta"    # 重新定义用于比对的reads文件，为以上的"blast.filter.sh"产生的文件名
                
                if filter_reads in ['N','NO']:
                    fasta1 = os.path.join(os.getcwd(), fastq1_prefix + ".fasta")
                    if not os.path.exists(fasta1) or not os.path.exists(record_file):
                        if file_format == 'FASTQ':                                   # 当输入文件为fastq格式，进行格式转换
                            logging.info("Converting fastq (5'-end) to fasta ...")
                            fasta_pair1 = fastq_to_fasta_multiprocess(fastq1, fastq1_prefix + ".fasta", num_processes=seqdepth_threads)
                            fastq2fasta = True
                        elif file_format == 'FASTA':
                            fasta_pair1 = fastq1    
                        fastq1 = fasta_pair1
                    
                fastq2_prefix = os.path.splitext(os.path.basename(fastq2))[0]         # 提取前缀
                file_format = check_file_format_efficient(fastq2)          # 检查输入文件的数据格式
                fasta_pair2_filter = f"{fasta1_prefix}_{project_id}_filtered.fasta"
                if filter_reads in ['Y','YES']:
                    fasta_pair2 = os.path.join(os.getcwd(), fastq2_prefix + ".fasta")
                    fasta2_prefix = os.path.splitext(os.path.basename(fasta_pair2))[0]           # 提取前缀
                    fasta_pair2_filter = f"{fasta2_prefix}_{project_id}_filtered.fasta"
                    if not os.path.exists(fasta_pair2_filter) or not os.path.exists(record_file):
                        if file_format == 'FASTQ':                         # 当输入文件为fastq格式，进行格式转换
                            logging.info("Converting fastq (3'-end) to fasta ...")
                            fasta_pair2 = fastq_to_fasta_multiprocess(fastq2, fastq2_prefix + ".fasta", num_processes=seqdepth_threads)
                            fastq2fasta = True
                        elif file_format == 'FASTA':
                            fasta_pair2 = fastq2                                # 为后面读取fasta设置路径

                        subprocess.run(["bash", os.path.join(dir_path, "blast.filter.sh"), inputfasta, fasta_pair2, f"{fastq2_prefix}_{project_id}", str(seqdepth_threads), str(blast_evalue)], check=True)     # 过滤reads
                        fasta_pair2_filter = f"{fasta2_prefix}_{project_id}_filtered.fasta"             # 重新定义用于比对的reads文件，为以上的"blast.filter.sh"产生的文件名
                
                if  filter_reads in ['N','NO']:
                    fasta2 = os.path.join(os.getcwd(), fastq2_prefix + ".fasta")
                    if not os.path.exists(fasta2) or not os.path.exists(record_file):
                        if file_format == 'FASTQ':                                      # 当输入文件为fastq格式，进行格式转换
                            logging.info("Converting fastq (3'-end) to fasta ...")
                            fasta_pair2 = fastq_to_fasta_multiprocess(fastq2, fastq2_prefix + ".fasta", num_processes=seqdepth_threads)
                            fastq2fasta = True
                        elif file_format == 'FASTA':
                            fasta_pair2 = fastq2    
                        fastq2 = fasta_pair2
                    
################################################################################
        if mode in ['A', 'C']:
            # 初始化colorama  
            init(autoreset=True)  

            #### subconfiguration - Process each FASTA file sequentially
            subcon_fasta_files = glob.glob(f"{extrsubcon_output_dir_prefix}_{project_id}/*.fasta")

            subcon_fasta_total = len(subcon_fasta_files)    # 次要构型中需要检测的重复单元个数
            print()
            logging.info(f"@@@@@@@@@@ Found {subcon_fasta_total} repeat pairs that may mediate genomic recombination! @@@@@@@@@@")
            time.sleep(5)

            print()
            logging.info(f"%%%%%%%%%%%%% Start to process {subcon_fasta_total} repeat units in the subconfiguration! %%%%%%%%%%%%%")
            logging.info(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            time.sleep(2)
            
            file_path_sub = os.path.join(extrsubcon_output_dir_prefix + "_" + project_id, "filter_subconfig_spanned_results.tsv")
            output_final_sub = os.path.dirname(file_path_sub)
            if not os.path.isdir(output_final_sub):
                os.makedirs(output_final_sub, exist_ok=True)

            header = "sub_trimmed_seq_name\tsub_trimmed_seq_length(bp)\tsub_spanning_read_number"
            # 检查文件是否存在或为空，然后写入表头
            if not os.path.isfile(file_path_sub) or os.path.getsize(file_path_sub) == 0:
                with open(file_path_sub, 'w') as file:
                    file.write(header + '\n')

            # 逐一处理截取的序列
            record_file = f"record_for_resume_{project_id}.txt"    #记录处理过的序列
            # 检查记录文件是否存在，如果存在则读取其中的内容
            if os.path.isfile(record_file):
                with open(record_file, "r") as f:
                    processed_files = f.read().splitlines()
            else:
                processed_files = []
                
            every_subcon_count = 1    # 初始化计数器  
            for fasta_file in subcon_fasta_files:    ##################### 这里开始处理次要构型 subconfiguration

                ordinal_suffix = get_ordinal_suffix(every_subcon_count)  
                ordinal_suffix_basename = os.path.basename(ordinal_suffix)
                # 使用colorama高亮显示计数器和序数词后缀  
                highlighted_text = f"{Fore.RED}{every_subcon_count}{ordinal_suffix_basename}{Style.RESET_ALL}"
                print()
                logging.info(f"**** Start to process the {highlighted_text} sequence: {os.path.basename(fasta_file)}. ****")  
                time.sleep(2)
                every_subcon_count += 1    # 计数器递增

                if os.path.basename(fasta_file) in processed_files and args.resume:
                    logging.info(f"Skipping already processed sequence: {fasta_file}")
                    # 在继续运行模式下(-resume)，如果发现序列之前已经处理过了，则删除掉
                    os.remove(fasta_file) if os.path.isfile(fasta_file) and redundant_intermediate_results == "D" else None  
                    continue
                    
                logging.info(f"Process trimmed seq {os.path.basename(fasta_file)} for subconfig identification!!")

                # 再次确定截取序列的长度，因为截取的时候，重复序列两端剩下的序列长度可能小于设定的截取序列
                extrsubcon_trimmed_length_reset = extract_trimmed_length_from_filename(fasta_file)
                if not extrsubcon_trimmed_length_reset:
                    logging.warning(f"Warning: Could not determine flanked length from file name {fasta_file}. Skipping.")
                    continue
                else:
                    extrsubcon_trimmed_length = extrsubcon_trimmed_length_reset

                if extrsubcon_trimmed_length < check_spanning_length:
                    logging.info("The reset 'flanked_sequence_length' is shorter than 'spanning_read_flanking_repeat_length' in subconfiguration.")
                    logging.info(f"The 'spanning_read_flanking_repeat_length' is also reset to {extrsubcon_trimmed_length}.")
                    check_spanning_length = extrsubcon_trimmed_length

                # Run seqdepth.py for each file
                map_output, _ = os.path.splitext(fasta_file)   
                if seqdepth_alignment_software != 'blast':
                    if provided_data_type == 'NGS_single_end':
                        if filter_reads in ['Y','YES']:        ############################################## 采用过滤后的reads，提高后面比对的速度。
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-single", seq_data_value_filter, "-output", map_output, "-threads", str(seqdepth_threads)] 
                        else:
                            # 不进行reads过滤，seq_data_value可以直接使用，不用区分 fasta/fastq
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-single", seq_data_value, "-output", map_output, "-threads", str(seqdepth_threads)] 
                        
                    if provided_data_type == 'NGS_pair_ends':
                        if filter_reads in ['Y','YES']:        ############################################## 采用过滤后的reads，提高后面比对的速度。
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-pair", fasta_pair1_filter, fasta_pair2_filter, "-output", map_output, "-threads", str(seqdepth_threads)]
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-pair", fastq1, fastq2, "-output", map_output, "-threads", str(seqdepth_threads)]

                    if provided_data_type == 'TGS' and seqdepth_type == 'ont':       # 处理 Nanopore 三代数据
                        if filter_reads in ['Y','YES']:        ############################################## 采用过滤后的reads，提高后面比对的速度。
                            # 使用过滤后的reads seq_data_value，
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", seq_data_value_filter, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", "ont"] 
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", seq_data_value, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", "ont"]

                    if provided_data_type == 'TGS' and seqdepth_type == 'pacbio':           # 处理 pacbio 三代数据
                        if filter_reads in ['Y','YES']:        ############################################## 采用过滤后的reads，提高后面比对的速度。
                            # 使用过滤后的reads seq_data_value，
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", seq_data_value_filter, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", "pacbio"] 
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", seq_data_value, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", "pacbio"]
                        
                ####### 采用 blast 软件 计算测序深度
                if seqdepth_alignment_software == 'blast':
                    if provided_data_type == 'NGS_single_end':           ##处理二代单端数据
                        if filter_reads in ['Y','YES']:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-single", seq_data_value_filter, "-output", map_output, "-threads", str(seqdepth_threads), "-spanning_length", str(check_spanning_length), "-seqdepth_type", seqdepth_type] 
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-single", fasta_one, "-output", map_output, "-threads", str(seqdepth_threads), "-spanning_length", str(check_spanning_length), "-seqdepth_type", seqdepth_type] 

                    if provided_data_type == 'TGS':              ##处理三代数据
                        if filter_reads in ['Y','YES']:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", seq_data_value_filter, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", seqdepth_type, "-spanning_length", str(check_spanning_length)] 
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", fasta_one, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", seqdepth_type, "-spanning_length", str(check_spanning_length)] 

                    if provided_data_type == 'NGS_pair_ends':           ##处理二代双端数据
                        if filter_reads in ['Y','YES']:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-pair", fasta_pair1_filter, fasta_pair2_filter, "-output", map_output, "-threads", str(seqdepth_threads), "-spanning_length", str(check_spanning_length), "-seqdepth_type", seqdepth_type]
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-pair", fasta_pair1, fasta_pair2, "-output", map_output, "-threads", str(seqdepth_threads), "-spanning_length", str(check_spanning_length), "-seqdepth_type", seqdepth_type]

                run_command(command)
           
                # Immediately after, run check_subconfig_spanning_reads.py for the output
                map_results_folder = f"{map_output}_results"
                if os.path.isdir(map_results_folder):
                    output_sam = os.path.join(map_results_folder, "output.sam")
                    run_command(["python", os.path.join(dir_path, "check_subconfig_spanning_reads.py"), "-s", output_sam, "-i", fasta_file, "-o", map_results_folder, "-cl", str(check_spanning_length), "-tl", str(extrsubcon_trimmed_length), "-rn", str(read_number), "-mr", redundant_intermediate_results, "-a", seqdepth_alignment_software])    # 产生repeat-spanning_read.txt，这是对每一个repeat的 spanning length 和read number 进行检查

                # 处理完一个trimmed seq后，将其添加到记录文件中，以备 -resume 参数时，让程序运行
                with open(record_file, "a") as f:
                    f.write(os.path.basename(fasta_file) + "\n")
            
            print()
            logging.info("Process recombination-mediated repeats in the subconfiguration! Waiting ... ")
            subcon_folders = glob.glob(f"{extrsubcon_output_dir_prefix}_{project_id}/*_results")
            for subcon_folder in subcon_folders:
                spanning_read_file = os.path.join(subcon_folder, "repeat-spanning_read.txt")    # 对每一个repeat的 repeat-spanning_read.txt 进行汇总
                subprocess.run(["python", os.path.join(dir_path, "count_subconfig_spanning_reads.py"), "-i", spanning_read_file, "-o", file_path_sub, "-cl", str(check_spanning_length), "-rn", str(read_number)], check=True, text=True, capture_output=True) 
            logging.info("Finished!")
            time.sleep(2)

################################################################################
        # 如果要求删除中间冗余结果，则可以根据剩下的 subconfiguration 提前删除一些 mainconfiguration，减少一部分比对运算
        # 预设文件名  
        output_filename = f"repeat_supp_subconfig_recomb_{project_id}.txt"  
        # 如果文件不存在，则创建文件并写入表头  
        if not os.path.isfile(output_filename):  
            with open(output_filename, 'w') as f:  
                f.write("fragment_id\tpaired_id\n")  # 写入表头，使用制表符分隔列  
        
        if mode in ['A', 'C'] and redundant_intermediate_results == "D":  
            # 初始化集合用于存储唯一行
            unique_lines = set()
            
            # 获取所有结果文件夹
            trimmed_seq_supp_recomb_folders = glob.glob(os.path.join(f"subconfig_repeat-spanned_results_{project_id}", '*_results'))  
      
            # 遍历每个文件夹，并提取相应的字段  
            for trimmed_seq_supp_recomb_folder in trimmed_seq_supp_recomb_folders:  
                folder_name = os.path.basename(trimmed_seq_supp_recomb_folder)  # 获取文件夹的名称  
                parts = folder_name.split('_')  # 使用 '_' 分割文件夹名称  
                if len(parts) >= 4:  
                    third_column = parts[2]  
                    fourth_column = parts[3]  
                    # 将提取的信息添加到集合中（自动去重）
                    unique_lines.add(f"{third_column}\t{fourth_column}\n")
            
            # 读取现有文件内容并添加到集合中（如果文件已存在）
            if os.path.exists(output_filename):
                with open(output_filename, 'r') as f:
                    unique_lines.update(f.readlines())

            # 将去重后的内容写回到文件中
            with open(output_filename, 'w') as f:
                f.writelines(unique_lines)

################################################################################
        if mode in ['A', 'C']:
            # 初始化colorama  
            init(autoreset=True)  

            maincon_fasta_files = glob.glob(f"{extrmaincon_output_dir_prefix}_{project_id}/*.fasta")
            maincon_fasta_total = len(maincon_fasta_files)    # 次要构型中需要检测的重复单元个数
            
            #### subconfiguration - Process each FASTA file sequentially
            logging.info(f"The {subcon_fasta_total} repeat units that may mediate genomic recombination have been processed completely!\n")
            time.sleep(2)
            
            print()
            logging.info(f"%%%%%%%%%%%%% Start to process {maincon_fasta_total} repeat units in the mainconfiguration! %%%%%%%%%%%%%")
            logging.info(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            time.sleep(2)
            
            file_path_main = os.path.join(extrmaincon_output_dir_prefix + "_" + project_id, "filter_mainconfig_spanned_results.tsv")

            output_final_main = os.path.dirname(file_path_main)    # 创建存放结果的文件夹
            # Ensure the directory for the results exists
            if not os.path.isdir(output_final_main):
                os.makedirs(output_final_main, exist_ok=True)

            header = "main_trimmed_seq\tmain_trimmed_seq_length(bp)\tmain_spanning_read_number"
            # Check if the file exists or is empty, and if so, create it and write the header
            if not os.path.isfile(file_path_main) or os.path.getsize(file_path_main) == 0:
                with open(file_path_main, 'w') as file:
                    file.write(header + '\n')

            # 逐一处理截取的序列
            record_file = f"record_for_resume_{project_id}.txt"    # 记录处理过的序列
            # 检查记录文件是否存在，如果存在则读取其中的内容
            if os.path.isfile(record_file):
                with open(record_file, "r") as f:
                    processed_files = f.read().splitlines()
            else:
                processed_files = []

            # 定义一个函数来检查third_column和fourth_column是否已存在于文件中  
            def check_columns_in_file(file_path, column1, column2):  
                with open(file_path, 'r') as file:  
                    for line in file:  
                        if column1 in line and column2 in line:  
                            return  False
                    return  True      # 不在记录文件中，则返回True

            output_filename = f"repeat_supp_subconfig_recomb_{project_id}.txt"
            every_maincon_count = 1      # 初始化计数器
            for fasta_file in maincon_fasta_files:

                ordinal_suffix = get_ordinal_suffix(every_maincon_count) 
                ordinal_suffix_basename = os.path.basename(ordinal_suffix)

                highlighted_text = f"{Fore.RED}{every_maincon_count}{ordinal_suffix_basename}{Style.RESET_ALL}" 
                print()
                logging.info(f"**** Start to process the {highlighted_text} sequence: {os.path.basename(fasta_file)}. ****")  
                time.sleep(5)
                every_maincon_count += 1      # 计数器递增

                if os.path.basename(fasta_file) in processed_files and args.resume:        # resume模式下，检测序列之前是否被处理过了，若果被记录处理过了，则直接跳过
                    logging.info(f"Skipping already processed sequence: {fasta_file}")
                    continue
                logging.info(f"Process trimmed seq {os.path.basename(fasta_file)} from mainconfig!!")

                if redundant_intermediate_results == "D":    # 根据对 subconfiguration 的结果解读，将不支持 subconfiguration 的重复序列提前去除，不再做mapping
                    parts = os.path.basename(fasta_file).split('_')  # 使用 '_' 分割文件名称  
                    if len(parts) >= 4:  
                        third_column = parts[2]  
                        fourth_column = parts[3]  
                    # 检查third_column和fourth_column是否已存在于文件中  
                    if check_columns_in_file(output_filename, third_column, fourth_column):  
                        logging.info(f"Skipping this trimmed seq. No reads support recombination in the subconfiguration!")
                        continue  # 如果不在记录文件中，则无需做mapping，无需做任何动作，继续下一个循环 

                extrmaincon_trimmed_length_reset = extract_trimmed_length_from_filename(fasta_file)
                logging.info(f"extrmaincon_trimmed_length_reset extrmaincon_trimmed_length_reset extrmaincon_trimmed_length_reset: {extrmaincon_trimmed_length_reset}")
                if extrmaincon_trimmed_length == 0:
                    logging.warning(f"Warning: Could not determine flanked length from file name {fasta_file}. Skipping.")
                    continue
                else:
                    extrmaincon_trimmed_length = extrmaincon_trimmed_length_reset
                
                if extrmaincon_trimmed_length < check_spanning_length:
                    logging.info(f"The reset 'flanked_sequence_length' is shorter than 'spanning_read_flanking_repeat_length' in mainconfiguration. The 'spanning_read_flanking_repeat_length' is also reset to {extrmaincon_trimmed_length}.")
                    check_spanning_length = extrmaincon_trimmed_length

                # Run seqdepth.py for each file
                map_output, _ = os.path.splitext(fasta_file)    ######################### 此处开始处理 mainconfiguration 
                if seqdepth_alignment_software != 'blast':
                    if provided_data_type == 'NGS_single_end':
                        if filter_reads in ['Y','YES']:        ############################################## 采用过滤后的reads，提高后面比对的速度。
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-single", seq_data_value_filter, "-output", map_output, "-threads", str(seqdepth_threads)] 
                        else:
                            # 不进行reads过滤，seq_data_value可以直接使用，不用区分 fasta/fastq
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-single", seq_data_value, "-output", map_output, "-threads", str(seqdepth_threads)] 
                    elif provided_data_type == 'NGS_pair_ends':
                        if filter_reads in ['Y','YES']:        ############################################## 采用过滤后的reads，提高后面比对的速度。
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-pair", fasta_pair1_filter, fasta_pair2_filter, "-output", map_output, "-threads", str(seqdepth_threads)]
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-pair", fastq1, fastq2, "-output", map_output, "-threads", str(seqdepth_threads)]
                    elif provided_data_type == 'TGS' and seqdepth_type == 'ont':
                        if filter_reads in ['Y','YES']:        ############################################## 采用过滤后的reads，提高后面比对的速度。
                            # 使用过滤后的reads seq_data_value，
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", seq_data_value_filter, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", "ont"] 
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", seq_data_value, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", "ont"]
                    elif provided_data_type == 'TGS' and seqdepth_type == 'pacbio':
                        if filter_reads in ['Y','YES']:        ############################################## 采用过滤后的reads，提高后面比对的速度。
                            # 使用过滤后的reads seq_data_value，
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", seq_data_value_filter, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", "pacbio"] 
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", seq_data_value, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", "pacbio"]
                
                ###########
                ########### 采用 blast 软件 计算测序深度
                if seqdepth_alignment_software == 'blast':
                    if provided_data_type == 'NGS_single_end':
                        if filter_reads in ['Y','YES']:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-single", seq_data_value_filter, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", seqdepth_type, "-spanning_length", str(check_spanning_length)] 
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-single", fasta_one, "-output", map_output, "-threads", str(seqdepth_threads), "-spanning_length", str(check_spanning_length), "-seqdepth_type", seqdepth_type] 

                    if provided_data_type == 'TGS':
                        if filter_reads in ['Y','YES']:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", seq_data_value_filter, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", seqdepth_type, "-spanning_length", str(check_spanning_length)] 
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-third", fasta_one, "-output", map_output, "-threads", str(seqdepth_threads), "-seqdepth_type", seqdepth_type, "-spanning_length", str(check_spanning_length)] 
                            
                    if provided_data_type == 'NGS_pair_ends':
                        if filter_reads in ['Y','YES']:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-pair", fasta_pair1_filter, fasta_pair2_filter, "-output", map_output, "-threads", str(seqdepth_threads), "-spanning_length", str(check_spanning_length), "-seqdepth_type", seqdepth_type]
                        else:
                            command = ["python", os.path.join(dir_path, "seqdepth.py"), "-alignment", seqdepth_alignment_software, "-reference", fasta_file, "-pair", fasta_pair1, fasta_pair2, "-output", map_output, "-threads", str(seqdepth_threads), "-spanning_length", str(check_spanning_length), "-seqdepth_type", seqdepth_type]

                run_command(command)

                # Immediately after, run check_mainconfig_spanning_reads.py for the output
                map_results_folder = f"{map_output}_results"
                if os.path.isdir(map_results_folder):
                    output_sam = os.path.join(map_results_folder, "output.sam")
                    run_command(["python", os.path.join(dir_path, "check_mainconfig_spanning_reads.py"), "-s", output_sam, "-i", fasta_file, "-o", map_results_folder, "-cl", str(check_spanning_length), "-tl", str(extrmaincon_trimmed_length), "-a", seqdepth_alignment_software])
                
                # 处理完一个文件后，将其添加到记录文件中
                with open(record_file, "a") as f:
                    f.write(os.path.basename(fasta_file) + "\n")
            
            print()
            logging.info("Process repeat information in the mainconfiguration! Waiting ... ")
            maincon_folders = glob.glob(f"{extrmaincon_output_dir_prefix}_{project_id}/*_results")
            for maincon_folder in maincon_folders:
                spanning_read_file = os.path.join(maincon_folder, "repeat-spanning_read.txt")    # 对每一个repeat的 repeat-spanning_read.txt 进行汇总
                subprocess.run(["python", os.path.join(dir_path, "count_mainconfig_spanning_reads.py"), "-i", spanning_read_file, "-o", file_path_main, "-cl", str(check_spanning_length)], check=True, text=True, capture_output=True)
            logging.info("Finished!")
            time.sleep(2)
                    
            ############################
            # mapping 之后，删除 由fastq转化来的 中间结果 fasta文件，或者是过滤后的reads，直接删掉 
            if fastq2fasta and provided_data_type == 'NGS_single_end':    
                os.remove(f"{fasta_one}") if os.path.isfile(f"{fasta_one}") else None
            if fastq2fasta and provided_data_type == 'TGS': 
                os.remove(f"{fasta_one}") if os.path.isfile(f"{fasta_one}") else None
            if fastq2fasta and provided_data_type == 'NGS_pair_ends':
                os.remove(f"{fasta_pair1}") if os.path.isfile(f"{fasta_pair1}") else None
                os.remove(f"{fasta_pair2}") if os.path.isfile(f"{fasta_pair2}") else None
                
            #对于过滤出来的reads，保存下来。
            # 确保目标目录存在
            if not os.path.exists(project_id):
                os.makedirs(project_id, exist_ok=True)
                
            # 获取当前目录下所有以 _filtered.fasta 结尾的文件
            filtered_files = [f for f in os.listdir('.') if f.endswith(f'_{project_id}_filtered.fasta')]
            # 遍历文件列表并移动文件
            for file_name in filtered_files:
                # 构造源文件路径和目标文件路径
                source_path = os.path.join('.', file_name)  # 当前目录下的文件路径
                destination_path = os.path.join(project_id, file_name)  # 目标目录下的文件路径
                # 移动文件
                shutil.move(source_path, destination_path)

            #### mainconfiguration - Process each FASTA file sequentially
            logging.info(f"The {maincon_fasta_total} repeat units in the mainconfiguration have been processed completely! \n")
            time.sleep(2)
            
            if redundant_intermediate_results == "D":    # 将不支持 subconfiguration 的 repeat 对应的 mainconfiguration 的 fasta 序列删除
                maincon_fasta_files = glob.glob(f"{extrmaincon_output_dir_prefix}_{project_id}/*.fasta")
                for fasta_file in maincon_fasta_files:
                    os.remove(fasta_file) if os.path.isfile(fasta_file) else None

################################################################################
        if mode == 'A': 
            # 汇总结果，用于绘图
            # subconfiguration
            run_command(["python", os.path.join(dir_path, "filter_spanning_results.py"), "-sr", extrsubcon_output_dir_prefix + "_" + project_id + "/filter_subconfig_spanned_results.tsv", "-mr", extrmaincon_output_dir_prefix + "_" + project_id + "/filter_mainconfig_spanned_results.tsv", "-rf", f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_5CT.tsv", "-rc", f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_calibration_table.txt", "-o", output_dir_prefix + "_" + project_id + "/", "-rn", str(read_number), "-cr", complementary_chain])

        if mode == 'C': 
            run_command(["python", os.path.join(dir_path, "filter_spanning_results.py"), "-sr", extrsubcon_output_dir_prefix + "_" + project_id + "/filter_subconfig_spanned_results.tsv", "-mr", extrmaincon_output_dir_prefix + "_" + project_id + "/filter_mainconfig_spanned_results.tsv", "-rf", f"ROUSFinder_results_{project_id}/" + prefix_manual + "_" + project_id + "_rep_5CT.tsv", "-rc", manually_calibrate_reset, "-o", output_dir_prefix + "_" + project_id + "/", "-rn", str(read_number), "-cr", complementary_chain])

        if mode in ['A', 'C']:
            run_command(["python", os.path.join(dir_path, "filter_correspondence.py"), "-mp", extrmaincon_output_dir_prefix + "_" + project_id, "-sr", extrsubcon_output_dir_prefix + "_" + project_id + "/filter_subconfig_spanned_results.tsv", "-mr", extrmaincon_output_dir_prefix + "_" + project_id + "/filter_mainconfig_spanned_results.tsv", "-rn", str(read_number), "-o", output_dir_prefix + "_" + project_id, "-ml", redundant_intermediate_results, "-t", str(seqdepth_threads)])

            logging.info(f"Please check the final confirmed paired repeats in the file 'paired_repeats_for_mapping.tsv'. They can be used to map mitogenomic recombination.")
            logging.info(f"Please check the recombination ratio in the file 'paired_repeats_recomb-supporting_ratio.tsv'. They are the ratio of mitogenomic recombination.")
            logging.info("Some repeat units from mainconfiguration (mcfg) may be with 'insufficient' spanning reads. If not keep, delete it/them before mapping.")
            logging.info("If only header row is displayed in the TSV file, it is because the program don't find paired repeat that can mediate genome recombination.")

########################################################################################################################################################################
        # 对产生的支持基因组重组的repeat的结果进行start和end的修正，也就是找出每个重复在每一条染色体的真实位置
        # 重新归档产生的结果
###############################################################################
        # 以下处理 ROUSFinder 查找的 repeat 结果，追溯 start 和 end 的具体位置，重新归档，位置在 ROUSFinder_results_project_id 文件夹中
        if mode in ['A', 'C', 'R'] and count == 1:
            pass
        if mode == 'A' and count > 1:
            if os.path.isfile(f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt"):  
                file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt"           # rous_output_file_prefix 为读入的fasta序列文件的名字
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path):   
                    default_headers = ['fragment_id','length','start','end','direction']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                df = pd.read_csv(file_path, sep='\t')  
                df = df.dropna(how='all')     # 去掉空白行
                
                # 检查df是否数据为空，排除表头  
                if len(df) == 0:  
                    logging.warning(f"The DataFrame {file_path} is empty!!")  
                    # 如果DataFrame为空，则复制原始文件并添加"pseudo_"前缀  
                    new_file_path = f"ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_table.txt"  
                    shutil.copy(file_path, new_file_path)  
                else: 
                    os.rename(f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt", f"ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_table.txt" ) 
                    file_path = f"ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_table.txt"  
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')     # 去掉空白行
                    # 指定新的表头  
                    new_header = ['fragment_id', 'length', 'start', 'end', 'direction', 'chromosome'] 
                    df = pd.read_csv(file_path, sep='\t', names=new_header, header=None, skiprows=1)    
                    # 以上从file_path指定的路径读取一个以制表符分隔的文件，跳过文件的第一行（可能是原始的表头或描述信息），不使用文件中的任何行作为列名，而是使用new_header列表中指定的名称作为DataFrame的列名。
                    # 读取的数据将被存储在一个新的DataFrame对象df中。
                
                    # 应用函数并分别更新列  
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = np.where(df['adjusted_start'] > 0, df['adjusted_start'], 1)
                    df['end'] = np.where(df['adjusted_end'] > 0, df['adjusted_end'], 1)
                    
                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])     # delete the two columns: 'adjusted_start', 'adjusted_end'
                    # 将chromosome列插入到direction列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('direction') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  

                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]
                    # 保存修改后的 DataFrame  
                    file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_table.txt"   
                    df.to_csv(file_path, sep='\t', index=False)  

            if os.path.isfile(f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt"):  
                file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt" 
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path):         
                    default_headers = ['fragment_id','length','start','end','direction']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                df = pd.read_csv(file_path, sep='\t')  
                df = df.dropna(how='all')     # 去掉空白行

                # 检查df数据是否为空，排除表头  
                if len(df) == 0:  
                    logging.warning(f"The DataFrame {file_path} is empty!!") 
                    # 如果DataFrame为空，则复制原始文件并添加"pseudo_"前缀  
                    new_file_path = f"ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt"  
                    shutil.copy(file_path, new_file_path)  
                else: 
                    os.rename(f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt", f"ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt" ) 
                    file_path = f"ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt"  
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')     # 去掉空白行

                    # 应用函数并分别更新列  
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = np.where(df['adjusted_start'] > 0, df['adjusted_start'], 1)
                    df['end'] = np.where(df['adjusted_end'] > 0, df['adjusted_end'], 1)

                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])     #删除两列
                    # 将chromosome列插入到 direction 列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('direction') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  

                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]
                    # 保存修改后的 DataFrame  
                    file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt"  
                    df.to_csv(file_path, sep='\t', index=False)  

            if os.path.isfile(f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_5CT.tsv"):  
                file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_5CT.tsv"  
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path):             #os.path.getsize(file_path) == 0: 
                    default_headers = ['fragment_id','start','end','type','paired_id']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                df = pd.read_csv(file_path, sep='\t')  
                df = df.dropna(how='all')     # 去掉空白行

                # 检查df数据是否为空，排除表头  
                if len(df) == 0:  
                    logging.warning(f"The DataFrame {file_path} is empty!!") 
                    # 如果DataFrame为空，则复制原始文件并添加"pseudo_"前缀  
                    new_file_path = f"ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_5CT.tsv"  
                    shutil.copy(file_path, new_file_path)  
                else: 
                    os.rename(f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_5CT.tsv", f"ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_5CT.tsv" ) 
                    file_path = f"ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_5CT.tsv"  
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')     # 去掉空白行

                    # 应用函数并分别更新列  
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = np.where(df['adjusted_start'] > 0, df['adjusted_start'], 1)
                    df['end'] = np.where(df['adjusted_end'] > 0, df['adjusted_end'], 1)

                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])     # 删除两列
                    # 将chromosome列插入到 paired_id 列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('paired_id') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  

                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]
                    # 保存修改后的 DataFrame  
                    file_path = f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_5CT.tsv" 
                    df.to_csv(file_path, sep='\t', index=False)  

            if os.path.isfile(f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep.fasta"):
                os.rename(f"ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep.fasta", f"ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep.fasta")                             # 对合并后的序列文件改名，标记为pseudo_

###############################################################################
        # 处理 mode == "C"，即用户提供 repeat 手工修正后的结果
        if mode == "C" and count >1:
            if os.path.isfile(f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_5CT.tsv"):
                file_path = f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_5CT.tsv"
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path):   
                    default_headers = ['fragment_id','start','end','type','paired_id']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                df = pd.read_csv(file_path, sep='\t')
                df = df.dropna(how='all')    # 去掉空白行

                if len(df) == 0:
                    logging.warning(f"The DataFrame {file_path} is empty!!")
                    new_file_path = f"ROUSFinder_results_{project_id}/pseudo_{prefix_manual}_{project_id}_rep_5CT.txt"  
                    shutil.copy(file_path, new_file_path)  
                else:
                    os.rename(f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_5CT.tsv", f"ROUSFinder_results_{project_id}/pseudo_{prefix_manual}_{project_id}_rep_5CT.tsv") 
                    file_path = f"ROUSFinder_results_{project_id}/pseudo_{prefix_manual}_{project_id}_rep_5CT.tsv"
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')     # 去掉空白行

                    # 应用函数并分别更新列  
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = np.where(df['adjusted_start'] > 0, df['adjusted_start'], 1)
                    df['end'] = np.where(df['adjusted_end'] > 0, df['adjusted_end'], 1)
                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])    # 删除两列
                    # 将chromosome列插入到 paired_id 列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('paired_id') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  

                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]
                    # 保存修改后的 DataFrame  
                    file_path = f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_5CT.tsv" 
                    df.to_csv(file_path, sep='\t', index=False)              
            if os.path.isfile(f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_calibration_table.txt"):  
                file_path = f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_calibration_table.txt" 
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path):         
                    default_headers = ['fragment_id','length','start','end','direction']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                df = pd.read_csv(file_path, sep='\t')  
                df = df.dropna(how='all')     # 去掉空白行

                # 检查df数据是否为空，排除表头  
                if len(df) == 0:  
                    logging.warning(f"The DataFrame {file_path} is empty!!") 
                    # 如果DataFrame为空，则复制原始文件并添加"pseudo_"前缀  
                    new_file_path = f"ROUSFinder_results_{project_id}/pseudo_{prefix_manual}_{project_id}_rep_calibration_table.txt"  
                    shutil.copy(file_path, new_file_path)  
                else: 
                    os.rename(f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_calibration_table.txt", f"ROUSFinder_results_{project_id}/pseudo_{prefix_manual}_{project_id}_rep_calibration_table.txt" ) 
                    file_path = f"ROUSFinder_results_{project_id}/pseudo_{prefix_manual}_{project_id}_rep_calibration_table.txt"  
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')     # 去掉空白行

                    # 应用函数并分别更新列  
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = np.where(df['adjusted_start'] > 0, df['adjusted_start'], 1)
                    df['end'] = np.where(df['adjusted_end'] > 0, df['adjusted_end'], 1)

                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])     #删除两列
                    # 将chromosome列插入到 direction 列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('direction') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  

                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]
                    # 保存修改后的 DataFrame  
                    file_path = f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_calibration_table.txt"  
                    df.to_csv(file_path, sep='\t', index=False)  

################################################################################
        # 以下处理 repeat 查询结果的重新归档与start和end具体位置的追溯，位置在 final_repeat-spanning_results_{project_id}文件夹中
        if mode in ['A', 'C', 'R']:    # 'R' 考虑仅仅运行ROUSFinfer时，将其结果重新归档
            if not os.path.isdir(f"{project_id}"):
                os.mkdir(f"{project_id}")
            shutil.move(f"{output_dir_prefix}_{project_id}", f"{project_id}") if os.path.isdir(f"{output_dir_prefix}_{project_id}") else None
            shutil.move(f"{extrsubcon_output_dir_prefix}_{project_id}", f"{project_id}") if os.path.isdir(f"{extrsubcon_output_dir_prefix}_{project_id}") else None
            shutil.move(f"{extrmaincon_output_dir_prefix}_{project_id}", f"{project_id}") if os.path.isdir(f"{extrmaincon_output_dir_prefix}_{project_id}") else None
            shutil.move(f"custom_log_{project_id}.txt", f"{project_id}") if os.path.isfile(f"custom_log_{project_id}.txt") else None
            shutil.move(f"ROUSFinder_results_{project_id}", f"{project_id}") if os.path.isdir(f"ROUSFinder_results_{project_id}") else None

        if mode in ['A', 'C'] and count > 1:
            if os.path.isfile(f"{project_id}/{output_dir_prefix}_{project_id}/one_chain_without_sufficient_spanning_reads.tsv"):  
                file_path = f"{project_id}/{output_dir_prefix}_{project_id}/one_chain_without_sufficient_spanning_reads.tsv"  
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path): 
                    default_headers = ['fragment_id', 'length', 'start', 'end', 'direction', 'chromosome', 'plus_ratio(s/m)', 'minus_ratio(s/m)', 'combined_ratio', 'type', 'paired_id', 'paired_length', 'paired_start', 'paired_end', 'paired_direction', 'paired_chromosome', 'paired_plus_ratio(s/m)', 'paired_minus_ratio(s/m)', 'paired_combined_ratio', 'spanning_read_mcfg']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 

                df = pd.read_csv(file_path, sep='\t') 
                df = df.dropna(how='all')  # 去掉空白行

                if len(df) == 0:
                    logging.warning(f"The DataFrame {file_path} is empty!!")
                    shutil.copy(f"{project_id}/{output_dir_prefix}_{project_id}/one_chain_without_sufficient_spanning_reads.tsv",   
                                f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_one_chain_without_sufficient_spanning_reads.tsv")  
                else:
                    os.rename(f"{project_id}/{output_dir_prefix}_{project_id}/one_chain_without_sufficient_spanning_reads.tsv",   
                              f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_one_chain_without_sufficient_spanning_reads.tsv")  
                    file_path = f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_one_chain_without_sufficient_spanning_reads.tsv"  
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')  # 去掉空白行

                    # 应用函数并分别更新列  
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = np.where(df['adjusted_start'] > 0, df['adjusted_start'], 1)
                    df['end'] = np.where(df['adjusted_end'] > 0, df['adjusted_end'], 1)
                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])     # 删除两列

                    # 将chromosome列插入到direction列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('direction') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  

                    # 处理 paired_start 和 paired_end  
                    df[['paired_chromosome', 'adjusted_paired_start', 'adjusted_paired_end']] = df.apply(lambda row: find_chromosome(row['paired_start'], row['paired_end'], cumulative_length), axis=1, result_type='expand')  
                    df['paired_start'] = np.where(df['adjusted_paired_start'] > 0, df['adjusted_paired_start'], 1)
                    df['paired_end'] = np.where(df['adjusted_paired_end'] > 0, df['adjusted_paired_end'], 1)
                    df = df.drop(columns=['adjusted_paired_start', 'adjusted_paired_end'])     # 删除这两列

                    # 将paired_chromosome列插入到paired_direction列之后  
                    paired_chromosome_index = cols.index('paired_direction') + 1  
                    cols.insert(paired_chromosome_index, 'paired_chromosome')  
                    df = df[cols]  

                    # 删除最后1列
                    if len(df.columns) > 1:
                        df = df.iloc[:, :-1]
                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]
            
                    # 保存修改后的 DataFrame  
                    file_path = f"{project_id}/{output_dir_prefix}_{project_id}/one_chain_without_sufficient_spanning_reads.tsv"  
                    df.to_csv(file_path, sep='\t', index=False)

            if os.path.isfile(f"{project_id}/{output_dir_prefix}_{project_id}/one_repeat_unit_without_spanning_reads.tsv"):
                file_path = f"{project_id}/{output_dir_prefix}_{project_id}/one_repeat_unit_without_spanning_reads.tsv"
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path):    
                    default_headers = ['fragment_id', 'length', 'start', 'end', 'direction', 'chromosome', 'plus_ratio(s/m)', 'minus_ratio(s/m)', 'combined_ratio', 'type', 'paired_id', 'paired_length', 'paired_start', 'paired_end', 'paired_direction', 'paired_chromosome', 'paired_plus_ratio(s/m)', 'paired_minus_ratio(s/m)', 'paired_combined_ratio', 'spanning_read_mcfg']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                df = pd.read_csv(file_path, sep='\t')
                df = df.dropna(how='all')  # 去掉空白行
                
                if len(df) == 0:
                    logging.warning(f"The DataFrame {file_path} is empty!!")
                    shutil.copy(f"{project_id}/{output_dir_prefix}_{project_id}/one_repeat_unit_without_spanning_reads.tsv",
                                f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_one_repeat_unit_without_spanning_reads.tsv")
                else:
                    os.rename(f"{project_id}/{output_dir_prefix}_{project_id}/one_repeat_unit_without_spanning_reads.tsv",
                              f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_one_repeat_unit_without_spanning_reads.tsv")
                    file_path = f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_one_repeat_unit_without_spanning_reads.tsv"
                    df = pd.read_csv(file_path, sep='\t')
                    df = df.dropna(how='all')  # 去掉空白行
            
                    # 应用函数并分别更新列
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')
                    df['start'] = np.where(df['adjusted_start'] > 0, df['adjusted_start'], 1)
                    df['end'] = np.where(df['adjusted_end'] > 0, df['adjusted_end'], 1)
                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])     # 删除这两列
            
                    # 将chromosome列插入到direction列之后
                    cols = df.columns.tolist()
                    chromosome_index = cols.index('direction') + 1
                    cols.insert(chromosome_index, 'chromosome')
                    df = df[cols]
            
                    # 处理 paired_start 和 paired_end
                    df[['paired_chromosome', 'adjusted_paired_start', 'adjusted_paired_end']] = df.apply(lambda row: find_chromosome(row['paired_start'], row['paired_end'], cumulative_length), axis=1, result_type='expand')
                    df['paired_start'] = np.where(df['adjusted_paired_start'] > 0, df['adjusted_paired_start'], 1)
                    df['paired_end'] = np.where(df['adjusted_paired_end'] > 0, df['adjusted_paired_end'], 1)
                    df = df.drop(columns=['adjusted_paired_start', 'adjusted_paired_end'])     # 删除两列
            
                    # 将paired_chromosome列插入到paired_direction列之后
                    paired_chromosome_index = cols.index('paired_direction') + 1
                    cols.insert(paired_chromosome_index, 'paired_chromosome')
                    df = df[cols]

                    # 删除最后1列
                    if len(df.columns) > 1:
                        df = df.iloc[:, :-1]
                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]
            
                    # 保存修改后的 DataFrame
                    file_path = f"{project_id}/{output_dir_prefix}_{project_id}/one_repeat_unit_without_spanning_reads.tsv"
                    df.to_csv(file_path, sep='\t', index=False)

            if os.path.isfile(f"{project_id}/{output_dir_prefix}_{project_id}/paired_repeats_recomb-supporting_ratio.tsv"):
                file_path = f"{project_id}/{output_dir_prefix}_{project_id}/paired_repeats_recomb-supporting_ratio.tsv"
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path):    
                    default_headers = ['fragment_id', 'length', 'start', 'end', 'direction', 'chromosome', 'plus_ratio(s/m)', 'minus_ratio(s/m)', 'combined_ratio', 'type', 'paired_id', 'paired_length', 'paired_start', 'paired_end', 'paired_direction', 'paired_chromosome', 'paired_plus_ratio(s/m)', 'paired_minus_ratio(s/m)', 'paired_combined_ratio', 'spanning_read_mcfg']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                df = pd.read_csv(file_path, sep='\t') 
                df = df.dropna(how='all')     # 去掉空白行
                
                if len(df) == 0:
                    logging.warning(f"The DataFrame {file_path} is empty!!")
                    shutil.copy(f"{project_id}/{output_dir_prefix}_{project_id}/paired_repeats_recomb-supporting_ratio.tsv",   
                                f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_paired_repeats_recomb-supporting_ratio.tsv")  
                else:
                    os.rename(f"{project_id}/{output_dir_prefix}_{project_id}/paired_repeats_recomb-supporting_ratio.tsv",   
                              f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_paired_repeats_recomb-supporting_ratio.tsv")  
                    file_path = f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_paired_repeats_recomb-supporting_ratio.tsv"  
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')     # 去掉空白行
            
                    # 应用函数并分别更新列  
                    temp_result = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = np.where(df['adjusted_start'] > 0, df['adjusted_start'], 1)
                    df['end'] = np.where(df['adjusted_end'] > 0, df['adjusted_end'], 1)
                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])      #删除两列
            
                    # 将chromosome列插入到direction列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('direction') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  
                    # 删除重复的chromosome列
                    df = df.loc[:,~df.columns.duplicated()]
            
                    # 处理 paired_start 和 paired_end  
                    df[['paired_chromosome', 'adjusted_paired_start', 'adjusted_paired_end']] = df.apply(lambda row: find_chromosome(row['paired_start'], row['paired_end'], cumulative_length), axis=1, result_type='expand')  
                    df['paired_start'] = np.where(df['adjusted_paired_start'] > 0, df['adjusted_paired_start'], 1)
                    df['paired_end'] = np.where(df['adjusted_paired_end'] > 0, df['adjusted_paired_end'], 1)
                    df = df.drop(columns=['adjusted_paired_start', 'adjusted_paired_end'])  
            
                    # 将paired_chromosome列插入到paired_direction列之后  
                    cols = df.columns.tolist()  
                    paired_chromosome_index = cols.index('paired_direction') + 1  
                    cols.insert(paired_chromosome_index, 'paired_chromosome')  
                    df = df[cols]  

                    # 删除多余的最后一列
                    if len(df.columns) > 1:
                        df = df.iloc[:, :-1]
            
                    # 保存修改后的 DataFrame  
                    file_path = f"{project_id}/{output_dir_prefix}_{project_id}/paired_repeats_recomb-supporting_ratio.tsv"  
                    df.to_csv(file_path, sep='\t', index=False)  

            if os.path.isfile(f"{project_id}/{output_dir_prefix}_{project_id}/paired_repeats_for_mapping.tsv"):
                file_path = f"{project_id}/{output_dir_prefix}_{project_id}/paired_repeats_for_mapping.tsv"
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path):  
                    default_headers = ['fragment_id', 'start', 'end', 'direction', 'chromosome', 'paired_id', 'paired_start', 'paired_end', 'paired_direction', 'paired_chromosome', 'spanning_read_mcfg']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                df = pd.read_csv(file_path, sep='\t') 
                df = df.dropna(how='all')  # 去掉空白行
    
                if len(df) == 0:
                    logging.warning(f"The DataFrame {file_path} is empty!!")
                    shutil.copy(f"{project_id}/{output_dir_prefix}_{project_id}/paired_repeats_for_mapping.tsv",   
                    f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_paired_repeats_for_mapping.tsv")  
                else:
                    os.rename(f"{project_id}/{output_dir_prefix}_{project_id}/paired_repeats_for_mapping.tsv",   
                              f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_paired_repeats_for_mapping.tsv")  
                    file_path = f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_paired_repeats_for_mapping.tsv"  
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')  # 去掉空白行

                    # 应用函数并分别更新列  
                    temp_result = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
        
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = np.where(df['adjusted_start'] > 0, df['adjusted_start'], 1)
                    df['end'] = np.where(df['adjusted_end'] > 0, df['adjusted_end'], 1)
                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])  

                    # 将chromosome列插入到direction列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('direction') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  
        
                    # 更新cols以避免重复列  
                    cols = df.columns.tolist()  

                    # 处理 paired_start 和 paired_end  
                    df[['paired_chromosome', 'adjusted_paired_start', 'adjusted_paired_end']] = df.apply(lambda row: find_chromosome(row['paired_start'], row['paired_end'], cumulative_length), axis=1, result_type='expand')  
                    df['paired_start'] = np.where(df['adjusted_paired_start'] > 0, df['adjusted_paired_start'], 1)
                    df['paired_end'] = np.where(df['adjusted_paired_end'] > 0, df['adjusted_paired_end'], 1)
                    df = df.drop(columns=['adjusted_paired_start', 'adjusted_paired_end'])     # 删除这两列

                    # 将paired_chromosome列插入到paired_direction列之后  
                    paired_chromosome_index = cols.index('paired_direction') + 1  
                    cols.insert(paired_chromosome_index, 'paired_chromosome')  
                    df = df[cols]  

                    # 删除最后1列
                    if len(df.columns) > 1:
                        df = df.iloc[:, :-1]
                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]
        
                    # 保存修改后的 DataFrame  
                    file_path = f"{project_id}/{output_dir_prefix}_{project_id}/paired_repeats_for_mapping.tsv"
                    df.to_csv(file_path, sep='\t', index=False)  

################################################################################
        if mode in ['A', 'C']: 
            # 删除形成的fasta文件，释放空间，这是为blast运行产生的fasta格式 
            if 'fasta_prefix' in locals() and os.path.isfile(f"{fasta_prefix}.fasta") and redundant_intermediate_results == "D":
                os.remove(f"{fasta_prefix}.fasta")
            if 'fasta1_prefix' in locals() and os.path.isfile(f"{fasta1_prefix}.fasta") and redundant_intermediate_results == "D":
                os.remove(f"{fasta1_prefix}.fasta")
            if 'fasta2_prefix' in locals() and os.path.isfile(f"{fasta2_prefix}.fasta") and redundant_intermediate_results == "D":
                os.remove(f"{fasta2_prefix}.fasta")
                
            # 重新归档 多条序列合并后的pseudo-genome
            shutil.move(f"cat_inputfasta_{project_id}.fasta", f"{project_id}") if os.path.isfile(f"cat_inputfasta_{project_id}.fasta") else None
            # 重新归档 支持基因组重组的重复序列单元
            shutil.move(f"repeat_supp_subconfig_recomb_{project_id}.txt", f"{project_id}") if os.path.isfile(f"repeat_supp_subconfig_recomb_{project_id}.txt") else None
            # 重新归档 断点继续计算的支持文件
            shutil.move(f"record_for_resume_{project_id}.txt", f"{project_id}") if os.path.isfile(f"record_for_resume_{project_id}.txt") else None
            
            if filter_reads in ['Y','YES']:        # 当过滤的时候，判断是否保留过滤的中间结果 
                # 重新归档筛选出来的fasta序列，由read_filter = YES 时过滤出来的reads，为fasta格式，文件名含有project_id
                if provided_data_type in ['NGS_single_end', 'TGS'] and redundant_intermediate_results == 'K':
                    shutil.move(f"{seq_data_value_filter}", f"{project_id}") if os.path.isfile(f"{seq_data_value_filter}") else None
                if fastq2fasta and provided_data_type in ['NGS_single_end', 'TGS'] and redundant_intermediate_results == 'D':
                    os.remove(f"{seq_data_value_filter}") if os.path.isfile(f"{seq_data_value_filter}") else None
                    
                # 以上筛选的fasta序列一定存在，单端和三代均能产生，但是双端数据必须是用的双端数据时才会产生
                if provided_data_type == 'NGS_pair_ends' and redundant_intermediate_results == 'K': 
                    shutil.move(f"{fasta_pair1_filter}", f"{project_id}") if os.path.isfile(f"{fasta_pair1_filter}") else None
                    shutil.move(f"{fasta_pair2_filter}", f"{project_id}") if os.path.isfile(f"{fasta_pair2_filter}") else None
                    
                if provided_data_type == 'NGS_pair_ends' and redundant_intermediate_results == 'D': 
                    os.remove(f"{fasta_pair1_filter}") if os.path.isfile(f"{fasta_pair1_filter}") else None
                    os.remove(f"{fasta_pair2_filter}") if os.path.isfile(f"{fasta_pair2_filter}") else None

            # 删除重复序列单元的对应关系  final_repeat-spanning_results
            if os.path.isfile(f"{project_id}/{output_dir_prefix}_{project_id}/recomb-supporting_paired_repeat_correspondence.tsv"):
                os.remove(f"{project_id}/{output_dir_prefix}_{project_id}/recomb-supporting_paired_repeat_correspondence.tsv") 
         
        if mode == 'R':
            os.remove(f"cat_inputfasta_{project_id}.fasta") if os.path.isfile(f"cat_inputfasta_{project_id}.fasta") else None
            os.remove(f"repeat_supp_subconfig_recomb_{project_id}.txt") if os.path.isfile(f"repeat_supp_subconfig_recomb_{project_id}.txt") else None
        
        if count > 1:      # 重新归档pseudo结果
            if not os.path.isdir(f"{project_id}/pseudo_{output_dir_prefix}_{project_id}"):
                os.mkdir(f"{project_id}/pseudo_{output_dir_prefix}_{project_id}")
                
            # 使用glob模块匹配通配符的文件
            for file_path in glob.glob(f"{project_id}/{output_dir_prefix}_{project_id}/pseudo_*"):
                # 获取目标路径
                destination = f"{project_id}/pseudo_{output_dir_prefix}_{project_id}"
                # 移动文件到目标目录
                shutil.move(file_path, destination)       ##已经实现伪结果的重新归档   
            
        if count > 1: 
            subprocess.run(["python", os.path.join(dir_path, "extract_fasta.py"), inputfasta, project_id], check=True)         # 将输入的fasta拆分为一条序列一个文件
            if mode == 'A':
                subprocess.run(["python", os.path.join(dir_path, "split_6CT_to_5CT.py"), f"{project_id}/ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_5CT.tsv"], check=True)      # 将输入的5CT文件，按照chrosome编号，一条染色体拆分为一个文件
            if mode == 'C':
                subprocess.run(["python", os.path.join(dir_path, "split_6CT_to_5CT.py"), f"{project_id}/ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_5CT.tsv"], check=True)      # 将输入的5CT文件，按照chrosome编号，一条染色体拆分为一个文件
            # 汇总重复序列在各条染色体上的位置，截取重复序列。dir_path：主程序所在的路径
            i = 1
            while i <= count:
                run_command(["python", os.path.join(dir_path, "process_repeats_from_5CT_OneByOne.py"), "-p", f"{project_id}", "-s", f"chr{i}.fasta", "-r", f"chr{i}_5CT.tsv", "-n", f"{i}"]) if os.path.isfile(f"chr{i}_5CT.tsv") else None
                
                os.remove(f"{project_id}/{output_dir_prefix}_{project_id}/repeat_positions_{project_id}_chr{i}.tsv") if os.path.isfile(f"{project_id}/{output_dir_prefix}_{project_id}/repeat_positions_{project_id}_chr{i}.tsv") else None
                os.remove(f"{project_id}/{output_dir_prefix}_{project_id}/repeat_sequences_{project_id}_chr{i}.fasta") if os.path.isfile(f"{project_id}/{output_dir_prefix}_{project_id}/repeat_sequences_{project_id}_chr{i}.fasta") else None
                shutil.move(f"repeat_positions_{project_id}_chr{i}.tsv", f"{project_id}/{output_dir_prefix}_{project_id}") if os.path.isfile(f"repeat_positions_{project_id}_chr{i}.tsv") else None
                shutil.move(f"repeat_sequences_{project_id}_chr{i}.fasta", f"{project_id}/{output_dir_prefix}_{project_id}") if os.path.isfile(f"repeat_sequences_{project_id}_chr{i}.fasta") else None
                
                os.remove(f"chr{i}.fasta") if os.path.isfile(f"chr{i}.fasta") else None          # 删除中间结果
                os.remove(f"chr{i}_5CT.tsv") if os.path.isfile(f"chr{i}_5CT.tsv") else None          # 删除中间结果
                i += 1         # 计时器
                
            if redundant_intermediate_results == "D":       # 删除中间结果 处理多条序列的 伪基因组的 中间结果
                shutil.rmtree(f"{project_id}/pseudo_{output_dir_prefix}_{project_id}") if os.path.isdir(f"{project_id}/pseudo_{output_dir_prefix}_{project_id}") else None
                
        if count == 1:
            # 汇总重复序列在基因组上的位置，截取重复序列。dir_path：主程序所在的路径
            if mode == 'A' and os.path.isfile(f"{project_id}/ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_5CT.tsv"):
                run_command(["python", os.path.join(dir_path, "process_repeats_from_5CT_OneByOne.py"), "-p", project_id, "-s", f"{inputfasta}", "-r", f"{project_id}/ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_5CT.tsv"]) 
            if mode == 'C' and os.path.isfile(f"{project_id}/ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_5CT.tsv"):
                run_command(["python", os.path.join(dir_path, "process_repeats_from_5CT_OneByOne.py"), "-p", project_id, "-s", f"{inputfasta}", "-r", f"{project_id}/ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_5CT.tsv"])
                
            os.remove(f"{project_id}/{output_dir_prefix}_{project_id}/repeat_positions_{project_id}.tsv") if os.path.isfile(f"{project_id}/{output_dir_prefix}_{project_id}/repeat_positions_{project_id}.tsv") else None
            os.remove(f"{project_id}/{output_dir_prefix}_{project_id}/repeat_sequences_{project_id}.fasta") if os.path.isfile(f"{project_id}/{output_dir_prefix}_{project_id}/repeat_sequences_{project_id}.fasta") else None
            shutil.move(f"repeat_positions_{project_id}.tsv", f"{project_id}/{output_dir_prefix}_{project_id}") if os.path.isfile(f"repeat_positions_{project_id}.tsv") else None
            shutil.move(f"repeat_sequences_{project_id}.fasta", f"{project_id}/{output_dir_prefix}_{project_id}") if os.path.isfile(f"repeat_sequences_{project_id}.fasta") else None
        
        modify_headers_in_place(inputfasta)    # 将输入的fasta中添加的chr{i}_删除，但之前将空格替换的"_"没有再替换为空格。

        logging.info("The computational process of searching for repeats that can mediate genome recombination has been completed.!\n")

################################################################################
    elif  mode in ['A', 'R', 'C'] and refilter_mode:    # 主程序的运行逻辑
        logging.error("Please reset 'refilter_mode = N' in the [refilter_params] section. Two modes can't be run at the same time!\n")
        sys.exit(1)

######################################################################################################################################################
######################################################################################################################################################
################################################################################
################################################################################
    #### 下面对前面查找的 repeat 结果，再次进行过滤。
    #### 可以独立于之前的程序而运行，但要基于之前程序的运行结果
    
    # 重新筛选的原理
    # 基于以下两个中间结果文件，重新计算repeat结果
    # redundant_intermediate_results == 'D'时，基于用户设置的 read number 和 spanning length，检索所有的结果，统计所有reads跨越repeat的情况
    # redundant_intermediate_results == 'K'时, 基于read number=1 和 spanning length=1，检索所有的结果，统计所有reads跨越repeat的情况
    # reads跨越repeat的情况全部写入了"repeat-spanning_read.txt"，重过滤的时候，只要反复读取"repeat-spanning_read.txt"即可
    # 目前可以实现 read number 和 spanning length 两个参数不断进行调整的重过滤

    if mode == "N" and refilter_mode:
        print("Current work is to RE-FILTER paired repeats that support genome recombination.\n")
        time.sleep(5)

        # 处理项目号
        project_id = config['general']['project_id']
        refilter_id = config['refilter_params']['refilter_id']
        redundant_intermediate_results = config['refilter_params'].get('redundant_intermediate_results', 'D').upper().strip() or 'D'
        
        if redundant_intermediate_results not in ['D','K']:
            logging.error(f"Error: Invalid 'redundant_intermediate_results' '{redundant_intermediate_results}' parameter in the [refilter_params] section. It must be one of 'D' or 'K'.")
            sys.exit(1)

        if not os.path.isfile(f"{project_id}/custom_log_{project_id}.txt"):
            logging.error(f"The 'custom_log_{project_id}.txt' file generated after the previous run is missing in the folder {project_id}.\n")
            print("The operation to re-filter previous results needs to be performed in the directory where the pipline was run previously, and do not modify the previously generated result files.\n")
            sys.exit(1)

        # Read old INI content from log.txt
        if os.path.isfile(f"{project_id}/custom_log_{project_id}.txt"):
            extract_ini_from_log(f"{project_id}/custom_log_{project_id}.txt", refilter_id)

        # Load external configuration file
        config_old = read_ini_file(f"old_log_{project_id}.ini")
        project_id_old = config_old['general']['project_id']
        mode_old = config_old['general']['mode']    # 用于辨别之前的运行模式是“A”还是“C”

        # 检验 project_id 是否有效。无效时，退出程序。
        if is_valid_project_id(project_id):
            logging.info("The 'project_id' is valid.")
        else:
            logging.error("Invalid 'project_id'. The 'project_id' must consist of letters, numbers, and underscores, and must start with a letter.\n")
            sys.exit(1)

        if is_valid_project_id(refilter_id):
            logging.info("The 'refilter_id' is valid.")
        else:
            logging.error("Invalid 'refilter_id'. The 'refilter_id' must consist of letters, numbers, and underscores, and must start with a letter.\n")
            sys.exit(1)

        if is_valid_project_id(project_id_old):
            logging.info("The 'project_id_old' is valid.")
        else:
            logging.error("The previous 'refilter_id' may be changed.\n")
            sys.exit(1)

        if project_id != project_id_old:
            logging.error("The previous 'project_id' in the [general] section may be changed.\n")
            print("The operation to re-filter previous results needs to be performed in the directory where the pipline was run previously, and do not modify the names of previously generated result files.\n")
            sys.exit(1)

        if project_id == refilter_id:
            logging.info(f"The 'refilter_id' in the [refilter_params] section cannot be the same as the 'project_id' in the [general] section.")

        log_file_name = f"refilter_log_{refilter_id}.txt"     # 设置配置log文件
        handlers = setup_refilter_logging(mode = 'a', refilter_id=refilter_id, log_file_name=log_file_name, include_file_handler = True)    # Append to log.txt for this case 

        signal.signal(signal.SIGINT, signal_handler)
        
        # 检测输入的resume参数，重过滤模式下，该参数不能用
        if args.resume:
            logging.warning("The '-resume' parameter is ineffective in Refiltering repeat info.\n")
            exit(1)
            
        # General parameters
        inputfasta = config['general']['inputfasta']
        seqdepth_alignment_software = config['sequencing_depth'].get('alignment_software').lower().strip()
        redundant_intermediate_results = config['general'].get('redundant_intermediate_results').upper().strip()

        # 从 [check_spanning_reads] 节获取参数
        check_spanning_length_old = config['check_spanning_reads'].getint('spanning_read_flanking_repeat_length')
        check_spanning_length_again = config['refilter_params'].getint('spanning_read_flanking_repeat_length', 5)            # 读取新旧跨越重复单元的序列长度，用于重新过滤结果
        try:      #检查跨越重复序列的长度参数是否合法
            validate_positive_integer(check_spanning_length_old, 'check_spanning_length')
        except ValueError as e:
            logging.error(f"An error occurred: {e}")
        try:      #检查跨越重复序列的长度参数是否合法
            validate_positive_integer(check_spanning_length_again, 'check_spanning_length')
        except ValueError as e:
            logging.error(f"An error occurred: {e}")
            
        # 处理线程数
        seqdepth_threads = config['sequencing_depth'].get('threads', '').strip()
        if seqdepth_threads.strip() == '':
            total_cores = os.cpu_count() or 1  # 处理 os.cpu_count() 可能返回 None 的情况
            seqdepth_threads = max(1, min(total_cores - 1, int(total_cores * 0.95)))
            logging.info(f"No thread parameter was provided, default parameter has been used instead: {seqdepth_threads} threads.")
        else:
            # 如果不是空字符串，则尝试转换为整数
            try:
                seqdepth_threads = int(float(seqdepth_threads)) if int(float(seqdepth_threads)) <= os.cpu_count() else None      # if 用于防止设置的进程数 > 计算机自身的进程数
            except ValueError:
                # 如果转换失败，则使用默认值，并记录错误信息
                total_cores = os.cpu_count() or 1
                seqdepth_threads = max(1, min(total_cores - 1, int(total_cores * 0.95)))
                logging.error(f"No thread parameter was provided, default parameter has been used instead: {seqdepth_threads} threads.")

        # 当mode为'A'时，检查 inputfasta 是否为空
        if not inputfasta:
            logging.error("When 'refilter_mode' in the [refilter_params] section is 'y/yes', the 'inputfasta' in section [general] cannot be empty.")
            sys.exit(1)
        count = check_fasta_sequence(inputfasta)
        if count == 0:
            logging.error(f"No sequences in the {inputfasta} file.")  
            sys.exit(1)                
        if count == 1:
            logging.info(f"There are {count} sequences in the input fasta file {inputfasta}.")
        if count > 1:
            cat_inputfasta_path, headers, lengths = concatenate_fasta(inputfasta, f"cat_inputfasta_{project_id}.fasta")
            logging.info(f"There are {count} sequences in the input fasta file {inputfasta}. They are concatenated into a pseudo-genome for further process.\n")
            # 使用itertools.accumulate进行  各条染色体长度的累积求和  
            cumulative_length = list(accumulate(lengths)) 
            time.sleep(5)

        # Define the keys to compare
        keys_to_compare = {
            'general': ['project_id', 'inputfasta', 'genome_type', 'complementary_chain', 'redundant_intermediate_results'],
            'ROUSFinder': ['repeat_length', 'reward', 'penalty'],
            'manually_calibrated_repeat_info': ['calibrated_repeat_file'],
            'mainconfiguration': ['flanked_sequence_length'],
            'subconfiguration': ['flanked_sequence_length'],
            'sequencing_depth': ['alignment_software', 'evalue', 'threads', 'NGS_single_end', 'NGS_pair_ends', 'TGS', 'TGS_type', 'filter_reads'],
            'check_spanning_reads': ['spanning_read_number'],
            }

        # Read old INI content from log.txt
        if os.path.isfile(f"{project_id}/custom_log_{project_id}.txt"):
            #extract_ini_from_log(f"{project_id}/custom_log_{project_id}.txt", refilter_id)
            # 前后两个配置的内容进行比较
            has_difference = compare_ini_contents(f"old_log_{project_id}.ini", args.config, keys_to_compare)
            if has_difference:
                logging.error(f"Previous '.ini' file {args.config} has been changed, except 'mode' and ['refilter_params'] section.")
                logging.info("The current work has been terminated.\n")
                sys.exit(1)     # 配置信息出现变化，则终止程序
        else:
            logging.error("The log file generated during the previous program run was not found.") 
            logging.info("The current work was terminated.\n")
            sys.exit(1)     # 之前的log文件不存在，终止程序

        if os.path.isfile(f"old_log_{project_id}.ini"):      # 删除用于检查配置文件有没有发生变化的临时配置文件
            os.remove(f"old_log_{project_id}.ini")

################################################################################
        # Get the directory of the current path
        dir_path = os.path.dirname(os.path.realpath(__file__))

        complementary_chain = config['general'].get('complementary_chain').upper().strip()
        if complementary_chain not in ['YES','Y','NO','N']:
            logging.error(f"Error: Invalid 'complementary_chain' value '{complementary_chain}'. It must be one of 'Yes', 'Y', 'No' or 'N'.")
            sys.exit(1)

        read_number_again = config['refilter_params'].getint('spanning_read_number', 5)
        try:    #检查read number 参数是否合法
            validate_positive_integer(read_number_again, 'read_number_again')
        except ValueError as e:
            logging.error(f"An error occurred: {e}")
            sys.exit(1)

        # read_number 在删除中间结果的模式下要不小于 之前设置的值
        read_number_old = config['check_spanning_reads'].getint('spanning_read_number')
        if float(read_number_again) < float(read_number_old) and redundant_intermediate_results == 'D':  
            logging.error(f"The 'spanning_read_number' in the [refilter_params] section should be >= that in the [check_spanning_reads] section.")
            sys.exit(1)

        # spanning length 在删除中间结果的模式下要不小于 之前设置的值
        if float(check_spanning_length_again) < float(check_spanning_length_old) and redundant_intermediate_results == 'D': 
            logging.error(f"The 'spanning_read_flanking_repeat_length' in the [refilter_params] section should be >= that in the [check_spanning_reads] section.")
            sys.exit(1)
            
        # 检测 spanning length 和 reads number 是否被重新设置
        is_refilter_logic = True          # 设置过滤的开关 is_spanning_length_same 
        if check_spanning_length_old == check_spanning_length_again and read_number_again == read_number_old:      # 再次筛选的条件未发生变化，无需再次筛选
            is_refilter_logic = False
            logging.error(f"Core re-filtering parameters 'spanning_read_number', 'spanning_read_flanking_repeat_length' and 'repeat_length' in the [refilter_params] section should be reset one or more!\n")
            exit(1)
        else:
            logging.info(f"Next step is to refilter repeat pairs supporting genomic recombination based on new parameters: 'spanning_read_number {read_number_again}' and 'spanning_read_flanking_repeat_length {check_spanning_length_again}' in the [refilter_params] section.")
            print()
        
################################################################################
################################################################################
        # 读取存放结果的路径
        extrsubcon_output_dir_prefix = config['subconfiguration'].get('output_directory_prefix', 'subconfig_repeat-spanned_results').strip() or 'subconfig_repeat-spanned_results'
        extrmaincon_output_dir_prefix = config['mainconfiguration'].get('output_directory_prefix', 'mainconfig_repeat-spanned_results').strip() or 'mainconfig_repeat-spanned_results'
        final_output_dir_sub = f"{extrsubcon_output_dir_prefix}_{project_id}/filter_subconfig_spanned_results.tsv"
        final_output_dir_main = f"{extrmaincon_output_dir_prefix}_{project_id}/filter_mainconfig_spanned_results.tsv"
        
        # 为再次过滤的结果建立新的存储文件夹 refilter_id
        output_dir_prefix_filter = config['refilter_params'].get('output_directory_prefix', 'refiltered_repeat-spanning_results').strip() or 'refiltered_repeat-spanning_results'
        output_dir_again = output_dir_prefix_filter + "_" + refilter_id       # + "/" 

        # 如果存在与 output_dir_again 同名的文件，需要删除该文件以便创建目录
        if os.path.exists(output_dir_again) and not os.path.isdir(output_dir_again):
            os.remove(output_dir_again)  # 删除与目录同名的文件        
        
        # 构造完整的文件夹路径
        full_output_dir_again = os.path.join(project_id, output_dir_again)
        pseudo_output_dir_again = f"pseudo_{output_dir_again}"
        full_pseudo_output_dir_again = os.path.join(project_id, f"pseudo_{output_dir_again}")

        # 检查文件夹是否存在，并根据 args.redo 决定是否删除
        def check_and_delete_dirs(args, output_dir_again, full_output_dir_again, pseudo_output_dir_again, full_pseudo_output_dir_again):
            # 检查并记录存在的文件夹
            existing_dirs = []
            if os.path.isdir(output_dir_again):
                existing_dirs.append(output_dir_again)
            if os.path.isdir(full_output_dir_again):
                existing_dirs.append(full_output_dir_again)
            if os.path.isdir(pseudo_output_dir_again):
                existing_dirs.append(pseudo_output_dir_again)
            if os.path.isdir(full_pseudo_output_dir_again):
                existing_dirs.append(full_pseudo_output_dir_again)

            # 如果存在文件夹，给出警告
            if existing_dirs and not args.redo:
                logging.warning(f"Previous refiltered results '{existing_dirs}' were found. Please use '-redo' parameter to continue or reset the 'refilter_id'.\n")
                sys.exit(1)

            # 如果使用 -redo 参数，则删除文件夹
            if existing_dirs and args.redo:
                logging.warning(f"Existing directories '{existing_dirs}' will be deleted. You can TERMINATE the program using 'Ctrl + C' within 10s to cancel.")
                try:
                    # 等待10秒，给用户一个取消的机会
                    time.sleep(10)
                    for dir in existing_dirs:
                        if os.path.isdir(dir):
                            shutil.rmtree(dir)
                    logging.info("Deleted existing directories.")
                except KeyboardInterrupt:
                    logging.info("Program terminated by the user.\n")
                    sys.exit(1)

        # 调用函数, 处理前一次运行留下来的结果
        check_and_delete_dirs(args, output_dir_again, full_output_dir_again, pseudo_output_dir_again, full_pseudo_output_dir_again)
        
        if not os.path.isdir(output_dir_again):
            os.makedirs(output_dir_again)      # 如果存储重新过滤结果的文件夹不存在，则重新建立

        if mode_old == 'A':
            inputfasta = config['general']['inputfasta']
            
            rous_output_file_prefix = config['ROUSFinder'].get('output_file_prefix')
            # Check if rous_output_file_prefix is empty, if so, use the prefix of inputfasta
            if not rous_output_file_prefix:
                basename_input = os.path.basename(inputfasta)  
                rous_output_file_prefix =  '.'.join(basename_input.split('.')[:-1]) 

            if count ==1: 
                check_files = [
                    f"{project_id}/ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_5CT.tsv",
                    f"{project_id}/ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt"]
                # 遍历文件列表，检查文件是否存在，如果不存在则报错退出
                for file_path in check_files:
                    if not os.path.isfile(file_path):
                        logging.error(f"Previous results {file_path} may be deleted or renamed!")
                        sys.exit(1)
            if count >1: 
                check_files = [
                    f"{project_id}/ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_5CT.tsv",
                    f"{project_id}/ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt",
                    ]
                # 遍历文件列表，检查文件是否存在，如果不存在则报错退出
                for file_path in check_files:
                    if not os.path.isfile(file_path):
                        logging.error(f"Previous results {file_path} may be deleted or renamed!")
                        sys.exit(1)

        if mode_old == 'C': 
            # 手工矫正的结果，修改了表头，通过检查格式，输出结果文件名为 adjusted_manual_calibrated_list_{project_id}.tsv
            # 并复制到 ROUSFinder_results 文件夹内，重新过滤结果时，ROUSFinder_results 被转移到了文件夹 project_id 内。
            manually_calibrate_adj = f"ROUSFinder_results_{project_id}/adjusted_manual_calibrated_list.tsv"
            # manually_calibrate = config['manually_calibrated_repeat_info'].get('calibrated_repeat_file', '').strip()
            # 以上 manually_calibrate_adj 对应的文件 adjusted_manual_calibrated_list.tsv
            # 以上 manually_calibrate_reset 对应的文件 adjusted_manual_calibrated_list_{project_id}_rep_calibration_table.txt
            
            dir_path_manual, file_name_manual = os.path.split(manually_calibrate_adj)
            prefix_manual = os.path.splitext(file_name_manual)[0]
            manually_calibrate_reset = f"ROUSFinder_results_{project_id}/{prefix_manual}_{project_id}_rep_calibration_table.txt"
            
            rous_output_file_prefix = config['ROUSFinder'].get('output_file_prefix')
            # Check if rous_output_file_prefix is empty, if so, use the prefix of inputfasta     # rous_output_file_prefix
            if not rous_output_file_prefix:
                basename_input = os.path.basename(manually_calibrate_adj)  
                rous_output_file_prefix =  '.'.join(basename_input.split('.')[:-1])         # rous_output_file_prefix 为文件 adjusted_manual_calibrated_list_{project_id}.tsv 的前缀
            
            # 当mode为'C'时，检查manually_calibrate是否为空
            if not manually_calibrate_adj:
                logging.error("The previous running mode was 'C', so, the 'inputfile' in the [manually_calibrated_repeat_info] section cannot be empty.")
                sys.exit(1)

            if count ==1: 
                check_files = [
                    f"{project_id}/ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_5CT.tsv",
                    f"{project_id}/ROUSFinder_results_{project_id}/{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt",
                    ]
                # 遍历文件列表，检查文件是否存在，如果不存在则报错退出
                for file_path in check_files:
                    if not os.path.isfile(file_path):
                        logging.error(f"Previous results {file_path} may be deleted or renamed!")
                        sys.exit(1)
            if count >1: 
                check_files = [
                    f"{project_id}/ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_5CT.tsv",
                    f"{project_id}/ROUSFinder_results_{project_id}/pseudo_{rous_output_file_prefix}_{project_id}_rep_calibration_table.txt"]
                # 遍历文件列表，检查文件是否存在，如果不存在则报错退出
                for file_path in check_files:
                    if not os.path.isfile(file_path):
                        logging.error(f"Previous results {file_path} may be deleted or renamed!")
                        sys.exit(1)

################################################################################
        # 当用户重新设置了 read number 和 spanning length时，进行重过滤
        if is_refilter_logic:      
            # 对次要构型进行重过滤
            map_results_folders = glob.glob(f"{project_id}/{extrsubcon_output_dir_prefix}_{project_id}/*_results")    # 找到存储所有结果的文件夹
            file_path_sub = os.path.join(project_id, extrsubcon_output_dir_prefix + "_" + project_id, "filter_subconfig_spanned_results.tsv")
            map_results_folders_count = len(map_results_folders)

            print()
            logging.info(f"@@@@@@@@@@ There are {map_results_folders_count} repeat pairs that need to be refiltered. @@@@@@@@@@\n\n")
            time.sleep(5)
            logging.info(f"%%%%%%%%%%%% Start to refilter {map_results_folders_count} repeat units in the subconfiguration! %%%%%%%%%%%%")
            logging.info(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            time.sleep(2)

            header = "sub_trimmed_seq\tsub_trimmed_seq_length(bp)\tsub_spanning_read_number"
            # Check if the file exists or is empty, and if so, create it and write the header
            if not os.path.isfile(file_path_sub) or os.path.getsize(file_path_sub) == 0:
                with open(file_path_sub, 'w') as file:
                    file.write(header + '\n')

            # Immediately after, run check_mainconfig_spanning_reads.py for the output
            every_subcon_count = 1
            for map_results_folder in map_results_folders:       # 针对每一个存放结果的文件夹\
                # 获取文件夹中的 .fasta 文件
                fasta_files = glob.glob(os.path.join(map_results_folder, "*.fasta"))
                # 检查文件数量是否为1
                if len(fasta_files) != 1:
                    raise ValueError(f"Error: There should be exactly one '.fasta' file in the folder '{map_results_folder}'. Found {len(fasta_files)} files.")
                
                # 如果满足条件，继续处理
                fasta_file = fasta_files[0]
                fasta_file_suffix_basename = os.path.basename(fasta_file)
                
                ordinal_suffix = get_ordinal_suffix(every_subcon_count)  
                ordinal_suffix_basename = os.path.basename(ordinal_suffix)
                
                # 使用colorama高亮显示计数器和序数词后缀  
                highlighted_text = f"{Fore.RED}{every_subcon_count}{ordinal_suffix_basename}{Style.RESET_ALL}"
                print()
                logging.info(f"**** Start to process the {highlighted_text} sequence: {fasta_file_suffix_basename}. ****")  
                time.sleep(2)
                every_subcon_count += 1      # 计数器递增
                
                ############# 重新过滤,所用的信息全部来自"repeat-spanning_read.txt"
                spanning_read_file = os.path.join(map_results_folder, "repeat-spanning_read.txt")
                run_command(["python", os.path.join(dir_path, "count_subconfig_spanning_reads.py"), "-i", spanning_read_file, "-o", file_path_sub, "-cl", str(check_spanning_length_again), "-rn", str(read_number_again)])

            print()
            logging.info(f"@@@@@@@@@@ The {map_results_folders_count} repeat units in subconfiguration have been refiltered completely! @@@@@@@@@@\n\n")

################################################################################
        if is_refilter_logic:      
            ################# 对主要构型进行重过滤
            map_results_folders = glob.glob(f"{project_id}/{extrmaincon_output_dir_prefix}_{project_id}/*_results")
            file_path_main = os.path.join(project_id, extrmaincon_output_dir_prefix + "_" + project_id, "filter_mainconfig_spanned_results.tsv")
            
            header = "main_trimmed_seq\tmain_trimmed_seq_length(bp)\tmain_spanning_read_number"
            # Check if the file exists or is empty, and if so, create it and write the header
            if not os.path.isfile(file_path_main) or os.path.getsize(file_path_main) == 0:
                with open(file_path_main, 'w') as file:
                    file.write(header + '\n')

            map_results_folders_count = len(map_results_folders)
            logging.info(f"%%%%%%%%%%%% Start to refilter {map_results_folders_count} repeat units in the mainconfiguration! %%%%%%%%%%%%")
            logging.info(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            time.sleep(2)

            # Immediately after, run check_mainconfig_spanning_reads.py for the output
            every_maincon_count = 1
            for map_results_folder in map_results_folders:

                # 获取文件夹中的 .fasta 文件
                fasta_files = glob.glob(os.path.join(map_results_folder, "*.fasta"))
                # 检查文件数量是否为1
                if len(fasta_files) != 1:
                    raise ValueError(f"Error: There should be exactly one '.fasta' file in the folder '{map_results_folder}'. Found {len(fasta_files)} files.")
                
                # 如果满足条件，继续处理
                fasta_file = fasta_files[0]
                fasta_file_suffix_basename = os.path.basename(fasta_file)
                
                ordinal_suffix = get_ordinal_suffix(every_maincon_count)  
                ordinal_suffix_basename = os.path.basename(ordinal_suffix)
                
                # 使用colorama高亮显示计数器和序数词后缀  
                highlighted_text = f"{Fore.RED}{every_maincon_count}{ordinal_suffix_basename}{Style.RESET_ALL}"
                print()
                logging.info(f"**** Start to process the {highlighted_text} sequence: {fasta_file_suffix_basename}. ****")  
                time.sleep(2)
                every_maincon_count += 1      # 计数器递增
                
                ############# 重新过滤,所用的信息全部来自"repeat-spanning_read.txt"
                spanning_read_file = os.path.join(map_results_folder, "repeat-spanning_read.txt")
                run_command(["python", os.path.join(dir_path, "count_mainconfig_spanning_reads.py"), "-i", spanning_read_file, "-o", file_path_main, "-cl", str(check_spanning_length_again)])

            print()
            logging.info(f"@@@@@@@@@@ The {map_results_folders_count} repeat units in mainconfiguration have been refiltered completely! @@@@@@@@@@\n\n")

################################################################################
             # 对sam文件重新过滤后，再进行一次结果的梳理，形成输出结果的表格
            if count == 1 and mode_old == "A":         # 当 spanning_length 相等时
                run_command(["python", os.path.join(dir_path, "filter_spanning_results.py"), "-sr", project_id + "/" + extrsubcon_output_dir_prefix + "_" + project_id + "/filter_subconfig_spanned_results.tsv", "-mr", project_id + "/" + extrmaincon_output_dir_prefix + "_" + project_id + "/filter_mainconfig_spanned_results.tsv", "-rf", project_id + "/" + f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_5CT.tsv", "-rc", project_id + "/" + f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_calibration_table.txt", "-o", output_dir_prefix_filter + "_" + refilter_id + "/", "-rn", str(read_number_again), "-cr", complementary_chain])
            elif count > 1 and mode_old == "A":
                run_command(["python", os.path.join(dir_path, "filter_spanning_results.py"), "-sr", project_id + "/" + extrsubcon_output_dir_prefix + "_" + project_id + "/filter_subconfig_spanned_results.tsv", "-mr", project_id + "/" + extrmaincon_output_dir_prefix + "_" + project_id + "/filter_mainconfig_spanned_results.tsv", "-rf", project_id + "/" + f"ROUSFinder_results_{project_id}/" + "pseudo_" + rous_output_file_prefix + "_" + project_id + "_rep_5CT.tsv", "-rc", project_id + "/" + f"ROUSFinder_results_{project_id}/" + "pseudo_" + rous_output_file_prefix + "_" + project_id + "_rep_calibration_table.txt", "-o", output_dir_prefix_filter + "_" + refilter_id + "/", "-rn", str(read_number_again), "-cr", complementary_chain])

            if count == 1 and mode_old == "C":         # 当 spanning_length 相等时
                run_command(["python", os.path.join(dir_path, "filter_spanning_results.py"), "-sr", project_id + "/" + extrsubcon_output_dir_prefix + "_" + project_id + "/filter_subconfig_spanned_results.tsv", "-mr", project_id + "/" + extrmaincon_output_dir_prefix + "_" + project_id + "/filter_mainconfig_spanned_results.tsv", "-rf", project_id + "/" + f"ROUSFinder_results_{project_id}/" + rous_output_file_prefix + "_" + project_id + "_rep_5CT.tsv", "-rc", project_id + "/" + manually_calibrate_reset, "-o", output_dir_prefix_filter + "_" + refilter_id + "/", "-rn", str(read_number_again), "-cr", complementary_chain])
            elif count > 1 and mode_old == "C":
                run_command(["python", os.path.join(dir_path, "filter_spanning_results.py"), "-sr", project_id + "/" + extrsubcon_output_dir_prefix + "_" + project_id + "/filter_subconfig_spanned_results.tsv", "-mr", project_id + "/" + extrmaincon_output_dir_prefix + "_" + project_id + "/filter_mainconfig_spanned_results.tsv", "-rf", project_id + "/" + f"ROUSFinder_results_{project_id}/" + "pseudo_" + rous_output_file_prefix + "_" + project_id + "_rep_5CT.tsv", "-rc", project_id + "/" + manually_calibrate_reset, "-o", output_dir_prefix_filter + "_" + refilter_id + "/", "-rn", str(read_number_again), "-cr", complementary_chain])
                

################################################################################
#########################################
        # 对过滤的结果进行 start 和 end 的溯源
        if count > 1:
            #重新归档 过滤后产生的新结果文件
            if os.path.isfile(f"refiltered_repeat-spanning_results_{refilter_id}/paired_repeats_for_mapping.tsv"):
                file_path = f"refiltered_repeat-spanning_results_{refilter_id}/paired_repeats_for_mapping.tsv"
                if is_file_empty_or_invisible_chars(file_path):    
                    default_headers = ['fragment_id', 'length', 'start', 'end', 'direction', 'chromosome', 'plus_ratio(s/m)', 'minus_ratio(s/m)', 'combined_ratio', 'type', 'paired_id', 'paired_length', 'paired_start', 'paired_end', 'paired_direction', 'paired_chromosome', 'paired_plus_ratio(s/m)', 'paired_minus_ratio(s/m)', 'paired_combined_ratio', 'spanning_read_mcfg']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                        
                df = pd.read_csv(file_path, sep='\t') 
                df = df.dropna(how='all')  # 去掉空白行
    
                if len(df) == 0:
                    logging.warning(f"The DataFrame {file_path} is empty!!")
                    shutil.copy(f"refiltered_repeat-spanning_results_{refilter_id}/paired_repeats_for_mapping.tsv",   
                    f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_paired_repeats_for_mapping.tsv")  
                else:
                    os.rename(f"refiltered_repeat-spanning_results_{refilter_id}/paired_repeats_for_mapping.tsv",   
                              f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_paired_repeats_for_mapping.tsv")  
                    file_path = f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_paired_repeats_for_mapping.tsv"  
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')  # 去掉空白行
    
                    # 应用函数并分别更新列  
                    temp_result = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  

                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = df['adjusted_start']  
                    df['end'] = df['adjusted_end']  
                    
                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])  

                    # 将chromosome列插入到direction列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('direction') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  

                    # 更新cols以避免重复列  
                    cols = df.columns.tolist()  

                    # 处理 paired_start 和 paired_end  
                    df[['paired_chromosome', 'adjusted_paired_start', 'adjusted_paired_end']] = df.apply(lambda row: find_chromosome(row['paired_start'], row['paired_end'], cumulative_length), axis=1, result_type='expand')  
                    df['paired_start'] = df['adjusted_paired_start']  
                    df['paired_end'] = df['adjusted_paired_end']  
                    df = df.drop(columns=['adjusted_paired_start', 'adjusted_paired_end'])  

                    # 将paired_chromosome列插入到paired_direction列之后  
                    paired_chromosome_index = cols.index('paired_direction') + 1  
                    cols.insert(paired_chromosome_index, 'paired_chromosome')  
                    df = df[cols]  

                    # 删除最后1列
                    if len(df.columns) > 1:
                        df = df.iloc[:, :-1]
                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]

                    # 保存修改后的 DataFrame  
                    file_path = f"refiltered_repeat-spanning_results_{refilter_id}/paired_repeats_for_mapping.tsv"
                    df.to_csv(file_path, sep='\t', index=False)  

            if os.path.isfile(f"refiltered_repeat-spanning_results_{refilter_id}/paired_repeats_recomb-supporting_ratio.tsv"):
                file_path = f"refiltered_repeat-spanning_results_{refilter_id}/paired_repeats_recomb-supporting_ratio.tsv"
                if is_file_empty_or_invisible_chars(file_path):    
                    default_headers = ['fragment_id', 'length', 'start', 'end', 'direction', 'chromosome', 'plus_ratio(s/m)', 'minus_ratio(s/m)', 'combined_ratio', 'type', 'paired_id', 'paired_length', 'paired_start', 'paired_end', 'paired_direction', 'paired_chromosome', 'paired_plus_ratio(s/m)', 'paired_minus_ratio(s/m)', 'paired_combined_ratio', 'spanning_read_mcfg']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                
                df = pd.read_csv(file_path, sep='\t') 
                df = df.dropna(how='all')     # 去掉空白行
                
                if len(df) == 0:
                    logging.warning(f"The DataFrame {file_path} is empty!!")
                    shutil.copy(f"refiltered_repeat-spanning_results_{refilter_id}/paired_repeats_recomb-supporting_ratio.tsv",   
                                f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_paired_repeats_recomb-supporting_ratio.tsv")  
                else:
                    os.rename(f"refiltered_repeat-spanning_results_{refilter_id}/paired_repeats_recomb-supporting_ratio.tsv",   
                              f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_paired_repeats_recomb-supporting_ratio.tsv")  
                    file_path = f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_paired_repeats_recomb-supporting_ratio.tsv"  
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')     # 去掉空白行
            
                    # 应用函数并分别更新列  
                    temp_result = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = df['adjusted_start']  
                    df['end'] = df['adjusted_end']  
                    
                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])  
            
                    # 将chromosome列插入到direction列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('direction') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  
            
                    # 删除重复的chromosome列
                    df = df.loc[:,~df.columns.duplicated()]
            
                    # 处理 paired_start 和 paired_end  
                    df[['paired_chromosome', 'adjusted_paired_start', 'adjusted_paired_end']] = df.apply(lambda row: find_chromosome(row['paired_start'], row['paired_end'], cumulative_length), axis=1, result_type='expand')  
                    df['paired_start'] = df['adjusted_paired_start']  
                    df['paired_end'] = df['adjusted_paired_end']  

                    df = df.drop(columns=['adjusted_paired_start', 'adjusted_paired_end'])  
            
                    # 将paired_chromosome列插入到paired_direction列之后  
                    cols = df.columns.tolist()  
                    paired_chromosome_index = cols.index('paired_direction') + 1  
                    cols.insert(paired_chromosome_index, 'paired_chromosome')  
                    df = df[cols]  
                    
                    # 删除多余的最后两列
                    df = df.iloc[:, :-2]
            
                    # 保存修改后的 DataFrame  
                    file_path = f"refiltered_repeat-spanning_results_{refilter_id}/paired_repeats_recomb-supporting_ratio.tsv"  
                    df.to_csv(file_path, sep='\t', index=False)  
            
            if os.path.isfile(f"refiltered_repeat-spanning_results_{refilter_id}/one_repeat_unit_without_spanning_reads.tsv"):
                file_path = f"refiltered_repeat-spanning_results_{refilter_id}/one_repeat_unit_without_spanning_reads.tsv"
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path):    
                    default_headers = ['fragment_id', 'length', 'start', 'end', 'direction', 'chromosome', 'plus_ratio(s/m)', 'minus_ratio(s/m)', 'combined_ratio', 'type', 'paired_id', 'paired_length', 'paired_start', 'paired_end', 'paired_direction', 'paired_chromosome', 'paired_plus_ratio(s/m)', 'paired_minus_ratio(s/m)', 'paired_combined_ratio', 'spanning_read_mcfg']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                
                df = pd.read_csv(file_path, sep='\t')
                df = df.dropna(how='all')  # 去掉空白行
                
                if len(df) == 0:
                    logging.warning(f"The DataFrame {file_path} is empty!!")
                    shutil.copy(f"refiltered_repeat-spanning_results_{refilter_id}/one_repeat_unit_without_spanning_reads.tsv", 
                                f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_one_repeat_unit_without_spanning_reads.tsv")
                else:
                    os.rename(f"refiltered_repeat-spanning_results_{refilter_id}/one_repeat_unit_without_spanning_reads.tsv",
                              f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_one_repeat_unit_without_spanning_reads.tsv")
                    file_path = f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_one_repeat_unit_without_spanning_reads.tsv"
                    df = pd.read_csv(file_path, sep='\t')
                    df = df.dropna(how='all')  # 去掉空白行
            
                    # 应用函数并分别更新列
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')
                    df['start'] = df['adjusted_start']
                    df['end'] = df['adjusted_end']

                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])
            
                    # 将chromosome列插入到direction列之后
                    cols = df.columns.tolist()
                    chromosome_index = cols.index('direction') + 1
                    cols.insert(chromosome_index, 'chromosome')
                    df = df[cols]
            
                    # 处理 paired_start 和 paired_end
                    df[['paired_chromosome', 'adjusted_paired_start', 'adjusted_paired_end']] = df.apply(lambda row: find_chromosome(row['paired_start'], row['paired_end'], cumulative_length), axis=1, result_type='expand')
                    df['paired_start'] = df['adjusted_paired_start']
                    df['paired_end'] = df['adjusted_paired_end']

                    df = df.drop(columns=['adjusted_paired_start', 'adjusted_paired_end'])
            
                    # 将paired_chromosome列插入到paired_direction列之后
                    paired_chromosome_index = cols.index('paired_direction') + 1
                    cols.insert(paired_chromosome_index, 'paired_chromosome')
                    df = df[cols]
            
                    # 删除最后1列
                    if len(df.columns) > 1:
                        df = df.iloc[:, :-1]
                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]
            
                    # 保存修改后的 DataFrame
                    file_path = f"refiltered_repeat-spanning_results_{refilter_id}/one_repeat_unit_without_spanning_reads.tsv"
                    df.to_csv(file_path, sep='\t', index=False)

            if os.path.isfile(f"refiltered_repeat-spanning_results_{refilter_id}/one_chain_without_sufficient_spanning_reads.tsv"):  
                file_path = f"refiltered_repeat-spanning_results_{refilter_id}/one_chain_without_sufficient_spanning_reads.tsv"  
                # 检查文件大小，确定是否为空文件  
                if is_file_empty_or_invisible_chars(file_path): 
                    default_headers = ['fragment_id', 'length', 'start', 'end', 'direction', 'chromosome', 'plus_ratio(s/m)', 'minus_ratio(s/m)', 'combined_ratio', 'type', 'paired_id', 'paired_length', 'paired_start', 'paired_end', 'paired_direction', 'paired_chromosome', 'paired_plus_ratio(s/m)', 'paired_minus_ratio(s/m)', 'paired_combined_ratio', 'spanning_read_mcfg']
                    # 文件为空，将默认表头写入文件  
                    with open(file_path, 'w', newline='', encoding='utf-8') as f:  
                        f.write('\t'.join(default_headers) + '\n') 
                
                df = pd.read_csv(file_path, sep='\t') 
                df = df.dropna(how='all')  # 去掉空白行
                
                if len(df) == 0:
                    logging.warning(f"The DataFrame {file_path} is empty!!")
                    shutil.copy(f"refiltered_repeat-spanning_results_{refilter_id}/one_chain_without_sufficient_spanning_reads.tsv",   
                                f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_one_chain_without_sufficient_spanning_reads.tsv")  
                else:
                    os.rename(f"refiltered_repeat-spanning_results_{refilter_id}/one_chain_without_sufficient_spanning_reads.tsv",   
                              f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_one_chain_without_sufficient_spanning_reads.tsv")  
                    file_path = f"refiltered_repeat-spanning_results_{refilter_id}/pseudo_one_chain_without_sufficient_spanning_reads.tsv"  
                    df = pd.read_csv(file_path, sep='\t')  
                    df = df.dropna(how='all')  # 去掉空白行
                  
                    # 应用函数并分别更新列  
                    df[['chromosome', 'adjusted_start', 'adjusted_end']] = df.apply(lambda row: find_chromosome(row['start'], row['end'], cumulative_length), axis=1, result_type='expand')  
                    df['start'] = df['adjusted_start']  
                    df['end'] = df['adjusted_end']  

                    df = df.drop(columns=['adjusted_start', 'adjusted_end'])  
                    
                    # 将chromosome列插入到direction列之后  
                    cols = df.columns.tolist()  
                    chromosome_index = cols.index('direction') + 1  
                    cols.insert(chromosome_index, 'chromosome')  
                    df = df[cols]  

                    # 处理 paired_start 和 paired_end  
                    df[['paired_chromosome', 'adjusted_paired_start', 'adjusted_paired_end']] = df.apply(lambda row: find_chromosome(row['paired_start'], row['paired_end'], cumulative_length), axis=1, result_type='expand')  
                    df['paired_start'] = df['adjusted_paired_start']  
                    df['paired_end'] = df['adjusted_paired_end']  

                    df = df.drop(columns=['adjusted_paired_start', 'adjusted_paired_end'])  

                    # 将paired_chromosome列插入到paired_direction列之后  
                    paired_chromosome_index = cols.index('paired_direction') + 1  
                    cols.insert(paired_chromosome_index, 'paired_chromosome')  
                    df = df[cols]  

                    # 删除最后1列
                    if len(df.columns) > 1:
                        df = df.iloc[:, :-1]
                    # 删除重复的列名
                    if df.columns.duplicated().sum() > 0:
                        df = df.loc[:, ~df.columns.duplicated()]

                    # 保存修改后的 DataFrame  
                    file_path = f"refiltered_repeat-spanning_results_{refilter_id}/one_chain_without_sufficient_spanning_reads.tsv"  
                    df.to_csv(file_path, sep='\t', index=False)

################################################################################
        if mode == "N" and refilter_mode:
            shutil.rmtree(f"{project_id}/{output_dir_again}") if os.path.isdir(f"{project_id}/{output_dir_again}") else None
            shutil.move(f"{output_dir_again}", f"{project_id}") if os.path.isdir(f"{output_dir_again}") else None
            if count > 1:
                if not os.path.isdir(f"{project_id}/pseudo_{output_dir_again}"):
                    os.mkdir(f"{project_id}/pseudo_{output_dir_again}")           # 创建文件夹，存放基于伪基因组的结果
                
                # 定义源文件模式和目标文件夹，移动pseudo_开头的文件
                source_pattern = f"{project_id}/{output_dir_again}/pseudo_*"
                destination_dir = f"{project_id}/pseudo_{output_dir_again}"
                # 查找匹配的文件并逐个移动
                for file in glob.glob(source_pattern):
                    shutil.move(file, destination_dir) 
                
            os.remove(f"cat_inputfasta_{project_id}.fasta") if os.path.isfile(f"cat_inputfasta_{project_id}.fasta") else None
            os.remove(f"custom_log_{project_id}.txt") if os.path.isfile(f"custom_log_{project_id}.txt") else None                   # 删除当前目录下的两个中间结果
            shutil.move(f"refilter_log_{refilter_id}.txt", f"refiltered_repeat-spanning_results_{refilter_id}") if os.path.isfile(f"refilter_log_{refilter_id}.txt") else None
            
            if redundant_intermediate_results == "D":       # 删除中间结果 处理多条序列的 伪基因组的 中间结果
                shutil.rmtree(f"{project_id}/pseudo_{output_dir_again}") if os.path.isdir(f"{project_id}/pseudo_{output_dir_again}") else None
            
            os.remove(f"{output_dir_again}") if os.path.isfile(f"{output_dir_again}") else None        # 不知道是哪里多了一个 空的 文件refiltered_repeat-spanning_results_flt，删除
            
            logging.info(f"Check the refiltered repeat info in '{output_dir_prefix_filter}_{refilter_id}/paired_repeats_for_mapping.tsv' about mitogenomic recombination.\n")

######################################################################################################################################################
######################################################################################################################################################
if __name__ == "__main__":
    start_time = time.time()     # Capture the start time
    main()
    end_time = time.time()  # Capture the end time when execution is about to finish
    
    execution_time = end_time - start_time  # Calculate the difference to get execution time
    print(f"Total execution time: {execution_time} seconds.\n")  # Log the execution time


