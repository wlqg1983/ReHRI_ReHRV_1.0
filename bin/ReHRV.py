#!/usr/bin/env python

import sys, os, configparser, re
import subprocess, shutil, time
import time, argparse
import pandas as pd
import logging, glob
import chardet
import codecs

#######################################################################################################################################################################
def convert_file_to_utf8(file_path):
    """
    Convert file encoding to UTF-8 (no return value, directly modifies the source file)
    Automatically handles BOM headers, line endings, and encoding exceptions
    
    Args:
        file_path (str): Path to the file to be converted
        
    Raises:
        UnicodeDecodeError: When the file cannot be decoded with any encoding
        IOError: When file read/write errors occur
    """
    # Phase 1: Detect file encoding
    with open(file_path, 'rb') as f:
        raw_data = f.read(10000)
        result = chardet.detect(raw_data)
        original_encoding = result['encoding']
        print(f"Detected original encoding: {original_encoding} (confidence: {result['confidence']:.2f})")

    # Phase 2: Read file content (try multiple encodings)
    content = None
    encodings_to_try = [
        'utf-8-sig',  # Try UTF-8 with BOM first
        original_encoding,  # Then try the detected encoding
        'gb18030',     # Common encoding in Chinese environments
        'latin-1'      # Final fallback
    ]
    
    for enc in encodings_to_try:
        try:
            with open(file_path, 'r', encoding=enc) as f:
                content = f.read()
            print(f"Successfully read file with encoding [{enc}]")
            original_encoding = enc
            break
        except UnicodeDecodeError:
            continue
    
    if content is None:
        raise UnicodeDecodeError(f"Unable to decode file: {file_path}")

    # Phase 3: Convert non-UTF-8 files
    if original_encoding.lower() not in ['utf-8', 'utf8', 'utf-8-sig']:
        print("\nWarning: Detected non-UTF-8 encoded file!")
        print(f"File path: {os.path.abspath(file_path)}")
        print(f"Current encoding: {original_encoding}")
        print("This operation will:")
        print("1. Permanently convert the file encoding to UTF-8")
        print("2. Overwrite the original file")
        print("\nPlease ensure you have backed up the original file!")
        
        # 10-second countdown confirmation
        try:
            for i in range(10, 0, -1):
                print(f"\rAuto-proceeding in {i} seconds... (Press Ctrl+C to cancel)", end='')
                sys.stdout.flush()
                time.sleep(1)
            print("\nStarting conversion...")
        except KeyboardInterrupt:
            print("\nOperation cancelled by user")
            return

        # Phase 4: Safely write to a temporary file
        temp_file = file_path + '.tmp'
        try:
            # Normalize line endings to Linux format (\n)
            normalized_content = content.replace('\r\n', '\n').replace('\r', '\n')
            
            with open(temp_file, 'w', encoding='utf-8', newline='\n') as f:
                f.write(normalized_content)
                
            # Atomically replace the original file
            os.replace(temp_file, file_path)
            print(f"File successfully converted to UTF-8 encoding (no BOM)")
            
        except Exception as e:
            if os.path.exists(temp_file):
                os.remove(temp_file)
            raise IOError(f"File conversion failed: {str(e)}")
    else:
        print(f"File {file_path} is already UTF-8 encoded, no conversion needed")

#######################################################################################################################################################
def is_valid_project_id(project_id):
    # 检查 project_id 是否有效：只能包含字母、数字、下划线， 必须以字母开头
    # 定义正则表达式
    pattern = r'^[A-Za-z][A-Za-z0-9_]*$'
    
    # 使用 re.match 检查 project_id 是否匹配正则表达式
    if re.match(pattern, project_id):
        return True
    else:
        return False
        
##########################################################################################################################################
def read_ini_file(ini_file_path):
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    config.optionxform = str
    # Read the INI file
    config.read(ini_file_path)
    return config

################################################################################
def write_ini_to_log(ini_file_path):
    # Create a ConfigParser object
    config = configparser.ConfigParser()
    config.optionxform = str
    # Read the INI file
    config.read(ini_file_path)

    # Convert the INI file content to a string
    ini_content = ""
    for section in config.sections():
        ini_content += f"[{section}]\n"
        for key, value in config.items(section):
            ini_content += f"{key} = {value}\n"
        ini_content += "\n"

    # Write the INI file content to the log file
    logging.info("INI file content:\n" + ini_content)
    
################################################################################
def run_command(command, output_file=None, append=False, add_header=False, header=None):
    try:
        # 记录即将执行的命令
        logging.info(f"Executing command: {' '.join(command)}")

        # 修改subprocess.run调用，捕获stdout和stderr
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

        # 记录标准输出
        if result.stdout:
            logging.error(f"Command stdout: {result.stdout}")

        # 记录标准错误（如果有的话）
        if result.stderr:
            logging.error(f"Command stderr: {result.stderr}")

        # 处理输出文件
        if output_file:
            mode = 'a' if append else 'w'
            with open(output_file, mode, encoding='utf-8') as f:
                if add_header and header and (not os.path.exists(output_file) or os.path.getsize(output_file) == 0):
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
    except KeyboardInterrupt:
        logging.info("User terminated the program using Ctrl+C.")

################################################################################
def validate_positive_integer(value, param_name):
    # Validate that the parameters are greater than zero and integers
    if not isinstance(value, int) or value <= 0:
        error_message = f"{param_name} must be a positive integer (>=1). Given value: {value}"
        # Log the error message
        logging.error(error_message)
        # Optionally, you can still raise the error to stop execution or signal that an error condition has occurred
        raise ValueError(error_message)
        sys.exit(1)

################################################################################
def check_file_format_efficient(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
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

##########################################################################################################################################
def check_fasta_sequence_length(fasta_file_path):  
    sequences = []  
    with open(fasta_file_path, 'r', encoding='utf-8') as fasta_file:  
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
        logging.error(f"There are multiple sequences in the file {fasta_file_path}.")
        sys.exit(1)  
    elif len(sequences) == 0:  
        logging.error(f"There is no sequence in the file {fasta_file_path}.")
        sys.exit(1)   
    else:  
        return len(sequences[0])  

################################################################################
def setup_logging(mode, include_file_handler=True):
    handlers = [logging.StreamHandler(sys.stdout)]
    if include_file_handler:
        handlers.append(logging.FileHandler("log.txt", mode=mode))
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=handlers)

################################################################################
def parse_log(value):  
    # 修改了isinstance函数的用法，需要两个参数：对象和类型  
    if isinstance(value, str):  
        normalized_value = value.strip().upper()  
        if normalized_value in ["Y", "YES"]:  
            return True  
        elif normalized_value in ["N", "NO"]:  
            return False  
        else:  
            # 记录错误日志  
            logging.error("Invalid log input for 'parse_log'. Valid input must be one of ['Y', 'YES', 'N', 'NO'].")  
            raise ValueError("Valid input must be one of ['Y', 'YES', 'N', 'NO']")  
    else:  
        # 记录错误日志  
        logging.error(f"Type error in parse_log: Input type is {type(value)}, but must be a string.")  
        raise TypeError("Input must be a string.")

################################################################################
def validate_file_format(file_path):
    expected_columns = 8
    expected_headers = [
        "fragment_id", "start", "end", "direction", 
        "paired_id", "paired_start", "paired_end", "paired_direction"
    ]

    # 在错误信息中提取文件名
    file_name = os.path.basename(file_path)

    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            # Read the header line
            header_line = file.readline().strip()
            headers = header_line.split('\t')
            
            # Check the number of columns
            if len(headers) != expected_columns:
                raise ValueError(f"Number of columns in the file {file_name} does not match the expected format. Expected: {expected_columns}, Found: {len(headers)}")
            
            # Check the column names
            if headers != expected_headers:
                raise ValueError(f"Column names in the file {file_name} do not match the expected format. Expected: {expected_headers}, Found: {headers}")
            
            # Check the data format for each row
            for line_number, line in enumerate(file, start=2):
                line = line.replace(' ', '').replace('\r', '').replace('\n', '')              # 这里保留/t，但去掉其他的空白字符
                if not line: 
                    continue  # Skip empty lines
                
                # 计算 /t 的个数
                tab_count = line.count('\t')

                # Check the number of columns based on tab count
                if tab_count != expected_columns - 1:
                    raise ValueError(f"Line {line_number} in the file {file_name} has an incorrect number of columns. Expected: {expected_columns}, Found: {tab_count + 1}")
                
                columns = line.split('\t')
                
                # Check if start, end, paired_start, and paired_end are integers
                try:
                    start = int(float(columns[1])) if columns[1] else None
                    end = int(float(columns[2])) if columns[2] else None
                    paired_start = int(float(columns[5])) if columns[5] else None
                    paired_end = int(float(columns[6])) if columns[6] else None
                except ValueError:
                    raise ValueError(f"Line {line_number} in the file {file_name} contains invalid values for start, end, paired_start, or paired_end. Expected integers.")
                
                # Check if direction and paired_direction are either "plus" or "minus"
                if columns[3] not in ["plus", "minus", ""] or columns[7] not in ["plus", "minus", ""]:
                    raise ValueError(f"Line {line_number} in the file {file_name} contains invalid values for direction or paired_direction. Expected: 'plus' or 'minus'")
    
    except FileNotFoundError:
        logging.error(f"File '{file_path}' not found.")
        sys.exit(1)
    except ValueError as e:
        logging.error(f"{e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Unknown error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("User terminated the program using Ctrl+C.")
        
        
################################################################################
def check_rp_ids_in_color_library(file_path, color_library):
    """
    Check if RU-prefixed repeat sequence IDs in the file exist in the color library.

    Args:
        file_path (str): Path to the input file.
        color_library (list): List of tuples containing RU IDs and their corresponding colors.

    Returns:
        None

    Raises:
        ValueError: If any RU ID is missing in the color library.
        FileNotFoundError: If the file does not exist.
    """
    # Convert color_library to a set of RU IDs for faster lookup
    color_library_ids = {rp_id for rp_id, _ in color_library}

    # Function to extract the base RU ID (e.g., 'RU61a' -> 'RU61')
    def extract_base_rp_id(fragment_id):
        if fragment_id.startswith('RU'):
            # Remove the suffix (e.g., 'a', 'b', 'aa', etc.)
            base_rp_id = ''.join([char for char in fragment_id if char.isdigit() or (char in ['R', 'P'])])
            return base_rp_id
        return None

    # Read the file and check for missing RU IDs
    missing_rp_ids = set()

    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            # Skip the header line
            next(file)
            
            # Process each line
            for line in file:
                columns = line.strip().split('\t')
                fragment_id = columns[0]
                
                # Extract the base RU ID
                base_rp_id = extract_base_rp_id(fragment_id)
                
                # Check if the base RU ID is in the color_library
                if base_rp_id and base_rp_id not in color_library_ids:
                    missing_rp_ids.add(base_rp_id)
        
        # Report missing RU IDs
        if missing_rp_ids:
            error_message = f"The following Repeat IDs are missing in the '[color_library]' section: {', '.join(missing_rp_ids)}"
            logging.error(error_message)
            print()
            sys.exit(1)
        else:
            logging.info("All Repeat IDs are present in the '[color_library]' section.")

    except FileNotFoundError:
        logging.error(f"File '{file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("User terminated the program using Ctrl+C.")

################################################################################
def config_has_required_sections(config_path: str) -> bool:
    """检查配置文件是否包含所有必需的节"""
    required_sections = {
        'general',
        'mainconfiguration', 
        'IR_mediated_reverse_recomb',
        'DR_mediated_recomb_1to2',
        'DR_mediated_recomb_2to1',
        'DR_mediated_recomb_2to2',
        'mapper_config',
        'color_library',
        'Arrange_map'
    }
    
    config = configparser.ConfigParser()
    config.read(config_path)
    
    # 检查所有必需的节是否存在
    return required_sections.issubset(config.sections())
    
##########################################################################################################################################
def main():
    parser = argparse.ArgumentParser(description="MiRIV: A tool to map the confgiure of your organelle genome.")
    parser.add_argument("-c", dest="config", help="Path to external configuration file.", required=True)
    parser.add_argument("-redo", help="Delete all previous results and start calculation anew.", action="store_true")
    parser.add_argument("-v", "--version", action="version", version="MiRIV 1.0", help="Show the version number and exit.")

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)
    except KeyboardInterrupt:
        logging.info("User terminated the program using Ctrl+C.")

    # 验证配置文件是否包含所有必需的节
    convert_file_to_utf8(args.config)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
    if not config_has_required_sections(args.config):
        logging.error(f"Config file '{args.config}' may not be required by this program.")
        sys.exit(1)
        
    # Set up logging to console only
    setup_logging(mode='w', include_file_handler=True)

    # Load external configuration file
    config = read_ini_file(args.config)
    
    # write *.ini content to "log.txt".
    write_ini_to_log(args.config)

################################################################################
    ######## [general] section
    project_id = config['general']['project_id']

    # 检验 project_id 是否有效。无效时，退出程序。logging
    if is_valid_project_id(project_id):
        logging.info(f"The 'project_id' {project_id} is valid.")
    else:
        logging.error("Invalid 'project_id'. The 'project_id' must consist of letters, numbers, and underscores, and must start with a letter.")
        sys.exit(1)

################################################################################
    try:
        # 提醒用户进行颜色方案的设置
        print()
        logging.info("ATTENTION: Ensure each repeat is colored in the '[color_library]' section of the '.ini' file.")
        logging.info("ATTENTION: Different repeats must have unique colors.")
        logging.info("ATTENTION: The program will resume execution after 10 seconds.")
        print()
        # 等待 10 秒
        time.sleep(10)
    except KeyboardInterrupt:
        logging.info("User terminated the program using Ctrl+C.")
    
################################################################################
    ######## [mainconfiguration] section
    auto_map_main = config['mainconfiguration'].get('auto_map', 'Y').upper()
    if auto_map_main in ["Y","YES"]:
        inputfile_main = config['mainconfiguration']['inputfile']
        convert_file_to_utf8(inputfile_main)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        
        genome_length_main = int(config['mainconfiguration']['genome_length'])
        genome_type_main = config['mainconfiguration'].get('genome_type', 'C').upper()
        if not inputfile_main:
            logging.error(f"Paired repeat infornation in [mainconfiguration] section was not provided!")
            sys.exit(1)
        if not genome_length_main:
            logging.error(f"The 'genome length' in [mainconfiguration] section was not provided!") 
            sys.exit(1)
        if not genome_type_main:
            logging.error(f"The 'genome type' in [mainconfiguration] section was not provided!")
            sys.exit(1)
    elif auto_map_main in ["N","NO"]:
        pass
    else:
        logging.error("'Auto_map' value in [mainconfiguration] section was invalid. Please provide one of 'Y/Yes/No/N'.")
        sys.exit(1)

################################################################################
    ######## [IR_mediated_reverse_recomb] section
    auto_map_inv = config['IR_mediated_reverse_recomb'].get('auto_map', 'Y').upper()
    output_dir_prefix_inv = config['IR_mediated_reverse_recomb'].get('output_directory_prefix', 'Inv_Rev')
    if auto_map_inv in ["Y","YES","M"]:
        inputfile_inv = config['IR_mediated_reverse_recomb']['inputfile']
        inputfasta_inv = config['IR_mediated_reverse_recomb']['inputfasta']
        convert_file_to_utf8(inputfile_inv)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        convert_file_to_utf8(inputfasta_inv)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        
        genome_type_inv = config['IR_mediated_reverse_recomb'].get('genome_type', 'C').upper()
        if not inputfile_inv:
            logging.error(f"Paired repeat infornation in [IR_mediated_reverse_recomb] section was not provided!")
            sys.exit(1)
        if not inputfasta_inv:
            logging.error(f"The 'genome sequence' in [IR_mediated_reverse_recomb] section was not provided!")
            sys.exit(1)
        if not genome_type_inv:
            logging.error(f"The 'genome_type' in [IR_mediated_reverse_recomb] section was not provided!")
            sys.exit(1)
        if genome_type_inv not in ["L", "C"]:
            logging.error(f"The 'genome_type' in [IR_mediated_reverse_recomb] section should be 'L/C'!")
            sys.exit(1)
        genome_length_inv = check_fasta_sequence_length(inputfasta_inv)
    elif auto_map_inv in ["N","NO"]:
        pass
    else:
        logging.error("The 'auto_map' value in the [IR_mediated_reverse_recomb] section should be one of 'Yes/Y/No/N/M'.")
        sys.exit(1)
        
################################################################################
    ######## [DR_mediated_recomb_1to2] section
    auto_map_dr_1to2 = config['DR_mediated_recomb_1to2'].get('auto_map', 'Y').upper()
    output_dir_prefix_dr_1to2 = config['DR_mediated_recomb_1to2'].get('output_directory_prefix', 'DR_1to2')
    if auto_map_dr_1to2 in ["Y","YES","M"]:
        inputfile_1to2 = config['DR_mediated_recomb_1to2']['inputfile']
        inputfasta_1to2 = config['DR_mediated_recomb_1to2']['inputfasta']
        convert_file_to_utf8(inputfile_1to2)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        convert_file_to_utf8(inputfasta_1to2)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        
        genome_type_1to2 = config['DR_mediated_recomb_1to2'].get('genome_type', 'C').upper()
        if not inputfile_1to2:
            logging.error(f"Paired repeat infornation in the [DR_mediated_recomb_1to2] section was not provided!")
            sys.exit(1)
        if not inputfasta_1to2:
            logging.error(f"The 'genome sequence' in the [DR_mediated_recomb_1to2] section was not provided!")
            sys.exit(1)
        if not genome_type_1to2:
            logging.error(f"The 'genome_type' in the [DR_mediated_recomb_1to1] section was not provided!")
            sys.exit(1)
        if genome_type_1to2 not in ["L", "C"]:
            logging.error(f"The 'genome_type' in the [DR_mediated_recomb_1to2] section should be 'L/C'!")
            sys.exit(1)
        genome_length_1to2 = check_fasta_sequence_length(inputfasta_1to2)
    elif auto_map_dr_1to2 in ["N","NO"]:
        pass
    else:
        logging.error("The 'auto_map' value in the [auto_map_dr_1to2] section should be one of 'Yes/Y/No/N/M'.")
        exit(1)

################################################################################
    ######## [DR_mediated_recomb_2to1] section
    auto_map_dr_2to1 = config['DR_mediated_recomb_2to1'].get('auto_map', 'N').upper()
    output_dir_prefix_dr_2to1 = config['DR_mediated_recomb_2to1'].get('output_directory_prefix', 'DR_2to1')
    if auto_map_dr_2to1 in ["Y","YES","M"]:
        comp_ch_2to1_log = config['DR_mediated_recomb_2to1'].get('flip_chain', 'Y').upper()
        chr1_file_2to1 = config['DR_mediated_recomb_2to1'].get('chr1_file', '')
        chr2_file_2to1 = config['DR_mediated_recomb_2to1'].get('chr2_file', '')
        chr1_fasta_2to1 = config['DR_mediated_recomb_2to1'].get('chr1_fasta', '')
        chr2_fasta_2to1 = config['DR_mediated_recomb_2to1'].get('chr2_fasta', '')
        chr1_type_2to1 = config['DR_mediated_recomb_2to1'].get('chr1_type', 'C').upper()
        chr2_type_2to1 = config['DR_mediated_recomb_2to1'].get('chr2_type', 'C').upper()
        convert_file_to_utf8(chr1_file_2to1)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        convert_file_to_utf8(chr2_file_2to1)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        convert_file_to_utf8(chr1_fasta_2to1)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        convert_file_to_utf8(chr2_fasta_2to1)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        
        if not chr1_file_2to1 or not chr2_file_2to1:
            logging.error(f"Paired repeat infornation in [DR_mediated_recomb_2to1] section was not provided!")
            sys.exit(1)
        if not chr1_fasta_2to1 or not chr2_fasta_2to1:
            logging.error(f"The 'chromsome sequences' in [DR_mediated_recomb_2to1] section were not provided completely!")
            sys.exit(1)
        if not chr1_type_2to1 or not chr2_type_2to1:
            logging.error(f"The 'chromsome type' in [DR_mediated_recomb_2to1] section were not provided completely!")
            sys.exit(1)
        if chr1_type_2to1 not in ["L", "C"]: 
            logging.error(f"The 'chr1_type' in [DR_mediated_recomb_2to1] section should be 'L/C'!")
            sys.exit(1)
        if chr2_type_2to1 != "C":
            logging.error(f"The 'chr2_type' in [DR_mediated_recomb_2to1] section must be 'C'!")
            sys.exit(1)
        if comp_ch_2to1_log not in ['Y','YES','N','NO']:
            logging.error(f"The 'flip_chain' parameter in the '[DR_mediated_recomb_2to1]' section should be one of '['Y','YES','N','NO']'.")
            sys.exit(1)
            
        chr1_len_2to1 = check_fasta_sequence_length(chr1_fasta_2to1)
        chr2_len_2to1 = check_fasta_sequence_length(chr2_fasta_2to1)  
        comp_ch_2to1_log = parse_log(comp_ch_2to1_log)
        
    elif auto_map_dr_2to1 in ["N","NO"]:
        pass
    else:
        logging.error("The 'auto_map' value in the [DR_mediated_recomb_2to1] section should be one of 'Yes/Y/No/N/M'.")
        sys.exit(1)

################################################################################
    ######## [DR_mediated_recomb_2to2] section
    auto_map_dr_2to2 = config['DR_mediated_recomb_2to2'].get('auto_map', 'N').upper()
    output_dir_prefix_dr_2to2 = config['DR_mediated_recomb_2to2'].get('output_directory_prefix', 'DR_2to2')
    if auto_map_dr_2to2 in ["Y","YES","M"]:
        comp_ch_2to2_log = config['DR_mediated_recomb_2to1'].get('flip_chain', 'Y').upper()
        chr1_file_2to2 = config['DR_mediated_recomb_2to2'].get('chr1_file', '')
        chr2_file_2to2 = config['DR_mediated_recomb_2to2'].get('chr2_file', '')
        chr1_fasta_2to2 = config['DR_mediated_recomb_2to2'].get('chr1_fasta', '')
        chr2_fasta_2to2 = config['DR_mediated_recomb_2to2'].get('chr2_fasta', '')
        convert_file_to_utf8(chr1_file_2to2)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        convert_file_to_utf8(chr2_file_2to2)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        convert_file_to_utf8(chr1_fasta_2to2)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        convert_file_to_utf8(chr2_fasta_2to2)         # 文本格式的转换，防止Win/Mac下的文档在linux下读取出问题，统一转换为utf-8
        
        if not chr1_file_2to2 or not chr2_file_2to2:
            logging.error(f"Paired repeats infornation in [DR_mediated_recomb_2to2] section was not provided!")
            sys.exit(1)
        if not chr1_fasta_2to2 or not chr2_fasta_2to2:
            logging.error(f"The 'chrosome sequences' in [DR_mediated_recomb_2to2] section were not provided completely!")
            sys.exit(1)
        if comp_ch_2to2_log not in ['Y','YES','N','NO']:
            logging.error(f"The 'flip_chain' parameter in the '[DR_mediated_recomb_2to2]' section should be one of '['Y','YES','N','NO']'.")
            sys.exit(1)
            
        chr1_len_2to2 = check_fasta_sequence_length(chr1_fasta_2to2)
        chr2_len_2to2 = check_fasta_sequence_length(chr2_fasta_2to2)  
        comp_ch_2to2_log = parse_log(comp_ch_2to2_log)
        
    elif auto_map_dr_2to2 in ["N","NO"]:
        pass
    else:
        logging.error("The 'auto_map' value in the [DR_mediated_recomb_2to2] section should be one of 'Yes/Y/No/N/M'.")
        sys.exit(1)
        
################################################################################
    #### [mapper_config] section
    picture_box = int(config['mapper_config'].get('picture_box', '450'))
    radius = int(config['mapper_config'].get('radius', '150'))
    arrow_radius = int(config['mapper_config'].get('arrow_radius', '170'))
    arrow_size = int(config['mapper_config'].get('arrow_size', '10'))
    arrow_thickness = int(config['mapper_config'].get('arrow_thickness', '2'))
    font_size = int(config['mapper_config'].get('font_size', '18'))
    tag_height = int(config['mapper_config'].get('tag_height', '20'))
    tag_line_width = int(config['mapper_config'].get('tag_line_width', '1'))
    
    
    #### [color_library] section
    # 获取color_library节的所有键值对
    color_library = config.items('color_library')
    # 动态设置颜色变量
    color_dict = {}

    # 定义一个正则表达式模式来匹配以字母开头，字母+数字的形式
    pattern = re.compile(r'^[A-Za-z]+[0-9]+$')

    for i, (key, value) in enumerate(color_library, start=1):
        if i > 30:
            raise ValueError("No more than 30 kinds of repeats were supported to draw subconfiguration.")
            sys.exit(1)
        if not value:
            print(f"The color of {key} was not set.")
            sys.exit(1)
        if key in color_dict:
            raise ValueError("There are Repeats with duplicate names.")
            sys.exit(1)
        color_dict[key] = config['color_library'].get(key, value)
        
        # 检查key的格式
        if not pattern.match(key):
            raise ValueError(f"The format of Repeat name '{unique_fragment_ids}' is invalid.")
            logging.info("It must start with letters and be followed by one or more numbers (e.g., 'A1', 'Ab2').")
            sys.exit(1)
        
    # 添加额外的颜色直到总数达到30
    n = 0
    while len(color_dict) < 30:
        n += 1
        color_key = f"RU{n}"
        if color_key in color_dict:
            continue  # 如果color_key已经存在，则跳过
        color_dict[color_key] = 'white'
        
    # 使用字典存储颜色和键，以数字为键
    color_dict_select = {}
    for i, (key, value) in enumerate(color_library, start=1):
        # 确保value是小写的，并将键值对添加到字典中
        color_dict_select[i] = (key, value)
    
    # 通过下标访问元素 
    RU_color = {}
    RU_key = {}
    for index in color_dict_select: 
        key, color = color_dict_select[index]

        RU_key[index] = key
        RU_color[index] = color

################################################################################
    #### [Arrange_map] section
    arrange_map = config['Arrange_map']
    arrange = arrange_map.get('arrange', 'NO').upper()
    output_dir_prefix_arr = config['Arrange_map'].get('output_directory_prefix', 'map_nine_squares')

    # Extract the values
    if arrange in ["YES", "Y"]:
        font_size_arr = arrange_map.get('font_size', '20')
        image_dpi_arr = arrange_map.get('image_dpi', 600)
    elif arrange in ["NO", "N"]:
        pass
    else:
        logging.error("The 'arrange' value in [Arrange_map] section should be 'Yes/Y/No/N'.")
        sys.exit(1)
        
##########################################################################################################################################
######主程序开始，处理与前一次运行结果的关系
    logging.info("Current work is to MAP genome configuration based on user provided repeat info!\n")
    time.sleep(2)
    
    current_dir = os.path.dirname(os.path.realpath(__file__))

    ckecklists = [f"mainconfig_{project_id}", 
                  f"mainconfig_map.tsv",
                  f"{output_dir_prefix_inv}_{project_id}",
                  f"{output_dir_prefix_dr_1to2}_{project_id}",
                  f"{output_dir_prefix_dr_2to1}_{project_id}",
                  f"{output_dir_prefix_arr}_{project_id}",
                  f"fig_layout_{project_id}.tsv",
                  f"{output_dir_prefix_dr_2to2}_{project_id}",
                  f"{project_id}"
                  ]
    for ckecklist in ckecklists:
        if os.path.exists(ckecklist) and not args.redo:
            logging.warning(f"Warning: Previous result(s) {ckecklist} are found.")
            logging.info("ATTENTION: Choose '-redo' to delete previous results and re-calculate!")
            logging.info(f"ATTENTION: Or reset the 'project_id' in {args.config} file!")
            sys.exit(1)
        if os.path.exists(ckecklist) and args.redo:
            if os.path.isdir(ckecklist):
                logging.info(f"Delete previous folder: {ckecklist}.")
                shutil.rmtree(ckecklist)
                print()
            elif os.path.isfile(ckecklist):
                logging.info(f"Delete previous file: {ckecklist}.")
                os.remove(ckecklist)

##########################################################################################################################################
    time.sleep(2)
    
    # 提前处理数据，生成5CT
    if auto_map_main in ["Y","YES"]:
        logging.info(f"%%%%%%%%%%%% Map the mainconfiguration! %%%%%%%%%%%%")
        logging.info(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        validate_file_format(inputfile_main)            # 检验读入的8CT是否符合程序要求的格式
        check_rp_ids_in_color_library(inputfile_main, color_library)           # 检验重复序列有没有被设置颜色
        
        prefix = os.path.splitext(os.path.basename(inputfile_main))[0]
        os.system(f"python {current_dir}/paired_info_to_5CT.py -i {inputfile_main} -o {prefix}_{project_id}_5CT.tsv -l {genome_length_main}")
        
        if not os.path.exists(f"{project_id}"):    # 创建结果存储文件夹，以project_id命名
            os.mkdir(f"{project_id}")

################################################################################
    #### 绘制主构型图谱
    if auto_map_main in ["Y","YES"]:
        if os.path.exists(inputfile_main) and genome_type_main == "C":
            os.system(f"python {current_dir}/map_recomb.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_main} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os mainconfig_{project_id}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
        elif os.path.exists(inputfile_main) and genome_type_main == "L":
            os.system(f"python {current_dir}/map_recomb_line.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_main} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os mainconfig_{project_id}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
        else:
            logging.error(f"File {inputfile_main} does not exist")
            sys.exit(1)

        # 移动产生的主要构型图谱至文件夹mainconfig内
        if not os.path.exists(f"mainconfig_{project_id}"):
            os.mkdir(f"mainconfig_{project_id}")
            shutil.move(f"mainconfig_{project_id}_map.svg", f"mainconfig_{project_id}")
        else:
            shutil.move(f"mainconfig_{project_id}_map.svg", f"mainconfig_{project_id}")
            os.remove(f"{prefix}_{project_id}_5CT.tsv")

        df = pd.read_csv(inputfile_main, sep='\t').iloc[:, :8]
    
        # 定义输出路径和文件名
        project_id = f"{project_id}"  # 替换为您的项目 ID
        output_path = f"mainconfig_{project_id}/"
        output_file = f"mainconfig_{project_id}_map.tsv"
    
        # 保存提取的数据到新的 TSV 文件
        df.to_csv(f"{output_path}{output_file}", sep='\t', index=False)
        
        if os.path.exists(f"mainconfig_{project_id}"):        #  结果重新归档
            shutil.move(f"mainconfig_{project_id}", f"{project_id}")
        
        logging.info(f"The mainconfiguration map has been completed. The results have been saved in folder '{project_id}/mainconfig_{project_id}'!")
        print()
    
##########################################################################################################################################
    # 提前处理数据，生成5CT
    if auto_map_inv in ["Y","YES","M"]:
        logging.info(f"%%%%%%%%%%%% Map the IR-mediated subconfiguration! %%%%%%%%%%%%")
        logging.info(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        validate_file_format(inputfile_inv)       # 检验读入的8CT是否符合程序要求的格式
        check_rp_ids_in_color_library(inputfile_inv, color_library)           # 检验重复序列有没有被设置颜色
        
        prefix = os.path.splitext(os.path.basename(inputfile_inv))[0]
        os.system(f"python {current_dir}/paired_info_to_5CT.py -i {inputfile_inv} -o {prefix}_{project_id}_5CT.tsv -l {genome_length_inv}")
        
        if not os.path.exists(f"{project_id}"):     # 创建结果存储文件夹，以project_id命名
            os.mkdir(f"{project_id}")

################################################################################
    #### 绘制反转重组图谱，创建文件夹
    if auto_map_inv in ["Y","YES","M"]:
        IR_folder = f"{output_dir_prefix_inv}_{project_id}"
        if not os.path.exists(IR_folder) and auto_map_inv in ["Y","YES","M"]:
            os.mkdir(IR_folder)

    if (auto_map_inv == "YES" or auto_map_inv == "Y") and genome_type_inv == "C":
        # 8CT转换来的5CT，将其绘制为mainconfiguration的图谱
        os.system(f"python {current_dir}/map_recomb.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {IR_folder}/mainconfig_{prefix}_{project_id}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")

        os.system(f"python {current_dir}/IR_inv_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_inv} -f {inputfasta_inv} -t {genome_type_inv} -o {IR_folder}/{project_id} -auto")   #产生用于绘图的5CT，以及5CT对应的8CT
        os.remove(f"{prefix}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix}_{project_id}_5CT.tsv") else None    # 删除8CT转换来的5CT，删除中间结果
        
        IR_files = glob.glob(f"{IR_folder}/{project_id}*_5CT.tsv")
        for IR_file in IR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {IR_file} -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(IR_file))[0]}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            IR_file_basename = f"{os.path.splitext(os.path.basename(IR_file))[0]}.svg"  # 获取文件名（不包含路径）
            IR_file_without_ext, ext = os.path.splitext(IR_file_basename)  # 分离文件名和扩展名
            IR_file_new_basename = IR_file_without_ext.replace("_5CT", "_map") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{IR_file_basename}", f"{IR_file_new_basename}")
            shutil.move(f"{IR_file_new_basename}", f"{IR_folder}")
            os.remove(IR_file) if os.path.exists(IR_file) else None
        if os.path.exists(f"{IR_folder}"):        #  结果重新归档
            shutil.move(f"{IR_folder}", f"{project_id}")

    if (auto_map_inv == "YES" or auto_map_inv == "Y") and genome_type_inv == "L":
        # 8CT转换来的5CT，将其绘制为mainconfiguration的图谱
        os.system(f"python {current_dir}/map_recomb_line.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {IR_folder}/mainconfig_{prefix}_{project_id}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")

        os.system(f"python {current_dir}/IR_inv_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_inv} -f {inputfasta_inv} -t {genome_type_inv} -o {IR_folder}/{project_id} -auto")
        os.remove(f"{prefix}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix}_{project_id}_5CT.tsv") else None    # 删除8CT转换来的5CT，删除中间结果
        
        IR_files = glob.glob(f"{IR_folder}/{project_id}*_5CT.tsv")
        for IR_file in IR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {IR_file} -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(IR_file))[0]}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            IR_file_basename = f"{os.path.splitext(os.path.basename(IR_file))[0]}.svg"  # 获取文件名（不包含路径）
            IR_file_without_ext, ext = os.path.splitext(IR_file_basename)  # 分离文件名和扩展名
            IR_file_new_basename = IR_file_without_ext.replace("_5CT", "_map") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{IR_file_basename}", f"{IR_file_new_basename}")
            shutil.move(f"{IR_file_new_basename}", f"{IR_folder}")
            os.remove(IR_file) if os.path.exists(IR_file) else None
        if os.path.exists(f"{IR_folder}"):        #  结果重新归档
            shutil.move(f"{IR_folder}", f"{project_id}")

    if auto_map_inv == "M" and genome_type_inv == "C":
        # 8CT转换来的5CT，将其绘制为mainconfiguration的图谱
        os.system(f"python {current_dir}/map_recomb.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {IR_folder}/mainconfig_{prefix}_{project_id}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")

        os.system(f"python {current_dir}/IR_inv_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_inv} -f {inputfasta_inv} -t {genome_type_inv} -o {IR_folder}/{project_id}")
        os.remove(f"{prefix}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix}_{project_id}_5CT.tsv")  else None     # 删除8CT转换来的5CT，删除中间结果
        
        IR_files = glob.glob(f"{IR_folder}/{project_id}*_5CT.tsv")
        for IR_file in IR_files:  #if os.path.exists(f"{IR_file}"):
            os.system(f"python {current_dir}/map_recomb.py -i {IR_file} -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(IR_file))[0]}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            IR_file_basename = f"{os.path.splitext(os.path.basename(IR_file))[0]}.svg"  # 获取文件名（不包含路径）
            IR_file_without_ext, ext = os.path.splitext(IR_file_basename)  # 分离文件名和扩展名
            IR_file_new_basename = IR_file_without_ext.replace("_5CT", "_map") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{IR_file_basename}", f"{IR_file_new_basename}")
            shutil.move(f"{IR_file_new_basename}", f"{IR_folder}")
            os.remove(IR_file) if os.path.exists(IR_file) else None
        if os.path.exists(f"{IR_folder}"):        #  结果重新归档
            shutil.move(f"{IR_folder}", f"{project_id}")
            
    if auto_map_inv == "M" and genome_type_inv == "L":
        # 8CT转换来的5CT，将其绘制为mainconfiguration的图谱
        os.system(f"python {current_dir}/map_recomb_line.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {IR_folder}/mainconfig_{prefix}_{project_id}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")

        os.system(f"python {current_dir}/IR_inv_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_inv} -f {inputfasta_inv} -t {genome_type_inv} -o {IR_folder}/{project_id}")
        os.remove(f"{prefix}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix}_{project_id}_5CT.tsv")  else None     # 删除8CT转换来的5CT，删除中间结果
        
        IR_files = glob.glob(f"{IR_folder}/{project_id}*_5CT.tsv")
        for IR_file in IR_files:      
            os.system(f"python {current_dir}/map_recomb_line.py -i {IR_file} -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(IR_file))[0]}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            IR_file_basename = f"{os.path.splitext(os.path.basename(IR_file))[0]}.svg"  # 获取文件名（不包含路径）
            IR_file_without_ext, ext = os.path.splitext(IR_file_basename)  # 分离文件名和扩展名
            IR_file_new_basename = IR_file_without_ext.replace("_5CT", "_map") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{IR_file_basename}", f"{IR_file_new_basename}")
            shutil.move(f"{IR_file_new_basename}", f"{IR_folder}")
            os.remove(IR_file) if os.path.exists(IR_file) else None
        if os.path.exists(f"{IR_folder}"):        #  结果重新归档
            shutil.move(f"{IR_folder}", f"{project_id}")
            
    if auto_map_inv in ["Y","YES","M"]:
        if os.path.exists(f"{prefix}_{project_id}_5CT.tsv"):
            os.remove(f"{prefix}_{project_id}_5CT.tsv")
        
        print(f"Results are saved in folder {project_id}/{IR_folder}.") if auto_map_inv in ['M', 'Y', 'YES'] and os.listdir(f"{project_id}/{IR_folder}") else None
        print(f"No results are saved in folder {project_id}/{IR_folder}.") if auto_map_inv in ['M', 'Y', 'YES'] and not os.listdir(f"{project_id}/{IR_folder}") else None
        print()
        time.sleep(2)
    
##########################################################################################################################################
    # 提前处理数据，生成5CT
    if auto_map_dr_1to2 in ["Y","YES","M"]:
        logging.info(f"%%%%%%%%%% Map the DR-mediated subconfiguration (1to2)! %%%%%%%%%%")
        logging.info(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

        validate_file_format(inputfile_1to2)       # 检验读入的8CT是否符合程序要求的格式
        check_rp_ids_in_color_library(inputfile_1to2, color_library)           # 检验重复序列有没有被设置颜色

        prefix = os.path.splitext(os.path.basename(inputfile_1to2))[0]
        os.system(f"python {current_dir}/paired_info_to_5CT.py -i {inputfile_1to2} -o {prefix}_{project_id}_5CT.tsv -l {genome_length_1to2}")
        
        if not os.path.exists(f"{project_id}"):    # 创建结果存储文件夹，以project_id命名
            os.mkdir(f"{project_id}")
    
################################################################################
    #### 绘制1to2重组图谱，创建文件夹
    if auto_map_dr_1to2 in ["Y","YES","M"]:
        DR_folder_1to2 = f"{output_dir_prefix_dr_1to2}_{project_id}"
        if not os.path.exists(DR_folder_1to2):
            os.mkdir(DR_folder_1to2)

    if (auto_map_dr_1to2 == "YES" or auto_map_dr_1to2 == "Y") and genome_type_1to2 == "C":
        # 8CT转换来的5CT，将其绘制为mainconfiguration的图谱
        os.system(f"python {current_dir}/map_recomb.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_1to2}/mainconfig_{prefix}_{project_id}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")

        os.system(f"python {current_dir}/DR_1to2_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_1to2} -f {inputfasta_1to2} -t {genome_type_1to2} -o {DR_folder_1to2}/{project_id} -auto")
        os.remove(f"{prefix}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix}_{project_id}_5CT.tsv")  else None     # 删除8CT转换来的5CT，删除中间结果
        
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        if os.path.exists(f"{DR_folder_1to2}"):        #  结果重新归档
            shutil.move(f"{DR_folder_1to2}", f"{project_id}")

    if (auto_map_dr_1to2 == "YES" or auto_map_dr_1to2 == "Y") and genome_type_1to2 == "L":
        # 8CT转换来的5CT，将其绘制为mainconfiguration的图谱
        os.system(f"python {current_dir}/map_recomb_line.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_1to2}/mainconfig_{prefix}_{project_id}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")

        os.system(f"python {current_dir}/DR_1to2_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_1to2} -f {inputfasta_1to2} -t {genome_type_1to2} -o {DR_folder_1to2}/{project_id} -auto")
        os.remove(f"{prefix}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix}_{project_id}_5CT.tsv")  else None     # 删除8CT转换来的5CT，删除中间结果
        
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_Chr1_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            os.remove(f"{DR_file}")
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_Chr2_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        if os.path.exists(f"{DR_folder_1to2}"):        #  结果重新归档
            shutil.move(f"{DR_folder_1to2}", f"{project_id}")

    if auto_map_dr_1to2 == "M" and genome_type_1to2 == "C":
        # 8CT转换来的5CT，将其绘制为mainconfiguration的图谱
        os.system(f"python {current_dir}/map_recomb.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_1to2}/mainconfig_{prefix}_{project_id}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")

        os.system(f"python {current_dir}/DR_1to2_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_1to2} -f {inputfasta_1to2} -t {genome_type_1to2} -o {DR_folder_1to2}/{project_id}")
        os.remove(f"{prefix}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix}_{project_id}_5CT.tsv")  else None       # 删除8CT转换来的5CT，删除中间结果
        
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        if os.path.exists(f"{DR_folder_1to2}"):        #  结果重新归档
            shutil.move(f"{DR_folder_1to2}", f"{project_id}")

    if auto_map_dr_1to2 == "M" and genome_type_1to2 == "L":
        # 8CT转换来的5CT，将其绘制为mainconfiguration的图谱
        os.system(f"python {current_dir}/map_recomb_line.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_1to2}/mainconfig_{prefix}_{project_id}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")

        os.system(f"python {current_dir}/DR_1to2_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_1to2} -f {inputfasta_1to2} -t {genome_type_1to2} -o {DR_folder_1to2}/{project_id}")
        os.remove(f"{prefix}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix}_{project_id}_5CT.tsv")  else None       # 删除8CT转换来的5CT，删除中间结果
        
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_Chr1_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            os.remove(f"{DR_file}")
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_Chr2_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        if os.path.exists(f"{DR_folder_1to2}"):        #  结果重新归档
            shutil.move(f"{DR_folder_1to2}", f"{project_id}")
            
    if auto_map_dr_1to2 in ["Y","YES","M"]:
        if os.path.exists(f"{prefix}_{project_id}_5CT.tsv"):
            os.remove(f"{prefix}_{project_id}_5CT.tsv")

        print(f"Results are saved in folder {project_id}/{DR_folder_1to2}.") if auto_map_dr_1to2 in ['M', 'Y', 'YES'] and os.listdir(f"{project_id}/{DR_folder_1to2}") else None
        print(f"No results are saved in folder {project_id}/{DR_folder_1to2}.") if auto_map_dr_1to2 in ['M', 'Y', 'YES'] and not os.listdir(f"{project_id}/{DR_folder_1to2}") else None
        print()
    time.sleep(2)
    
##########################################################################################################################################
    #### 绘制2to1重组图谱，创建文件夹
    if auto_map_dr_2to1 in ["Y","YES","M"]:
        logging.info(f"%%%%%%%%%% Map the DR-mediated subconfiguration (2to1)! %%%%%%%%%%")
        logging.info(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        DR_folder_2to1 = f"{output_dir_prefix_dr_2to1}_{project_id}"
        if not os.path.exists(DR_folder_2to1):
            os.mkdir(DR_folder_2to1)
            
        if not os.path.exists(f"{project_id}"):    # 创建结果存储文件夹，以project_id命名
            os.mkdir(f"{project_id}")

        validate_file_format(chr1_file_2to1)       # 检验读入的8CT是否符合程序要求的格式
        validate_file_format(chr2_file_2to1)       # 检验读入的8CT是否符合程序要求的格式
        check_rp_ids_in_color_library(chr1_file_2to1, color_library)           # 检验重复序列有没有被设置颜色
        check_rp_ids_in_color_library(chr2_file_2to1, color_library)           # 检验重复序列有没有被设置颜色

        prefix_chr1 = os.path.splitext(os.path.basename(chr1_file_2to1))[0]
        os.system(f"python {current_dir}/paired_info_to_5CT.py -i {chr1_file_2to1} -o {prefix_chr1}_{project_id}_5CT.tsv -l {chr1_len_2to1}")
        prefix_chr2 = os.path.splitext(os.path.basename(chr2_file_2to1))[0]
        os.system(f"python {current_dir}/paired_info_to_5CT.py -i {chr2_file_2to1} -o {prefix_chr2}_{project_id}_5CT.tsv -l {chr2_len_2to1}")

################################################################################
    if (auto_map_dr_2to1 == "YES" or auto_map_dr_2to1 == "Y") and (chr1_type_2to1 == "C" and chr2_type_2to1 == "C"):
        os.system(f"python {current_dir}/map_recomb.py -i {prefix_chr1}_{project_id}_5CT.tsv -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to1}/{prefix_chr1}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
         
        os.system(f"python {current_dir}/map_recomb.py -i {prefix_chr2}_{project_id}_5CT.tsv -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to1}/{prefix_chr2}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")

        os.system(f"python bin/DR_2to1_tsv.py -i {chr1_file_2to1} -c1 {chr1_type_2to1} -j {chr2_file_2to1} -c2 {chr2_type_2to1} -l {chr1_fasta_2to1} -s {chr2_fasta_2to1} -o {DR_folder_2to1}/{project_id} -auto -g {comp_ch_2to1_log}")
        os.remove(f"{prefix_chr1}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr1}_{project_id}_5CT.tsv")  else None       # 删除8CT转换来的5CT，删除中间结果
        os.remove(f"{prefix_chr2}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr2}_{project_id}_5CT.tsv")  else None       # 删除8CT转换来的5CT，删除中间结果
        
        DR_files = glob.glob(f"{DR_folder_2to1}/{project_id}_DR_*_2to1_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            DR_2to1_basename = f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg"  # 获取文件名（不包含路径）
            DR_2to1_without_ext, ext = os.path.splitext(DR_2to1_basename)  # 分离文件名和扩展名
            DR_2to1_new_basename = DR_2to1_without_ext.replace("_5CT", "") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{DR_2to1_basename}", f"{DR_2to1_new_basename}")
            shutil.move(f"{DR_2to1_new_basename}", f"{DR_folder_2to1}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        DR_2to1_tsv = glob.glob("DR_*_2to1_map.tsv")
        for DR_2to1 in DR_2to1_tsv:
            shutil.move(f"{DR_2to1}", f"{DR_folder_2to1}/{project_id}_{DR_2to1}")
        if os.path.exists(f"{DR_folder_2to1}"):        #  结果重新归档
            shutil.move(f"{DR_folder_2to1}", f"{project_id}")
            
    if (auto_map_dr_2to1 == "YES" or auto_map_dr_2to1 == "Y") and (chr1_type_2to1 == "L" and chr2_type_2to1 == "C"):
        os.system(f"python {current_dir}/map_recomb_line.py -i {prefix_chr1}_{project_id}_5CT.tsv -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to1}/{prefix_chr1}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
         
        os.system(f"python {current_dir}/map_recomb.py -i {prefix_chr2}_{project_id}_5CT.tsv -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to1}/{prefix_chr2}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
        
        os.system(f"python bin/DR_2to1_tsv.py -i {chr1_file_2to1} -c1 {chr1_type_2to1} -j {chr2_file_2to1} -c2 {chr2_type_2to1} -l {chr1_fasta_2to1} -s {chr2_fasta_2to1} -o {DR_folder_2to1}/{project_id} -auto -g {comp_ch_2to1_log}")
        os.remove(f"{prefix_chr1}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr1}_{project_id}_5CT.tsv")  else None       # 删除8CT转换来的5CT，删除中间结果
        os.remove(f"{prefix_chr2}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr2}_{project_id}_5CT.tsv")  else None        # 删除8CT转换来的5CT，删除中间结果
        
        DR_files = glob.glob(f"{DR_folder_2to1}/{project_id}_DR_*_2to1_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            DR_2to1_basename = f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg"  # 获取文件名（不包含路径）
            DR_2to1_without_ext, ext = os.path.splitext(DR_2to1_basename)  # 分离文件名和扩展名
            DR_2to1_new_basename = DR_2to1_without_ext.replace("_5CT", "") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{DR_2to1_basename}", f"{DR_2to1_new_basename}")
            shutil.move(f"{DR_2to1_new_basename}", f"{DR_folder_2to1}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        DR_2to1_tsv = glob.glob("DR_*_2to1_map.tsv")
        for DR_2to1 in DR_2to1_tsv:
            shutil.move(f"{DR_2to1}", f"{DR_folder_2to1}/{project_id}_{DR_2to1}")
        if os.path.exists(f"{DR_folder_2to1}"):        #  结果重新归档
            shutil.move(f"{DR_folder_2to1}", f"{project_id}")
            
    if auto_map_dr_2to1 == "M" and (chr1_type_2to1 == "C" and chr2_type_2to1 == "C"):
        os.system(f"python {current_dir}/map_recomb.py -i {prefix_chr1}_{project_id}_5CT.tsv -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to1}/{prefix_chr1}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
         
        os.system(f"python {current_dir}/map_recomb.py -i {prefix_chr2}_{project_id}_5CT.tsv -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to1}/{prefix_chr2}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
        
        os.system(f"python bin/DR_2to1_tsv.py -i {chr1_file_2to1} -c1 {chr1_type_2to1} -j {chr2_file_2to1} -c2 {chr2_type_2to1} -l {chr1_fasta_2to1} -s {chr2_fasta_2to1} -o {DR_folder_2to1}/{project_id} -g {comp_ch_2to1_log}")
        os.remove(f"{prefix_chr1}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr1}_{project_id}_5CT.tsv")  else None        # 删除8CT转换来的5CT，删除中间结果
        os.remove(f"{prefix_chr2}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr2}_{project_id}_5CT.tsv")  else None        # 删除8CT转换来的5CT，删除中间结果
        
        DR_files = glob.glob(f"{DR_folder_2to1}/{project_id}_DR_*_2to1_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            DR_2to1_basename = f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg"  # 获取文件名（不包含路径）
            DR_2to1_without_ext, ext = os.path.splitext(DR_2to1_basename)  # 分离文件名和扩展名
            DR_2to1_new_basename = DR_2to1_without_ext.replace("_5CT", "") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{DR_2to1_basename}", f"{DR_2to1_new_basename}")
            shutil.move(f"{DR_2to1_new_basename}", f"{DR_folder_2to1}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        DR_2to1_tsv = glob.glob("DR_*_2to1_map.tsv")
        for DR_2to1 in DR_2to1_tsv:
            shutil.move(f"{DR_2to1}", f"{DR_folder_2to1}/{project_id}_{DR_2to1}")
        if os.path.exists(f"{DR_folder_2to1}"):        #  结果重新归档
            shutil.move(f"{DR_folder_2to1}", f"{project_id}")
            
    if auto_map_dr_2to1 == "M" and (chr1_type_2to1 == "L" and chr2_type_2to1 == "C"):
        os.system(f"python {current_dir}/map_recomb_line.py -i {prefix_chr1}_{project_id}_5CT.tsv -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to1}/{prefix_chr1}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
         
        os.system(f"python {current_dir}/map_recomb.py -i {prefix_chr2}_{project_id}_5CT.tsv -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to1}/{prefix_chr2}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
    
        os.system(f"python bin/DR_2to1_tsv.py -i {chr1_file_2to1} -c1 {chr1_type_2to1} -j {chr2_file_2to1} -c2 {chr2_type_2to1} -l {chr1_fasta_2to1} -s {chr2_fasta_2to1} -o {DR_folder_2to1}/{project_id} -g {comp_ch_2to1_log}")
        os.remove(f"{prefix_chr1}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr1}_{project_id}_5CT.tsv")  else None      # 删除8CT转换来的5CT，删除中间结果
        os.remove(f"{prefix_chr2}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr2}_{project_id}_5CT.tsv")  else None      # 删除8CT转换来的5CT，删除中间结果
        
        DR_files = glob.glob(f"{DR_folder_2to1}/{project_id}_DR_*_2to1_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {chr1_len_2to1+chr2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            DR_2to1_basename = f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg"  # 获取文件名（不包含路径）
            DR_2to1_without_ext, ext = os.path.splitext(DR_2to1_basename)  # 分离文件名和扩展名
            DR_2to1_new_basename = DR_2to1_without_ext.replace("_5CT", "") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{DR_2to1_basename}", f"{DR_2to1_new_basename}")
            shutil.move(f"{DR_2to1_new_basename}", f"{DR_folder_2to1}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        DR_2to1_tsv = glob.glob("DR_*_2to1_map.tsv")
        for DR_2to1 in DR_2to1_tsv:
            shutil.move(f"{DR_2to1}", f"{DR_folder_2to1}/{project_id}_{DR_2to1}")
        if os.path.exists(f"{DR_folder_2to1}"):        #  结果重新归档
            shutil.move(f"{DR_folder_2to1}", f"{project_id}")
        
    if auto_map_dr_2to1 in ["Y","YES","M"]:
        print(f"Results are saved in folder {project_id}/{DR_folder_2to1}.") if auto_map_dr_2to1 in ['M', 'Y', 'YES'] and os.listdir(f"{project_id}/{DR_folder_2to1}") else None
        print(f"No results are saved in folder {project_id}/{DR_folder_2to1}.") if auto_map_dr_2to1 in ['M', 'Y', 'YES'] and not os.listdir(f"{project_id}/{DR_folder_2to1}") else None
        print()
        time.sleep(2)
    
##########################################################################################################################################
    #### 绘制2to2重组图谱，创建文件夹
    if auto_map_dr_2to2 in ["Y","YES","M"]:
        logging.info(f"%%%%%%%%%% Map the DR-mediated subconfiguration (2to2)! %%%%%%%%%%")
        logging.info(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        DR_folder_2to2 = f"{output_dir_prefix_dr_2to2}_{project_id}"
        if not os.path.exists(DR_folder_2to2):
            os.mkdir(DR_folder_2to2)
            
        if not os.path.exists(f"{project_id}"):    # 创建结果存储文件夹，以project_id命名
            os.mkdir(f"{project_id}")

        validate_file_format(chr1_file_2to2)       # 检验读入的8CT是否符合程序要求的格式
        validate_file_format(chr2_file_2to2)       # 检验读入的8CT是否符合程序要求的格式
        check_rp_ids_in_color_library(chr1_file_2to2, color_library)           # 检验重复序列有没有被设置颜色
        check_rp_ids_in_color_library(chr2_file_2to2, color_library)           # 检验重复序列有没有被设置颜色

        prefix_chr1 = os.path.splitext(os.path.basename(chr1_file_2to2))[0]
        os.system(f"python {current_dir}/paired_info_to_5CT.py -i {chr1_file_2to2} -o {prefix_chr1}_{project_id}_5CT.tsv -l {chr1_len_2to2}")
        prefix_chr2 = os.path.splitext(os.path.basename(chr2_file_2to2))[0]
        os.system(f"python {current_dir}/paired_info_to_5CT.py -i {chr2_file_2to2} -o {prefix_chr2}_{project_id}_5CT.tsv -l {chr2_len_2to2}")
        
################################################################################
    if auto_map_dr_2to2 == "YES" or auto_map_dr_2to2 == "Y":           # 这种情况下，两条染色体只能都是线性的
        os.system(f"python {current_dir}/map_recomb_line.py -i {prefix_chr1}_{project_id}_5CT.tsv -l {chr1_len_2to2+chr2_len_2to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to2}/{prefix_chr1}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
         
        os.system(f"python {current_dir}/map_recomb_line.py -i {prefix_chr2}_{project_id}_5CT.tsv -l {chr1_len_2to2+chr2_len_2to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to2}/{prefix_chr2}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
        
        # DR_2to2_tsv.py  释放绘图用的5CT，由对应的8CT产生
        os.system(f"python bin/DR_2to2_tsv.py -i {chr1_file_2to2} -j {chr2_file_2to2} -l {chr1_fasta_2to2} -s {chr2_fasta_2to2} -o {DR_folder_2to2}/{project_id} -auto -g {comp_ch_2to2_log}")
        os.remove(f"{prefix_chr1}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr1}_{project_id}_5CT.tsv")  else None       # 删除8CT转换来的5CT，删除中间结果
        os.remove(f"{prefix_chr2}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr2}_{project_id}_5CT.tsv")  else None       # 删除8CT转换来的5CT，删除中间结果
        
        DR_files = glob.glob(f"{DR_folder_2to2}/{project_id}_DR_*_2to2_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {chr1_len_2to2+chr2_len_2to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            DR_2to2_basename = f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg"  #  DR_2to2_basename为绘制的图片文件的名字   具体为 MiRIV_DR_RU3a_2to2_5CT_map.svg
            DR_2to2_without_ext, ext = os.path.splitext(DR_2to2_basename)    # 分离文件名和扩展名
            DR_2to2_new_basename = DR_2to2_without_ext.replace("_5CT", "") + ext    # 去掉 "_5CT" 并添加回扩展名   # 更换后缀 _5CT_map.svg 为 _map.svg， ext为“svg”
            os.rename(f"{DR_2to2_basename}", f"{DR_2to2_new_basename}")    # 前面是准备名字，仅仅是字符串的操作，os.rename是真正的对文件改名字
            shutil.move(f"{DR_2to2_new_basename}", f"{DR_folder_2to2}")    # 处理产生的svg图片文件的名字，并移动位置
            if os.path.exists(DR_file):    # 删除绘图用的5CT
                os.remove(f"{DR_file}")
        DR_2to2_tsv = glob.glob("{DR_folder_2to2}/{project_id}_DR_*_5CT.tsv")    # 移动8CT 至 对应的文件夹内
        for DR_2to2 in DR_2to2_tsv:
            os.remove(f"{DR_2to2}", f"{DR_folder_2to2}/{DR_folder_2to2}")
        if os.path.exists(f"{DR_folder_2to2}"):        #  结果重新归档
            shutil.move(f"{DR_folder_2to2}", f"{project_id}")
        
    if auto_map_dr_2to2 == "M": 
        os.system(f"python {current_dir}/map_recomb_line.py -i {prefix_chr1}_{project_id}_5CT.tsv -l {chr1_len_2to2+chr2_len_2to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to2}/{prefix_chr1}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
         
        os.system(f"python {current_dir}/map_recomb_line.py -i {prefix_chr2}_{project_id}_5CT.tsv -l {chr1_len_2to2+chr2_len_2to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {DR_folder_2to2}/{prefix_chr2}_{project_id}.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
    
        os.system(f"python bin/DR_2to2_tsv.py -i {chr1_file_2to2} -j {chr2_file_2to2} -l {chr1_fasta_2to2} -s {chr2_fasta_2to2} -o {DR_folder_2to2}/{project_id} -g {comp_ch_2to2_log} -g {comp_ch_2to2_log}")
        os.remove(f"{prefix_chr1}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr1}_{project_id}_5CT.tsv")  else None       # 删除8CT转换来的5CT，删除中间结果
        os.remove(f"{prefix_chr2}_{project_id}_5CT.tsv") if os.path.isfile(f"{prefix_chr2}_{project_id}_5CT.tsv")  else None       # 删除8CT转换来的5CT，删除中间结果
        
        DR_files = glob.glob(f"{DR_folder_2to2}/{project_id}_DR_*_2to2_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {chr1_len_2to2+chr2_len_2to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RU_key[1]}:{RU_color[1]}, {RU_key[2]}:{RU_color[2]}, {RU_key[3]}:{RU_color[3]}, {RU_key[4]}:{RU_color[4]}, {RU_key[5]}:{RU_color[5]}, {RU_key[6]}:{RU_color[6]}, {RU_key[7]}:{RU_color[7]}, {RU_key[8]}:{RU_color[8]}, {RU_key[9]}:{RU_color[9]}, {RU_key[10]}:{RU_color[10]}, {RU_key[11]}:{RU_color[11]}, {RU_key[12]}:{RU_color[12]}, {RU_key[13]}:{RU_color[13]}, {RU_key[14]}:{RU_color[14]}, {RU_key[15]}:{RU_color[15]}, {RU_key[16]}:{RU_color[16]}, {RU_key[17]}:{RU_color[17]}, {RU_key[18]}:{RU_color[18]}, {RU_key[19]}:{RU_color[19]}, {RU_key[20]}:{RU_color[20]}, {RU_key[21]}:{RU_color[21]}, {RU_key[22]}:{RU_color[22]}, {RU_key[23]}:{RU_color[23]}, {RU_key[24]}:{RU_color[24]}, {RU_key[25]}:{RU_color[25]}, {RU_key[26]}:{RU_color[26]}, {RU_key[27]}:{RU_color[27]}, {RU_key[28]}:{RU_color[28]}, {RU_key[29]}:{RU_color[29]}, {RU_key[30]}:{RU_color[30]}'")
            DR_2to2_basename = f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg"  # 获取文件名（不包含路径）
            DR_2to2_without_ext, ext = os.path.splitext(DR_2to2_basename)  # 分离文件名和扩展名
            DR_2to2_new_basename = DR_2to2_without_ext.replace("_5CT", "") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{DR_2to2_basename}", f"{DR_2to2_new_basename}")
            shutil.move(f"{DR_2to2_new_basename}", f"{DR_folder_2to2}")
            if os.path.exists(DR_file):    # 删除绘图用的5CT
                os.remove(f"{DR_file}")
        DR_2to2_tsv = glob.glob("{DR_folder_2to2}/{project_id}_DR_*_5CT.tsv")    # 移动8CT 至 对应的文件夹内
        for DR_2to2 in DR_2to2_tsv:
            os.remove(f"{DR_2to2}", f"{DR_folder_2to2}/{DR_folder_2to2}")
        if os.path.exists(f"{DR_folder_2to2}"):        #  结果重新归档
            shutil.move(f"{DR_folder_2to2}", f"{project_id}")

    if auto_map_dr_2to2 in ["Y","YES","M"]:
        print(f"Results are saved in folder {project_id}/{DR_folder_2to2}.") if auto_map_dr_2to2 in ['M', 'Y', 'YES'] and os.listdir(f"{project_id}/{DR_folder_2to2}") else None
        print(f"No results are saved in folder {project_id}/{DR_folder_2to2}.") if auto_map_dr_2to2 in ['M', 'Y', 'YES'] and not os.listdir(f"{project_id}/{DR_folder_2to2}") else None
        print()
        time.sleep(2)
    
##########################################################################################################################################
##########################################################################################################################################
    ####  绘制九宫格 
    # Extract the values
    if arrange in ["YES", "Y"]:
        logging.info(f"%%%%%%%%%%%% Arrange maps into a grid of nine squares! %%%%%%%%%%%")
        logging.info(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        #### 创建文件夹
        arr_folder = f"{output_dir_prefix_arr}_{project_id}"
        if not os.path.exists(arr_folder):
            os.mkdir(arr_folder)
            
        if not os.path.exists(f"{project_id}"):    # 创建结果存储文件夹，以project_id命名
            os.mkdir(f"{project_id}")

        # Create a list of tuples to store the data
        data = [
            ('center', arrange_map.get('center_font'), arrange_map.get('center_path')),
            ('left_middle', arrange_map.get('left_middle_font'), arrange_map.get('left_middle_path')),
            ('right_middle', arrange_map.get('right_middle_font'), arrange_map.get('right_middle_path')),
            ('top_middle', arrange_map.get('top_middle_font'), arrange_map.get('top_middle_path')),
            ('bottom_middle', arrange_map.get('bottom_middle_font'), arrange_map.get('bottom_middle_path')),
            ('top_left', arrange_map.get('top_left_font'), arrange_map.get('top_left_path')),
            ('top_right', arrange_map.get('top_right_font'), arrange_map.get('top_right_path')),
            ('bottom_left', arrange_map.get('bottom_left_font'), arrange_map.get('bottom_left_path')),
            ('bottom_right', arrange_map.get('bottom_right_font'), arrange_map.get('bottom_right_path'))]

        # Write the data to a TSV file
        with open(f'fig_layout_{project_id}.tsv', 'w', encoding='utf-8') as file:
            file.write('fig_id\tlabel\tfigpath\n')
            for row in data:
                # 过滤掉None值，将非None值转换为字符串后组成新的列表
                valid_row = [str(item) for item in row if item is not None]
                if valid_row:  # 如果过滤后还有元素，就进行写入操作
                    file.write('\t'.join(valid_row) + '\n')
    
        if arrange == "Y" or arrange == "YES":
            os.system(f"python bin/collage.py -i fig_layout_{project_id}.tsv -o {arr_folder}/nine_squares_{project_id} -fs {font_size_arr} -dpi {image_dpi_arr}")
        
        if os.path.exists(f"fig_layout_{project_id}.tsv"):
            os.remove(f"fig_layout_{project_id}.tsv")
            
        if os.path.exists(f"{arr_folder}") and os.path.exists(f"{project_id}"):        #  结果重新归档
            shutil.move(f"{arr_folder}", f"{project_id}")
        
        logging.info(f"The Arranged map has been completed, saved in folder {project_id}/map_nine_squares_{project_id}!")
        time.sleep(2)
    
    #归档log.txt
    shutil.move("log.txt", f"{project_id}") if os.path.isfile("log.txt") else None 

##########################################################################################################################################
if __name__ == "__main__":

    start_time = time.time()     # Capture the start time
    main()
    end_time = time.time()    # Capture the end time when execution is about to finish
    
    execution_time = end_time - start_time  # Calculate the difference to get execution time
    print()
    print(f"Total execution time: {round(execution_time,1)} seconds.\n")  # Log the execution time


