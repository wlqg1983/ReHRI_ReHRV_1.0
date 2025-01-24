#!/usr/bin/env python

import sys, os, configparser, re
import subprocess, shutil, time
import time, argparse
import pandas as pd
import logging, glob

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
        #result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
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
            with open(output_file, mode) as f:
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

##########################################################################################################################################
def check_fasta_sequence_length(fasta_file_path):  
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
        
##########################################################################################################################################
def main():
    parser = argparse.ArgumentParser(description="OGRecomMapper: A tool to map the confgiure of your organelle genome.")
    parser.add_argument("-c", dest="config", help="Path to external configuration file.", required=True)
    parser.add_argument("-redo", help="Delete all previous results and start calculation anew.", action='store_true')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
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

    # 检验 project_id 是否有效。无效时，退出程序。
    if is_valid_project_id(project_id):
        logging.info(f"The 'project_id' {project_id} is valid.")
    else:
        logging.error("Invalid 'project_id'. The 'project_id' must consist of letters, numbers, and underscores, and must start with a letter.")
        sys.exit(1)

################################################################################
    ######## [mainconfiguration] section
    auto_map_main = config['mainconfiguration'].get('auto_map', 'Y').upper()
    if auto_map_main in ["Y","YES"]:
        inputfile_main = config['mainconfiguration']['inputfile']
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
        comp_ch_2to1_log = config['DR_mediated_recomb_2to1'].get('complementary_chain', 'Y').upper()
        chrom1_file_2to1 = config['DR_mediated_recomb_2to1'].get('chrom1_file', '')
        chrom2_file_2to1 = config['DR_mediated_recomb_2to1'].get('chrom2_file', '')
        chrom1_fasta_2to1 = config['DR_mediated_recomb_2to1'].get('chrom1_fasta', '')
        chrom2_fasta_2to1 = config['DR_mediated_recomb_2to1'].get('chrom2_fasta', '')
        chrom1_type_2to1 = config['DR_mediated_recomb_2to1'].get('chrom1_type', 'C').upper()
        chrom2_type_2to1 = config['DR_mediated_recomb_2to1'].get('chrom2_type', 'C').upper()
        if not chrom1_file_2to1 or not chrom2_file_2to1:
            logging.error(f"Paired repeat infornation in [DR_mediated_recomb_2to1] section was not provided!")
            sys.exit(1)
        if not chrom1_fasta_2to1 or not chrom2_fasta_2to1:
            logging.error(f"The 'chromsome sequences' in [DR_mediated_recomb_2to1] section were not provided completely!")
            sys.exit(1)
        if not chrom1_type_2to1 or not chrom2_type_2to1:
            logging.error(f"The 'chromsome type' in [DR_mediated_recomb_2to1] section were not provided completely!")
            sys.exit(1)
        if chrom1_type_2to1 not in ["L", "C"]: 
            logging.error(f"The 'chrom1_type' in [DR_mediated_recomb_2to1] section should be 'L/C'!")
            sys.exit(1)
        if chrom2_type_2to1 != "C":
            logging.error(f"The 'chrom2_type' in [DR_mediated_recomb_2to1] section must be 'C'!")
            sys.exit(1)
        if comp_ch_2to1_log not in ['Y','YES','N','NO']:
            logging.error(f"The 'complementary_chain' parameter in the '[DR_mediated_recomb_2to1]' section should be one of '['Y','YES','N','NO']'.")
            sys.exit(1)
            
        chrom1_len_2to1 = check_fasta_sequence_length(chrom1_fasta_2to1)
        chrom2_len_2to1 = check_fasta_sequence_length(chrom2_fasta_2to1)  
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
        comp_ch_2to2_log = config['DR_mediated_recomb_2to1'].get('complementary_chain', 'Y').upper()
        chrom1_file_2to2 = config['DR_mediated_recomb_2to2'].get('chrom1_file', '')
        chrom2_file_2to2 = config['DR_mediated_recomb_2to2'].get('chrom2_file', '')
        chrom1_fasta_2to2 = config['DR_mediated_recomb_2to2'].get('chrom1_fasta', '')
        chrom2_fasta_2to2 = config['DR_mediated_recomb_2to2'].get('chrom2_fasta', '')
        if not chrom1_file_2to2 or not chrom2_file_2to2:
            logging.error(f"Paired repeats infornation in [DR_mediated_recomb_2to2] section was not provided!")
            sys.exit(1)
        if not chrom1_fasta_2to2 or not chrom2_fasta_2to2:
            logging.error(f"The 'chrosome sequences' in [DR_mediated_recomb_2to2] section were not provided completely!")
            sys.exit(1)
        if comp_ch_2to2_log not in ['Y','YES','N','NO']:
            logging.error(f"The 'complementary_chain' parameter in the '[DR_mediated_recomb_2to2]' section should be one of '['Y','YES','N','NO']'.")
            sys.exit(1)
            
        chrom1_len_2to2 = check_fasta_sequence_length(chrom1_fasta_2to2)
        chrom2_len_2to2 = check_fasta_sequence_length(chrom2_fasta_2to2)  
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
    #print(f"color_library: {color_library}")
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
        color_key = f"RP{n}"
        if color_key in color_dict:
            continue  # 如果color_key已经存在，则跳过
        color_dict[color_key] = 'white'
        
    # 使用字典存储颜色和键，以数字为键
    color_dict_select = {}
    for i, (key, value) in enumerate(color_library, start=1):
        #print(f"{key}: {value}")
        # 确保value是小写的，并将键值对添加到字典中
        color_dict_select[i] = (key, value)
    
    # 输出存储的字典内容
    #print("Selected Color dictionary:", color_dict_select)
    
    # 通过下标访问元素 
    RP_color = {}
    RP_key = {}
    for index in color_dict_select: 
        key, color = color_dict_select[index]

        RP_key[index] = key
        RP_color[index] = color
    #print(f"Index {index} - Key: {RP_key}, Color: {RP_color}")
    #print(RP_key[1],RP_color[1])



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
                logging.info(f"Delete previous folder: {ckecklist}")
                shutil.rmtree(ckecklist)
                print()
            elif os.path.isfile(ckecklist):
                logging.info(f"Delete previous file: {ckecklist}")
                os.remove(ckecklist)

##########################################################################################################################################
    time.sleep(2)
    
    # 提前处理数据，生成5CT
    if auto_map_main in ["Y","YES"]:
        logging.info(f"#### Map the mainconfiguration! ####")
        prefix = os.path.splitext(os.path.basename(inputfile_main))[0]
        os.system(f"python {current_dir}/paired_info_to_5CT.py -i {inputfile_main} -o {prefix}_{project_id}_5CT.tsv -l {genome_length_main}")

        if not os.path.exists(f"{project_id}"):    # 创建结果存储文件夹，以project_id命名
            os.mkdir(f"{project_id}")

################################################################################
    #### 绘制主构型图谱
    if auto_map_main in ["Y","YES"]:
        if os.path.exists(inputfile_main) and genome_type_main == "C":
            os.system(f"python {current_dir}/map_recomb.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_main} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os mainconfig_{project_id}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
        elif os.path.exists(inputfile_main) and genome_type_main == "L":
            os.system(f"python {current_dir}/map_recomb_line.py -i {prefix}_{project_id}_5CT.tsv -l {genome_length_main} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os mainconfig_{project_id}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
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
        
        logging.info(f"The mainconfiguration map has been completed. The results have been saved in folder '{project_id}/mainconfig_{project_id}'!\n")
    
##########################################################################################################################################
    # 提前处理数据，生成5CT
    if auto_map_inv in ["Y","YES","M"]:
        logging.info(f"#### Map the IR-mediated subconfiguration! ####")
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
        os.system(f"python {current_dir}/IR_inv_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_inv} -f {inputfasta_inv} -t {genome_type_inv} -o {IR_folder}/{project_id} -auto")   #产生用于绘图的5CT，以及5CT对应的8CT
        os.remove(f"{prefix}_{project_id}_5CT.tsv")    #删除8CT转换来的5CT，删除中间结果
        IR_files = glob.glob(f"{IR_folder}/{project_id}*_5CT.tsv")
        for IR_file in IR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {IR_file} -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(IR_file))[0]}.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            IR_file_basename = f"{os.path.splitext(os.path.basename(IR_file))[0]}.svg"  # 获取文件名（不包含路径）
            IR_file_without_ext, ext = os.path.splitext(IR_file_basename)  # 分离文件名和扩展名
            IR_file_new_basename = IR_file_without_ext.replace("_5CT", "_map") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{IR_file_basename}", f"{IR_file_new_basename}")
            shutil.move(f"{IR_file_new_basename}", f"{IR_folder}")
            os.remove(f"{IR_file}")
        if os.path.exists(f"{IR_folder}"):        #  结果重新归档
            shutil.move(f"{IR_folder}", f"{project_id}")

    if (auto_map_inv == "YES" or auto_map_inv == "Y") and genome_type_inv == "L":
        os.system(f"python {current_dir}/IR_inv_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_inv} -f {inputfasta_inv} -t {genome_type_inv} -o {IR_folder}/{project_id} -auto")
        IR_files = glob.glob(f"{IR_folder}/{project_id}*_5CT.tsv")
        for IR_file in IR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {IR_file} -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(IR_file))[0]}.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            IR_file_basename = f"{os.path.splitext(os.path.basename(IR_file))[0]}.svg"  # 获取文件名（不包含路径）
            IR_file_without_ext, ext = os.path.splitext(IR_file_basename)  # 分离文件名和扩展名
            IR_file_new_basename = IR_file_without_ext.replace("_5CT", "_map") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{IR_file_basename}", f"{IR_file_new_basename}")
            shutil.move(f"{IR_file_new_basename}", f"{IR_folder}")
            os.remove(f"{IR_file}")
        if os.path.exists(f"{IR_folder}"):        #  结果重新归档
            shutil.move(f"{IR_folder}", f"{project_id}")

    if auto_map_inv == "M" and genome_type_inv == "C":
        os.system(f"python {current_dir}/IR_inv_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_inv} -f {inputfasta_inv} -t {genome_type_inv} -o {IR_folder}/{project_id}")
        IR_files = glob.glob(f"{IR_folder}/{project_id}*_5CT.tsv")
        for IR_file in IR_files:  #if os.path.exists(f"{IR_file}"):
            os.system(f"python {current_dir}/map_recomb.py -i {IR_file} -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(IR_file))[0]}.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            IR_file_basename = f"{os.path.splitext(os.path.basename(IR_file))[0]}.svg"  # 获取文件名（不包含路径）
            IR_file_without_ext, ext = os.path.splitext(IR_file_basename)  # 分离文件名和扩展名
            IR_file_new_basename = IR_file_without_ext.replace("_5CT", "_map") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{IR_file_basename}", f"{IR_file_new_basename}")
            shutil.move(f"{IR_file_new_basename}", f"{IR_folder}")
            os.remove(f"{IR_file}")
        if os.path.exists(f"{IR_folder}"):        #  结果重新归档
            shutil.move(f"{IR_folder}", f"{project_id}")
            
    if auto_map_inv == "M" and genome_type_inv == "L":
        os.system(f"python {current_dir}/IR_inv_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_inv} -f {inputfasta_inv} -t {genome_type_inv} -o {IR_folder}/{project_id}")
        IR_files = glob.glob(f"{IR_folder}/{project_id}*_5CT.tsv")
        for IR_file in IR_files:  #if os.path.exists(f"{IR_file}"):
            os.system(f"python {current_dir}/map_recomb_line.py -i {IR_file} -l {genome_length_inv} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(IR_file))[0]}.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            IR_file_basename = f"{os.path.splitext(os.path.basename(IR_file))[0]}.svg"  # 获取文件名（不包含路径）
            IR_file_without_ext, ext = os.path.splitext(IR_file_basename)  # 分离文件名和扩展名
            IR_file_new_basename = IR_file_without_ext.replace("_5CT", "_map") + ext  # 去掉 "_5CT" 并添加回扩展名
            os.rename(f"{IR_file_basename}", f"{IR_file_new_basename}")
            shutil.move(f"{IR_file_new_basename}", f"{IR_folder}")
            os.remove(f"{IR_file}")
        if os.path.exists(f"{IR_folder}"):        #  结果重新归档
            shutil.move(f"{IR_folder}", f"{project_id}")
            
    if auto_map_inv in ["Y","YES","M"]:
        if os.path.exists(f"{prefix}_{project_id}_5CT.tsv"):
            os.remove(f"{prefix}_{project_id}_5CT.tsv")
        
    print(f"Results are saved in folder {project_id}/{IR_folder}.\n") if auto_map_inv in ['M', 'Y', 'YES'] and os.listdir(f"{project_id}/{IR_folder}") else None
    print(f"No results are saved in folder {project_id}/{IR_folder}.\n") if auto_map_inv in ['M', 'Y', 'YES'] and not os.listdir(f"{project_id}/{IR_folder}") else None
    time.sleep(2)
    
##########################################################################################################################################
    # 提前处理数据，生成5CT
    if auto_map_dr_1to2 in ["Y","YES","M"]:
        logging.info(f"#### Map the DR-mediated subconfiguration (1to2)! ####")
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
        os.system(f"python {current_dir}/DR_1to2_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_1to2} -f {inputfasta_1to2} -t {genome_type_1to2} -o {DR_folder_1to2}/{project_id} -auto")
        os.remove(f"{prefix}_{project_id}_5CT.tsv")    #删除8CT转换来的5CT，删除中间结果
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        if os.path.exists(f"{DR_folder_1to2}"):        #  结果重新归档
            shutil.move(f"{DR_folder_1to2}", f"{project_id}")

    if (auto_map_dr_1to2 == "YES" or auto_map_dr_1to2 == "Y") and genome_type_1to2 == "L":
        os.system(f"python {current_dir}/DR_1to2_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_1to2} -f {inputfasta_1to2} -t {genome_type_1to2} -o {DR_folder_1to2}/{project_id} -auto")
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_Chr1_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            os.remove(f"{DR_file}")
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_Chr2_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        if os.path.exists(f"{DR_folder_1to2}"):        #  结果重新归档
            shutil.move(f"{DR_folder_1to2}", f"{project_id}")

    if auto_map_dr_1to2 == "M" and genome_type_1to2 == "C":
        os.system(f"python {current_dir}/DR_1to2_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_1to2} -f {inputfasta_1to2} -t {genome_type_1to2} -o {DR_folder_1to2}/{project_id}")
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        if os.path.exists(f"{DR_folder_1to2}"):        #  结果重新归档
            shutil.move(f"{DR_folder_1to2}", f"{project_id}")

    if auto_map_dr_1to2 == "M" and genome_type_1to2 == "L":
        os.system(f"python {current_dir}/DR_1to2_tsv.py -i {prefix}_{project_id}_5CT.tsv -j {inputfile_1to2} -f {inputfasta_1to2} -t {genome_type_1to2} -o {DR_folder_1to2}/{project_id}")
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_Chr1_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            os.remove(f"{DR_file}")
        DR_files = glob.glob(f"{DR_folder_1to2}/{project_id}*_Chr2_1to2.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {genome_length_1to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            shutil.move(f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg", f"{DR_folder_1to2}")
            os.system(f"python {current_dir}/FiveCT_to_paired_info.py -i {DR_file} -o {DR_folder_1to2}")
            if os.path.exists(DR_file):
                os.remove(f"{DR_file}")
        if os.path.exists(f"{DR_folder_1to2}"):        #  结果重新归档
            shutil.move(f"{DR_folder_1to2}", f"{project_id}")
            
    if auto_map_dr_1to2 in ["Y","YES","M"]:
        if os.path.exists(f"{prefix}_{project_id}_5CT.tsv"):
            os.remove(f"{prefix}_{project_id}_5CT.tsv")

    print(f"Results are saved in folder {project_id}/{DR_folder_1to2}.\n") if auto_map_dr_1to2 in ['M', 'Y', 'YES'] and os.listdir(f"{project_id}/{DR_folder_1to2}") else None
    print(f"No results are saved in folder {project_id}/{DR_folder_1to2}.\n") if auto_map_dr_1to2 in ['M', 'Y', 'YES'] and not os.listdir(f"{project_id}/{DR_folder_1to2}") else None
    time.sleep(2)
    
##########################################################################################################################################
    #### 绘制2to1重组图谱，创建文件夹
    if auto_map_dr_2to1 in ["Y","YES","M"]:
        logging.info(f"#### Map the DR-mediated subconfiguration (2to1)! ####")
        DR_folder_2to1 = f"{output_dir_prefix_dr_2to1}_{project_id}"
        if not os.path.exists(DR_folder_2to1):
            os.mkdir(DR_folder_2to1)
            
        if not os.path.exists(f"{project_id}"):    # 创建结果存储文件夹，以project_id命名
            os.mkdir(f"{project_id}")

    if (auto_map_dr_2to1 == "YES" or auto_map_dr_2to1 == "Y") and (chrom1_type_2to1 == "C" and chrom2_type_2to1 == "C"):
        os.system(f"python bin/DR_2to1_tsv.py -i {chrom1_file_2to1} -c1 {chrom1_type_2to1} -j {chrom2_file_2to1} -c2 {chrom2_type_2to1} -l {chrom1_fasta_2to1} -s {chrom2_fasta_2to1} -o {DR_folder_2to1}/{project_id} -auto -g {comp_ch_2to1_log}")
        DR_files = glob.glob(f"{DR_folder_2to1}/{project_id}_DR_*_2to1_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {chrom1_len_2to1+chrom2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
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
            

    if (auto_map_dr_2to1 == "YES" or auto_map_dr_2to1 == "Y") and (chrom1_type_2to1 == "L" and chrom2_type_2to1 == "C"):
        os.system(f"python bin/DR_2to1_tsv.py -i {chrom1_file_2to1} -c1 {chrom1_type_2to1} -j {chrom2_file_2to1} -c2 {chrom2_type_2to1} -l {chrom1_fasta_2to1} -s {chrom2_fasta_2to1} -o {DR_folder_2to1}/{project_id} -auto -g {comp_ch_2to1_log}")
        DR_files = glob.glob(f"{DR_folder_2to1}/{project_id}_DR_*_2to1_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {chrom1_len_2to1+chrom2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
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
            
            
    if auto_map_dr_2to1 == "M" and (chrom1_type_2to1 == "C" and chrom2_type_2to1 == "C"):
        os.system(f"python bin/DR_2to1_tsv.py -i {chrom1_file_2to1} -c1 {chrom1_type_2to1} -j {chrom2_file_2to1} -c2 {chrom2_type_2to1} -l {chrom1_fasta_2to1} -s {chrom2_fasta_2to1} -o {DR_folder_2to1}/{project_id} -g {comp_ch_2to1_log}")
        DR_files = glob.glob(f"{DR_folder_2to1}/{project_id}_DR_*_2to1_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb.py -i {DR_file} -l {chrom1_len_2to1+chrom2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
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
            
            
    if auto_map_dr_2to1 == "M" and (chrom1_type_2to1 == "L" and chrom2_type_2to1 == "C"):
        os.system(f"python bin/DR_2to1_tsv.py -i {chrom1_file_2to1} -c1 {chrom1_type_2to1} -j {chrom2_file_2to1} -c2 {chrom2_type_2to1} -l {chrom1_fasta_2to1} -s {chrom2_fasta_2to1} -o {DR_folder_2to1}/{project_id} -g {comp_ch_2to1_log}")
        DR_files = glob.glob(f"{DR_folder_2to1}/{project_id}_DR_*_2to1_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {chrom1_len_2to1+chrom2_len_2to1} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
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
        
    print(f"Results are saved in folder {project_id}/{DR_folder_2to1}.\n") if auto_map_dr_2to1 in ['M', 'Y', 'YES'] and os.listdir(f"{project_id}/{DR_folder_2to1}") else None
    print(f"No results are saved in folder {project_id}/{DR_folder_2to1}.\n") if auto_map_dr_2to1 in ['M', 'Y', 'YES'] and not os.listdir(f"{project_id}/{DR_folder_2to1}") else None
    time.sleep(2)
    
##########################################################################################################################################
    #### 绘制2to2重组图谱，创建文件夹
    if auto_map_dr_2to2 in ["Y","YES","M"]:
        logging.info(f"#### Map the DR-mediated subconfiguration (2to2)! ####")
        DR_folder_2to2 = f"{output_dir_prefix_dr_2to2}_{project_id}"
        if not os.path.exists(DR_folder_2to2):
            os.mkdir(DR_folder_2to2)
            
        if not os.path.exists(f"{project_id}"):    # 创建结果存储文件夹，以project_id命名
            os.mkdir(f"{project_id}")

    if auto_map_dr_2to2 == "YES" or auto_map_dr_2to2 == "Y":  
        os.system(f"python bin/DR_2to2_tsv.py -i {chrom1_file_2to2} -j {chrom2_file_2to2} -l {chrom1_fasta_2to2} -s {chrom2_fasta_2to2} -o {DR_folder_2to2}/{project_id} -auto -g {comp_ch_2to2_log}")
        # DR_2to2_tsv.py  释放绘图用的5CT，还有5CT对应的8CT
        
        DR_files = glob.glob(f"{DR_folder_2to2}/{project_id}_DR_*_2to2_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {chrom1_len_2to2+chrom2_len_2to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
            DR_2to2_basename = f"{os.path.splitext(os.path.basename(DR_file))[0]}_map.svg"  #  DR_2to2_basename为绘制的图片文件的名字   具体为 MiRIV_DR_RP3a_2to2_5CT_map.svg
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
        os.system(f"python bin/DR_2to2_tsv.py -i {chrom1_file_2to2} -j {chrom2_file_2to2} -l {chrom1_fasta_2to2} -s {chrom2_fasta_2to2} -o {DR_folder_2to2}/{project_id} -g {comp_ch_2to2_log} -g {comp_ch_2to2_log}")
        DR_files = glob.glob(f"{DR_folder_2to2}/{project_id}_DR_*_2to2_5CT.tsv")
        for DR_file in DR_files:
            os.system(f"python {current_dir}/map_recomb_line.py -i {DR_file} -l {chrom1_len_2to2+chrom2_len_2to2} -pb {picture_box} -r {radius} -ar {arrow_radius} -as {arrow_size} -at {arrow_thickness} -fs {font_size} -th {tag_height} -tl {tag_line_width} -os {os.path.splitext(os.path.basename(DR_file))[0]}_map.svg -c '{RP_key[1]}:{RP_color[1]}, {RP_key[2]}:{RP_color[2]}, {RP_key[3]}:{RP_color[3]}, {RP_key[4]}:{RP_color[4]}, {RP_key[5]}:{RP_color[5]}, {RP_key[6]}:{RP_color[6]}, {RP_key[7]}:{RP_color[7]}, {RP_key[8]}:{RP_color[8]}, {RP_key[9]}:{RP_color[9]}, {RP_key[10]}:{RP_color[10]}, {RP_key[11]}:{RP_color[11]}, {RP_key[12]}:{RP_color[12]}, {RP_key[13]}:{RP_color[13]}, {RP_key[14]}:{RP_color[14]}, {RP_key[15]}:{RP_color[15]}, {RP_key[16]}:{RP_color[16]}, {RP_key[17]}:{RP_color[17]}, {RP_key[18]}:{RP_color[18]}, {RP_key[19]}:{RP_color[19]}, {RP_key[20]}:{RP_color[20]}, {RP_key[21]}:{RP_color[21]}, {RP_key[22]}:{RP_color[22]}, {RP_key[23]}:{RP_color[23]}, {RP_key[24]}:{RP_color[24]}, {RP_key[25]}:{RP_color[25]}, {RP_key[26]}:{RP_color[26]}, {RP_key[27]}:{RP_color[27]}, {RP_key[28]}:{RP_color[28]}, {RP_key[29]}:{RP_color[29]}, {RP_key[30]}:{RP_color[30]}'")
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
        
    print(f"Results are saved in folder {project_id}/{DR_folder_2to2}.\n") if auto_map_dr_2to2 in ['M', 'Y', 'YES'] and os.listdir(f"{project_id}/{DR_folder_2to2}") else None
    print(f"No results are saved in folder {project_id}/{DR_folder_2to2}.\n") if auto_map_dr_2to2 in ['M', 'Y', 'YES'] and not os.listdir(f"{project_id}/{DR_folder_2to2}") else None
    time.sleep(2)
    
##########################################################################################################################################
##########################################################################################################################################
    ####  绘制九宫格 
    # Extract the values
    if arrange in ["YES", "Y"]:
        logging.info(f"#### Arrange maps into a grid of nine squares! ####")
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
        with open(f'fig_layout_{project_id}.tsv', 'w') as file:
            file.write('fig_id\tlabel\tfigpath\n')
            for row in data:
                # 过滤掉None值，将非None值转换为字符串后组成新的列表
                valid_row = [str(item) for item in row if item is not None]
                if valid_row:  # 如果过滤后还有元素，就进行写入操作
                    file.write('\t'.join(valid_row) + '\n')
    
        if arrange == "Y" or arrange == "YES":
            os.system(f"python bin/collage.py -i fig_layout_{project_id}.tsv -o {arr_folder}/nine_squares_{project_id} -fs {font_size_arr} -dpi {image_dpi_arr}")
        
        #if os.path.exists(f"fig_layout_{project_id}.tsv"):
        #    os.remove(f"fig_layout_{project_id}.tsv")
            
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
    print(f"\nTotal execution time: {round(execution_time,1)} seconds.\n")  # Log the execution time


