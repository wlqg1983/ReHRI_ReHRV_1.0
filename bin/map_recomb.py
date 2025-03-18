#!/usr/bin/env python

import subprocess
import sys
import argparse
import os
import pandas as pd

##################################################################################################################################
default_values = {'-pb':450.0, '--picture_box':450, '-r':150, '--radius':150, '-ar':170.0, '--arrow_radius':170, '-as':10, '--arrow_size':10, '-at':2, '--arrow_thickness':2.5, '-fs':18, '--font_size':18, '-th':20, '--tag_height':20, '-tl':1, '--tag_line_width': 1}
##################################################################################################################################

def separate_args(external_args):
    five_ct_args = []
    plasmidrender_args = []
    genome_length = None

    has_five_ct = False
    has_plasmidrender = False

    i = 0
    while i < len(external_args):
        arg = external_args[i]

        if arg in ['-i', '-c', '-cf']:
            five_ct_args.append(arg)
            has_five_ct = True
            if i + 1 < len(external_args):
                five_ct_args.append(external_args[i + 1])
                i += 1
        elif arg == '-l':
            if i + 1 < len(external_args):
                genome_length = float(external_args[i + 1])
                i += 1
        elif arg in ['--picture_box', '-pb', '--radius', '-r', '--arrow_radius', '-ar', '--arrow_size', '-as', '--arrow_thickness', '-at', '--font_size', '-fs', '--tag_height', '-th', '--tag_line_width', '-tl']:
            plasmidrender_args.append(arg)
            has_plasmidrender = True
            if i + 1 < len(external_args):
                plasmidrender_args.append(external_args[i + 1])
                i += 1
        elif arg in ['--output_svg_file', '-os', '--output_png_file', '-op', '-cp']:
            if i + 1 < len(external_args):
                plasmidrender_args.append(arg)
                plasmidrender_args.append(external_args[i + 1])
                i += 1
        i += 1

    if '-cp' in plasmidrender_args and has_plasmidrender:
        print("The '-cp' needs to be used separately.")
        sys.exit(1)

    if not has_five_ct:
        print("The '-i' option is required. Please provide an input file.")
        sys.exit(1)

    if has_five_ct and ('-c' not in five_ct_args and '-cf' not in five_ct_args):
        print("Either -c or -cf option is required. Please provide one of them.")
        sys.exit(1)

    if has_plasmidrender and ('-os' not in plasmidrender_args and '-op' not in plasmidrender_args):
        print("Either -os or -op option is required for plasmidrender. Please provide one or two of them.")
        sys.exit(1)

    if not genome_length:
        print("You have not provided the length of the genome. The data will be plotted as the main configuration of the genome.")

    return five_ct_args, plasmidrender_args, genome_length

##########################################################################################################################################
def adjust_plasmidrender_args(genome_length, plasmidrender_args, args):
    plasmidrender_args = replace_short_options_with_long(plasmidrender_args) 

    # Read the input TSV file into a DataFrame
    input_df = pd.read_csv(args.input_file, sep='\t')

    ratio = None
    if genome_length and not args.config_path:
        # Get the maximum value between the last values in the second and third columns
        subtype_length = max(input_df.iloc[-1]['start'], input_df.iloc[-1]['end'])
        
        # Calculate adjustment based on the ratio of lengths
        ratio = subtype_length / genome_length

        if ratio < 0.35:            # 为了避免绘制的图形过小，对下限做出限制
            ratio = 0.35
            
        # Adjust parameters
        plasmidrender_args_dict = {plasmidrender_args[i]: plasmidrender_args[i + 1] for i in range(0, len(plasmidrender_args), 2)}

        # 使用一个函数来简化逻辑
        def adjust_parameter(param_name, default_key, min_val=None, max_val=None):
            if hasattr(args, param_name) and getattr(args, param_name) is not None:
                value = getattr(args, param_name)
            else:
                value = default_values[default_key]
            adjusted_value = int(value) * ratio
            if min_val is not None:
                adjusted_value = max(min_val, adjusted_value)
            if max_val is not None:
                adjusted_value = min(max_val, adjusted_value)
            return str(adjusted_value)

        plasmidrender_args_dict['--picture_box'] = adjust_parameter('picture_box', '--picture_box')
        plasmidrender_args_dict['-r'] = adjust_parameter('radius', '-r')
        plasmidrender_args_dict['--arrow_radius'] = adjust_parameter('arrow_radius', '--arrow_radius')
        plasmidrender_args_dict['--arrow_size'] = adjust_parameter('arrow_size', '--arrow_size', min_val=5)
        plasmidrender_args_dict['--arrow_thickness'] = adjust_parameter('arrow_thickness', '--arrow_thickness', min_val=1.5)
        plasmidrender_args_dict['--font_size'] = adjust_parameter('font_size', '--font_size', min_val=10)
        plasmidrender_args_dict['--tag_height'] = adjust_parameter('tag_height', '--tag_height', min_val=10)
        plasmidrender_args_dict['--tag_line_width'] = adjust_parameter('tag_line_width', '--tag_line_width')
        
        plasmidrender_args = [item for pair in plasmidrender_args_dict.items() for item in pair]

    return plasmidrender_args

################################################################################
    if args.config_path and genome_length:
        # Get the maximum value between the last values in the second and third columns
        subtype_length = max(input_df.iloc[-1]['start'], input_df.iloc[-1]['end'])
        # Calculate adjustment based on the ratio of lengths
        ratio = subtype_length / genome_length

        if ratio < 0.35:            # 为了避免绘制的图形过小，对下限做出限制
            ratio = 0.35

        # Adjust parameters
        if args.config_path is not None:
            plasmidrender_args_dict = read_config_file(args.config_path)
     
        plasmidrender_args_dict['--picture_box'] = adjust_parameter('picture_box', '--picture_box')
        plasmidrender_args_dict['-r'] = adjust_parameter('radius', '-r')
        plasmidrender_args_dict['--arrow_radius'] = adjust_parameter('arrow_radius', '--arrow_radius')
        plasmidrender_args_dict['--arrow_size'] = adjust_parameter('arrow_size', '--arrow_size', min_val=5)
        plasmidrender_args_dict['--arrow_thickness'] = adjust_parameter('arrow_thickness', '--arrow_thickness', min_val=1.5)
        plasmidrender_args_dict['--font_size'] = adjust_parameter('font_size', '--font_size', min_val=10)
        plasmidrender_args_dict['--tag_height'] = adjust_parameter('tag_height', '--tag_height', min_val=10)
        plasmidrender_args_dict['--tag_line_width'] = adjust_parameter('tag_line_width', '--tag_line_width')
     
        plasmidrender_args = [item for pair in plasmidrender_args_dict.items() for item in pair]
     
        return plasmidrender_args

################################################################################
    if not genome_length and args.config_path:
        plasmidrender_args_dict = read_config_file(args.config_path)
        plasmidrender_args = [item for pair in plasmidrender_args_dict.items() for item in pair]

        return plasmidrender_args

##########################################################################################################################################

def read_config_file(config_path):
    config_dict = {}

    with open(config_path, 'r') as file:
        for line in file:
            # Remove leading and trailing whitespaces, and convert to lowercase
            line = line.strip().lower()

            # Skip empty lines and comments
            if not line or line.startswith('#'):
                continue

            # Split the line into key and value
            key, value = map(str.strip, line.split('='))

            # Remove invisible characters (e.g., '\t', '\n', '\r')
            key = ''.join(char for char in key if char.isprintable())
            value = ''.join(char for char in value if char.isprintable())

            # Map keys to corresponding argparse format
            key_mapping = {
                'picture_box': '--picture_box',
                'radius': '-r',
                'arrow_radius': '--arrow_radius',
                'arrow_size': '--arrow_size',
                'arrow_thickness': '--arrow_thickness',
                'font_size': '--font_size',
                'tag_height': '--tag_height',
                'tag_line_width': '--tag_line_width',
                'output_svg_file': '--output_svg_file',
                'output_png_file': '--output_png_file'
                # Add more mappings as needed
            }

            # Convert key to argparse format using the mapping
            arg_key = key_mapping.get(key, key)

            # If the value is missing or empty, use the default value
            if not value:
                value = str(default_values.get(arg_key, ''))

            # Add the key-value pair to the dictionary
            config_dict[arg_key] = value

    return config_dict
    
##################################################################################################################################

def replace_short_options_with_long(plasmidrender_args):
    """替换 plasmidrender_args 中的短选项为长选项"""
    # 映射字典，将短选项映射到对应的长选项
    short_to_long_mapping = {
        '-pb': '--picture_box',
        '-ar': '--arrow_radius',
        '-as': '--arrow_size',
        '-at': '--arrow_thickness',
        '-fs': '--font_size',
        '-l': '--genome_length',
        '--radius': '-r',
        '-th': '--tag_height',
        '-tl': '--tag_line_width'
    }

    # 复制一份 plasmidrender_args，避免在原列表上进行修改
    replaced_args = plasmidrender_args.copy()

    for i in range(len(replaced_args)):
        if replaced_args[i] in short_to_long_mapping:
            replaced_args[i] = short_to_long_mapping[replaced_args[i]]

    # 返回替换后的列表
    return replaced_args

##################################################################################################################################

def main():
    # Create command-line argument parser
    parser = argparse.ArgumentParser(description='A script for mapping repeats in organellar genome.')
    parser.add_argument('-i', dest='input_file', required=True, help='Input file (TSV format, Containing 5 columns).')
    parser.add_argument('-c', dest='user_color_input', help='User-defined color settings.')
    parser.add_argument('-cf', dest='custom_color_file', help='User-defined color settings file.')
    parser.add_argument('-gt', dest='genome_type', help='Point out the linear or circular configuration of input genome.')

    # Plasmidrender-related arguments
    parser.add_argument('-l', '--genome_length', dest='genome_length', type=float, help='Main genome length for adjusting parameters for drawing pictures.')
    parser.add_argument('-pb', '--picture_box', dest='picture_box', type=float, help='Picture box size (default: {})'.format(default_values['-pb']))
    parser.add_argument('-r', '--radius', dest='radius', type=float, help='Radius (default: {})'.format(default_values['-r']))
    parser.add_argument('-ar', '--arrow_radius', dest='arrow_radius', type=float, help='Arrow radius (default: {})'.format(default_values['-ar']))
    parser.add_argument('-as', '--arrow_size', dest='arrow_size', type=int, help='Arrow size (default: {})'.format(default_values['-as']))
    parser.add_argument('-at', '--arrow_thickness', dest='arrow_thickness', type=int, help='Arrow thickness (default: {})'.format(default_values['-at']))
    parser.add_argument('-fs', '--font_size', dest='font_size', type=int, help='Font size (default: {})'.format(default_values['-fs']))
    parser.add_argument('-th', '--tag_height', dest='tag_height', type=int, help='Tag height (default: {})'.format(default_values['-th']))
    parser.add_argument('-tl', '--tag_line_width', dest='tag_line_width', type=int, help='Tag line width (default: {})'.format(default_values['-tl']))
    parser.add_argument('-os', '--output_svg_file', dest='output_svg_file', help='Output file prefix.')
    parser.add_argument('-op', '--output_png_file', dest='output_png_file', help='Output file prefix.')
    parser.add_argument('-cp', '--config_path', help='User-defined map settings.')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)
    
    # Separate arguments for FiveCT_to_json, plasmidrender, and others
    five_ct_args, plasmidrender_args, genome_length = separate_args(sys.argv[1:])

    # Assuming both recan.py and FiveCT_to_json.py are in the same directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    five_ct_script = os.path.join(script_dir, 'FiveCT_to_json.py')


    # Call FiveCT_to_json.py with subprocess and its arguments
    subprocess.run(['python', five_ct_script] + five_ct_args)

    # Assuming the JSON file is generated with the same prefix as specified in FiveCT_to_json
    json_file = 'Media_5CT_to_json.json'

    # Adjust plasmidrender_args based on parameters
    plasmidrender_args = adjust_plasmidrender_args(genome_length, plasmidrender_args, args)

    subprocess.run(['plasmidrender', '-i', json_file] + plasmidrender_args)
    os.remove(json_file) if os.path.isfile(json_file) else None

##################################################################################################################################

if __name__ == "__main__":
    main()

