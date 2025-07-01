#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import json
import warnings
import re, os, sys

#######################################################################################################################################################################
# Filter out MatplotlibDeprecationWarning
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Filter out FutureWarning for DataFrame concatenation
warnings.filterwarnings("ignore", category=FutureWarning)

# 预定义颜色字典，将颜色名称映射到对应的十六进制值
color_dict = {
'black': '#000000',
'red': '#ff0000',
'orange': '#ffa500',
'yellow': '#ffff00',
'green': '#008000',
'cyan': '#00ffff',
'blue': '#0000ff',
'purple': '#800080',
'brown': '#a52a2a',
'gray': '#808080',
'darkslategray': '#2f4f4f',
'dimgray': '#696969',
'navy': '#000080',
'indigo': '#4b0082',
'darkgreen': '#006400',
'darkred': '#8b0000',
'firebrick': '#b22222',
'crimson': '#dc143c',
'chocolate': '#d2691e',
'olive': '#808000',
'yellowgreen': '#9acd32',
'lawngreen': '#7cfc00',
'limegreen': '#32cd32',
'greenyellow': '#adff2f',
'lightseagreen': '#20b2aa',
'seagreen': '#2e8b57',
'darkseagreen': '#8fbc8f',
'lightgreen': '#90ee90',
'forestgreen': '#228b22',
'darkcyan': '#008b8b',
'mediumturquoise': '#48d1cc',
'turquoise': '#40e0d0',
'aquamarine': '#7fffd4',
'mediumaquamarine': '#66cdaa',
'aqua': '#00ffff',
'deepskyblue': '#00bfff',
'skyblue': '#87ceeb',
'steelblue': '#4682b4',
'cadetblue': '#5f9ea0',
'royalblue': '#4169e1',
'mediumblue': '#0000cd',
'darkviolet': '#9400d3',
'plum': '#dda0dd',
'deeppink': '#ff1493',
'hotpink': '#ff69b4',
'pink': '#ffc0cb',
'palevioletred': '#db7093',
'mediumvioletred': '#c71585',
'coral': '#ff7f50',
'orangered': '#ff4500',
'darkorange': '#ff8c00',
'goldenrod': '#daa520',
'gold': '#ffd700',
'khaki': '#f0e68c',
'darkkhaki': '#bdb76b',
'wheat': '#f5deb3',
'lightgrey': '#d3d3d3',
'lightslategray': '#778899',
'slategray': '#708090',
'darkgray': '#a9a9a9',
}

#######################################################################################################################################################################
def validate_color_format(color_str):
    # Remove leading and trailing whitespaces from the color
    color = color_str.strip()

    # Check if the color is in hexadecimal format
    if re.match(r'^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$', color):
        return 'hexadecimal'  # Indicates that the input is in hexadecimal format

    # Check if the color matches a key in the color_dict
    elif color.lower() in color_dict:
        return 'color_word'  # Indicates that the input matches a color_dict key

    else:
        # Check RGB format (R.G.B)
        rgb_values = color.split('.')
        if len(rgb_values) == 3 and all(re.match(r'^\d+(\.\d+)?$', value) for value in rgb_values):
            return 'RGB'  # Indicates that the input is in RGB format

        # If none of the above conditions are met, the format is unknown or not recognized
        print("ERROR: Custom color is not recognised. Please check the spelling of words, the range of hexadecimal numbers, or the RGB format.")
        print(f"ATTENTION: Problematic color value: {color}")
        sys.exit(1)

    # If the function reaches this point, the color format is unknown
    return 'unknown'

#######################################################################################################################################################################
# Import the function definition
def rgb_to_hex(rgb):
    r, g, b = rgb
    return "#{:02x}{:02x}{:02x}".format(r, g, b)

#######################################################################################################################################################################
def process_data_user_color(data, genome_length, output_file, user_color_input):
    # Initialize an empty dictionary for custom color settings
    custom_color_settings = {}

    # Parse user color input into the custom_color_settings dictionary
    for color_setting in user_color_input.split(','):
        tag, color = color_setting.split(':')
        # Use regex to extract the prefix (e.g., RU1 from RU1az)
        tag_prefix = re.match(r'([A-Za-z]+[0-9]+)', tag.strip()).group(1)
        color_format = validate_color_format(color)
        if color_format == 'RGB':
            color = rgb_to_hex([int(c) for c in color.split('.')])
        custom_color_settings[tag_prefix] = color.strip()
        #print(custom_color_settings)

    # Create an intermediate table
    result = pd.DataFrame(columns=['fragment_id', 'length', 'type', 'paired_id', 'angle', 'color'])

    # Iterate through the rows of the DataFrame
    for index, row in data.iterrows():
        fragment_id = row['fragment_id']
        length = row['end'] - row['start']
        fragment_type = 'space' if row['type'] == 'inter_region' else 'tag'
        paired_id = row['paired_id']
        angle = (length / genome_length) * 360

        if fragment_type != 'space':
            # Extract the prefix of the fragment_id
            fragment_prefix = re.match(r'([A-Za-z]+[0-9]+)', fragment_id)
            if fragment_prefix:
                fragment_prefix = fragment_prefix.group(1)
                # Assign color if the prefix exists in custom_color_settings
                if fragment_prefix in custom_color_settings:
                    result_color = custom_color_settings[fragment_prefix]
                else:
                    print()
                    print(f"Error: The color for displaying Repeat {fragment_prefix} has not defined.")
                    print(f"Modification: Open the '.ini' file. Revise and Define the color of Repeat {fragment_prefix} in the '[color_library]' section.")
                    print("Attention: Avoid setting colors duplicately for the same Repeat.")
                    print()
                    sys.exit(1)
            else:
                result_color = '#ffffff'  # Set default color to white if no prefix found
        else:
            result_color = '#ffffff'  # Set space color to empty

        new_row = pd.DataFrame({'fragment_id': [fragment_id], 'length': [length], 'type': [fragment_type], 'paired_id': [paired_id], 'angle': [angle], 'color': [result_color]}, index=[index])
        result = pd.concat([result, new_row], ignore_index=True)

    # Save the intermediate table as a CSV file
    result.to_csv(output_file, sep='\t', index=False)

#######################################################################################################################################################################
def process_data_custom_color_file(data, genome_length, output_file, color_file_path):
    '''对每一个重复进行上色'''
    # Initialize an empty set for sorted tags
    sorted_tags = set()

    custom_color_settings = {}  # Store user-defined color settings for repeated sequence pairs
    inconsistent_tags = []  # Lists to track inconsistencies

    # Read color settings from the color_file_path
    with open(color_file_path, 'r', encoding='utf-8') as color_file:
        for line in color_file:
            line = line.strip()
            if not line:
                continue
                
            # Split the line into individual tag-color pairs separated by commas
            color_settings = line.split(',')
            
            for color_setting in color_settings:
                color_protocol = validate_color_format(color_setting)
                parts = color_setting.split(':')
                if len(parts) != 2:
                    print(f"Error: Custom color format in the file is incorrect. Use 'tag:color' format. Or check that the word describing the colour is spelled correctly.")
                    print(f"Example: R03b:red, R01a:255.0.0, R03a:#FF0000")
                    sys.exit(1)

                tag, color = parts[0].strip(), parts[1].strip()
                # Extract the prefix from the tag
                tag_prefix = re.match(r'([A-Za-z]+[0-9]+)', tag)
                print(tag_prefix)
                if tag_prefix:
                    tag_prefix = tag_prefix.group(1)
                    # Store the user-defined colors in the custom_color_settings dictionary
                    if color_protocol == "RGB":
                        custom_color_settings[tag_prefix] = rgb_to_hex([int(c) for c in color.split('.')])
                    elif color_protocol == "color_word" or color_protocol == "hexadecimal":
                        custom_color_settings[tag_prefix] = color

    # Create an intermediate table
    result = pd.DataFrame(columns=['fragment_id', 'length', 'type', 'paired_id', 'angle', 'color'])

    # Iterate through the rows of the DataFrame
    for index, row in data.iterrows():
        fragment_id = row['fragment_id']
        length = row['end'] - row['start']
        fragment_type = 'space' if row['type'] == 'inter_region' else 'tag'
        paired_id = row['paired_id']
        angle = (length / genome_length) * 360

        if fragment_type != 'space':
            # Extract the prefix of the fragment_id
            fragment_prefix = re.match(r'([A-Za-z]+[0-9]+)', fragment_id)
            if fragment_prefix:
                fragment_prefix = fragment_prefix.group(1)
                # Assign color if the prefix exists in custom_color_settings
                if fragment_prefix in custom_color_settings:
                    result_color = custom_color_settings[fragment_prefix]
                else:
                    print(f"Error: The color for displaying Repeat {fragment_prefix} has not defined.")
                    print(f"Open the '.ini' file. Revise and Define the color of Repeat {fragment_prefix} in the '[color_library]' section.")
                    print("Avoid setting colors duplicately for the same Repeat.")
                    sys.exit(1)
            else:
                result_color = '#ffffff'  # Set default color to white if no prefix found
        else:
            result_color = '#ffffff'  # Set space color to empty

        new_row = pd.DataFrame({'fragment_id': [fragment_id], 'length': [length], 'type': [fragment_type], 'paired_id': [paired_id], 'angle': [angle], 'color': [result_color]}, index=[index])
        result = pd.concat([result, new_row], ignore_index=True)

    # Save the intermediate table as a CSV file
    result.to_csv(output_file, sep='\t', index=False)

#####################################################################################################################################################################

def convert_media_to_json(input_file):
    # Load the 'Media' DataFrame
    media_df = pd.read_csv(input_file, delimiter='\t')
    
    # Adjust angles to be in the range [-6, 6]
    media_df['angle'] = media_df['angle'].apply(lambda angle: 5 if angle >= 0 and angle <= 5 else (-5 if angle < 0 and angle >= -5 else angle))
    
    # Calculate the total angle sum using absolute values
    total_angle = media_df['angle'].abs().sum()
    
    # Calculate the adjustment needed
    adjustment = 350 - total_angle       ###########################################################3# 此处经摸索发现，被减数越小，map末端月不容易overlap，先设置为350，实际应该是360
    
    # Find 'space' areas with lengths greater than 6 degrees
    large_space_areas = media_df[(media_df['type'] == 'space') & (media_df['angle'].abs() > 5)]

    # Calculate the total length of long 'space' areas
    total_space_length = large_space_areas['angle'].abs().sum()
    if total_space_length < 0.001:            ################## 当出现为0 的情况，设定 total_space_length = 1
        total_space_length = 0.001

    # Calculate the adjustment per degree
    adjustment_per_degree = adjustment / total_space_length

    # Update 'space' areas, including those with angles less than 6 degrees
    for index, row in large_space_areas.iterrows():
        space_length = abs(row['angle'])
        new_space_length = space_length + (space_length * adjustment_per_degree)
        media_df.at[index, 'angle'] = new_space_length if row['angle'] > 0 else -new_space_length

    # Create an empty list to store the JSON data
    json_data = []

    # Iterate through the rows of the DataFrame
    for index, row in media_df.iterrows():
        media_type = row['type']
        angle = row['angle']
        label = row['fragment_id']
        paired_id = row['paired_id']

        # Append JSON data based on 'type' and using the color from the 'color' column
        if abs(angle) <= 10:
            if media_type == 'space':
                if angle < 0:
                    json_data.append({
                        "type": "space",
                        "angle": str(0.5),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(abs(angle+1)),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "arrow",
                        "angle": str(angle+1),
                        "label": label,
                        "color": "#000000"  # Black color for "arrow"
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(0.5),
                        "label": "",
                        "color": "#000000"
                    })
                elif angle > 0:
                    json_data.append({
                        "type": "space",
                        "angle": str(0.5),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "arrow",
                        "angle": str(angle-1),
                        "label": label,
                        "color": "#000000"  # Black color for "arrow"
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(angle-1),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(0.5),
                        "label": "",
                        "color": "#000000"
                    })
            elif media_type == 'tag':
                # Get the color from the 'color' column and strip any surrounding double quotes
                tag_color = row['color'].strip('"') if not pd.isna(row['color']) else "#000000"
                if angle < 0:
                    json_data.append({
                        "type": "space",
                        "angle": str(0.5),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "tag",
                        "angle": str(abs(angle+1)),
                        "label": "",
                        "color": tag_color
                    })
                    json_data.append({
                        "type": "arrow",
                        "angle": str(angle+1),
                        "label": label,
                        "color": "#000000"  # Black color for "arrow"
                    })
                    json_data.append({
                        "type": "tag",
                        "angle": str(0.5),
                        "label": "",
                        "color": tag_color
                    })
                elif angle > 0:
                    json_data.append({
                        "type": "tag",
                        "angle": str(0.5),
                        "label": "",
                        "color": tag_color
                    })
                    json_data.append({
                        "type": "arrow",
                        "angle": str(angle-1),
                        "label": label,
                        "color": "#000000"  # Black color for "arrow"
                    })
                    json_data.append({
                        "type": "tag",
                        "angle": str(angle-1),
                        "label": "",
                        "color": tag_color
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(0.5),
                        "label": "",
                        "color": "#000000"
                    })
                    
                    
        if abs(angle) > 10:
            if media_type == 'space':
                if angle < 0:
                    json_data.append({
                        "type": "space",
                        "angle": str(abs((angle+10)/2)),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(10),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "arrow",
                        "angle": str(-10),
                        "label": label,
                        "color": "#000000"  # Black color for "arrow"
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(abs((angle+10)/2)),
                        "label": "",
                        "color": "#000000"
                    })
                    
                elif angle > 0:
                    json_data.append({
                        "type": "space",
                        "angle": str((angle-10)/2),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "arrow",
                        "angle": str(10),
                        "label": label,
                        "color": "#000000"  # Black color for "arrow"
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(10),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str((angle-10)/2),
                        "label": "",
                        "color": "#000000"
                    })
            elif media_type == 'tag':
                # Get the color from the 'color' column and strip any surrounding double quotes
                tag_color = row['color'].strip('"') if not pd.isna(row['color']) else "#000000"
                if angle < 0:
                    tag_angle = angle if angle > -180 else -180         # plasmidrender 可能存在一个bug，当tag的角度小于 -180 度时，tag 的位置就会出问题。做出限制，不让其小于 -180度
                    json_data.append({
                        "type": "tag",
                        "angle": str(abs(tag_angle)),
                        "label": "",
                        "color": tag_color
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str((angle+10)/2),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "arrow",
                        "angle": str(-10),
                        "label": label,
                        "color": "#000000"  # Black color for "arrow"
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(abs((angle+10)/2)),
                        "label": "",
                        "color": "#000000"
                    })

                elif angle > 0:
                    tag_angle = angle if angle < 180 else 180         # plasmidrender 可能存在一个bug，当tag的角度大于180度时，tag的位置就会出问题。做出限制，不让其大于 180度
                    json_data.append({
                        "type": "tag",
                        "angle": str(tag_angle),
                        "label": "",
                        "color": tag_color
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(-10-(angle-10)/2),
                        "label": "",
                        "color": "#000000"
                    })
                    json_data.append({
                        "type": "arrow",
                        "angle": str(10),
                        "label": label,
                        "color": "#000000"  # Black color for "arrow"
                    })
                    json_data.append({
                        "type": "space",
                        "angle": str(10+(angle-10)/2),
                        "label": "",
                        "color": "#000000"
                    })

    # Convert the list of dictionaries to a JSON string
    json_string = json.dumps(json_data, indent=4)
    
    output_file = "Media_5CT_to_json.json"
    # Save the JSON data to the output file
    with open(output_file, 'w') as json_file:
        json_file.write(json_string)


###################################################################################################################################################################
if __name__ == "__main__":
    # Create command-line argument parser
    parser = argparse.ArgumentParser(description='Process data and convert to JSON')
    parser.add_argument('-i', dest='input_file', required=True, help='Input file (TSV format)')
    parser.add_argument('-c', dest='user_color_input', help='User-defined color settings')
    parser.add_argument('-cf', dest='custom_color_file', help='User-defined color settings file')

    # Parse command-line arguments
    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)
    
    # Check if the -i is provided
    if not args.input_file:
        print("Error: The -i must be provided.")
        sys.exit(1)
    else:
        # Check if either -c or -cf is provided
        if not (args.user_color_input or args.custom_color_file):
            print("Error: Please provide either -c or -cf.")
            sys.exit(1)
        else:
            if args.user_color_input and args.custom_color_file:
                print("Error: Please provide either -c or -cf, but not both.")
                parser.print_help()
                sys.exit(1)

    # Read the raw data
    data = pd.read_csv(args.input_file, delimiter='\t')
    data.columns = data.columns.str.lower()
    data['type'] = data['type'].str.replace(' ', '').str.lower().str.replace(r'\([^)]*\)', '', regex=True)
    genome_length = abs(data['end'] - data['start']).sum()

    # Create ArgumentParser object
    parser = argparse.ArgumentParser(description="Your program description here")


    # Check if either -c or -cf flag is provided
    if args.user_color_input and args.custom_color_file:
        print("Error: Please provide either -c or -cf, but not both.")
        sys.exit(1)
    elif args.user_color_input:
        process_data_user_color(data, genome_length, 'MediaTableEnd-Start.tsv', args.user_color_input)
    elif args.custom_color_file:
        process_data_custom_color_file(data, genome_length, 'MediaTableEnd-Start.tsv', args.custom_color_file)
    #else:
    #    process_data_no_color(data, genome_length, '/tmp/MediaTableEnd-Start.tsv')

    convert_media_to_json('MediaTableEnd-Start.tsv')
    os.remove('MediaTableEnd-Start.tsv')

    # Reset warning filters to default
    #warnings.resetwarnings()

