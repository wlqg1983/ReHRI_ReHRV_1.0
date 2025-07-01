import os
import sys
import argparse
import pandas as pd
import subprocess
import shutil

def main():
    parser = argparse.ArgumentParser(description='Process SVG files based on a configuration file.')
    parser.add_argument('-i', '--input', required=True, help='Path to the configuration file')
    parser.add_argument('-o', '--output', required=True, help='Output folder for grid images')
    parser.add_argument('-dpi', '--dpi', type=int, default=300, help='Output resolution in DPI (e.g., 600)')
    parser.add_argument('-fs', '--font_size', type=int, default=20, help='Font size for the added text')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # Check if input and output arguments are provided
    if not (args.input and args.output):
        sys.exit(1)

    # Read the configuration file into a DataFrame
    df = pd.read_csv(args.input, delimiter='\t', comment='#')

    # Check for missing values in the 'figpath' column
    if df['figpath'].isna().any():
        print("Warning: Some rows have missing 'figpath' values. These rows will be skipped.")
        df = df.dropna(subset=['figpath'])  # Drop rows with missing 'figpath'

    # Create the output directory if it doesn't exist
    output_directory = 'svgAddText'
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
    
    fon_size = args.font_size

    # Iterate through rows in the DataFrame
    for index, row in df.iterrows():
        fig_id = row['fig_id'].lower()
        fig_label = row['label']
        fig_path = row['figpath']

        # Construct the command for svgAddText.py
        input_svg = fig_path
        text_to_add = fig_label
        output_svg = os.path.join(output_directory, os.path.basename(input_svg))
        command_svgAddText = [
            "python", os.path.join(os.path.dirname(__file__), "svgAddText.py"),
            "-i", str(input_svg),
            "-o", str(output_svg),
            "-t", str(text_to_add),
            "-fs", str(fon_size)
        ]

        # Use subprocess to call svgAddText.py
        subprocess.run(command_svgAddText)
        
    # Modify the third column (figpath)
    df['figpath'] = 'svgAddText/' + df['figpath'].apply(lambda x: os.path.basename(x))

    # Save the modified configuration file
    modified_config_file = os.path.join(output_directory, os.path.basename(args.input))
    df.to_csv(modified_config_file, sep='\t', index=False)

    # Construct the command for form_grid.py
    input_file_form_grid = modified_config_file
    output_folder_form_grid = args.output
    command_form_grid = [
        "python", os.path.join(os.path.dirname(__file__), "form_grid.py"),
        "-c", str(input_file_form_grid),
        "-o", str(output_folder_form_grid),
        "-dpi", str(args.dpi)]
    
    # Use subprocess to call form_grid.py
    subprocess.run(command_form_grid)
    
    # Delete the svgAddText directory
    shutil.rmtree(output_directory)

if __name__ == "__main__":
    main()
    