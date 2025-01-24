#!/usr/bin/env python

import argparse
import csv
from PIL import Image
import cairosvg
from io import BytesIO
import sys
import subprocess
import xml.etree.ElementTree as ET

#######################################################################################################################
def get_svg_size(svg_file_path):
    # Parse the SVG file
    tree = ET.parse(svg_file_path)
    root = tree.getroot()

    # Get the width and height from the SVG root element
    width = root.attrib.get('width')
    height = root.attrib.get('height')

    return width, height

def svg_to_pil_with_cairosvg(svg_file_path, dpi=300):
    # Get the size of the SVG
    width, height = get_svg_size(svg_file_path)

    # Convert SVG to PNG with specified DPI and output dimensions
    png_data = cairosvg.svg2png(url=svg_file_path, dpi=dpi, output_width=4*round(float(width)), output_height=4*round(float(height)), write_to=None, background_color="white")

    # Load the PNG data into a PIL image
    pil_image = Image.open(BytesIO(png_data))

    return pil_image

#######################################################################################################################
def validate_config_data(config_data):
    valid_fig_ids = {"center", "left_middle", "right_middle", "top_middle", "bottom_middle", "top_left", "top_right", "bottom_left", "bottom_right"}
    if any(item["fig_id"].lower() not in valid_fig_ids for item in config_data):
        raise ValueError(f" The postion of each figure must be one of the predefined position names {valid_fig_ids}, and it is case-insensitive.")
    if len(config_data) > 9:
        raise ValueError("The total number of images must not exceed 9. If more than 9 images are needed, please call this program again.")

#######################################################################################################################
def create_grid_images(config_data, output_dpi, output_path_template):
    validate_config_data(config_data)

    # 找出最大的图片尺寸
    max_width, max_height = 0, 0
    for item in config_data:
        image = svg_to_pil_with_cairosvg(item["figpath"], output_dpi)
        max_width = max(max_width, image.size[0])
        max_height = max(max_height, image.size[1])

    grid_size = (3 * max_width, 3 * max_height)

    # 创建 images 字典时，确保键都是小写
    images = {item["fig_id"].lower(): svg_to_pil_with_cairosvg(item["figpath"], output_dpi) for item in config_data}
    
    positions = {
        "center": (1, 1),
        "left_middle": (0, 1),
        "right_middle": (2, 1),
        "top_middle": (1, 0),
        "bottom_middle": (1, 2),
        "top_left": (0, 0),
        "top_right": (2, 0),
        "bottom_left": (0, 2),
        "bottom_right": (2, 2)
    }

    grid_number = 0
    for i in range(0, len(config_data), 9):
        grid_image = Image.new('RGB', grid_size, color='white')

        for item in config_data[i:i + 9]:
            fig_id = item["fig_id"].lower()  # 将 fig_id 转换为小写
            image = images[fig_id]

            column, row = positions[fig_id]
            x = column * max_width + (max_width - image.size[0]) // 2
            y = row * max_height + (max_height - image.size[1]) // 2 

            grid_image.paste(image, (x, y))

        grid_image.save(output_path_template.format(grid_number), dpi=(output_dpi, output_dpi))
        grid_number += 1

#######################################################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create multiple grid images from SVG files.')
    parser.add_argument('-c', '--config', type=str, required=True, help='Path to the configuration file')
    parser.add_argument('-dpi', '--dpi', type=int, default=300, help='Output resolution in DPI (e.g., 1200)')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file path template (e.g., output_image_{}.png)')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    config_data = []
    with open(args.config, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            config_data.append({'fig_id': row['fig_id'], 'label': row['label'], 'figpath': row['figpath']})

    # 示例：确保输出文件路径包含 .png 文件扩展名
    output_file = args.output
    if not output_file.lower().endswith('.png'):
        output_file += '.png'  # 确保扩展名是 .png

    create_grid_images(config_data, args.dpi, output_file)


