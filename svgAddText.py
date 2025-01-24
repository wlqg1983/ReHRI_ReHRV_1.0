#!/usr/bin/env python

import argparse
import xml.etree.ElementTree as ET
import sys

###############################################################################################################################################################################
def get_svg_circle_center(svg_content):
    root = ET.fromstring(svg_content)
    ns = {'ns0': 'http://www.w3.org/2000/svg'}
    circle = root.find(".//ns0:circle[@id='circle1']", namespaces=ns)
    if circle is not None:
        cx = float(circle.attrib['cx'])
        cy = float(circle.attrib['cy'])
        return cx, cy
    else:
        raise ValueError("No circle with id 'circle1' found in the SVG.")

###############################################################################################################################################################################
def add_text_to_existing_svg(existing_svg_filename, text, output_svg_filename, font_size):
    with open(existing_svg_filename, 'r') as f:
        svg_content = f.read()

    cx, cy = get_svg_circle_center(svg_content)

    # 设置文本属性
    text_style = f"font-size:{font_size}px; font-family:Arial; fill:black;"
    
    # 在SVG中添加文本到圆心位置
    text_x = cx + 30  # 微调，x轴“+”使文字向右移动
    text_y = cy + 30  # 微调，y轴“+”使文字向上移动

    # 注册SVG命名空间
    svg_namespace = 'http://www.w3.org/2000/svg'
    ET.register_namespace('', svg_namespace)

    # 创建文本元素，确保使用正确的命名空间
    text_element = ET.Element('{http://www.w3.org/2000/svg}text', attrib={
        'x': str(text_x), 
        'y': str(text_y), 
        'style': text_style, 
        'text-anchor': 'middle', 
        'dominant-baseline': 'central'
    })
    text_element.text = text

    # 将现有SVG内容添加到新的ElementTree对象
    root = ET.fromstring(svg_content)
    root.append(text_element)

    # 保存修改后的SVG文件
    tree = ET.ElementTree(root)
    tree.write(output_svg_filename)

###############################################################################################################################################################################
def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='Add text to an existing SVG file.')
    parser.add_argument('-i', '--input', required=True, help='Input SVG filename')
    parser.add_argument('-o', '--output', required=True, help='Output SVG filename')
    parser.add_argument('-t', '--text', required=True, help='Text to be added')
    parser.add_argument('-fs', '--font_size', type=int, default=20, help='Font size for the added text')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as e:
        parser.print_help()
        sys.exit(1)

    # 调用函数并传入参数
    add_text_to_existing_svg(args.input, args.text, args.output, args.font_size)

###############################################################################################################################################################################
if __name__ == "__main__":
    main()


