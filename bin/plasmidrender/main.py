import json
from . import svg_draw
from . import option


#input_file=None,input_json=None,output_svg_file=None,output_png_file=None,config_path=None,picture_box=450.0,radius=150.0,plasmid_width=2.5,tag_height=20,tag_line_width=1.0,cut_line_length=20.0,cut_line_thickness=2.5,arrow_size=10.0,arrow_radius=170.0,arrow_thickness=2.0,font=None,font_size=16.0,rotation_angle=-90.0
def draw(input_file=None,input_json=None,output_svg_file=None,output_png_file=None,config_path=None,picture_box=450.0,radius=150.0,plasmid_width=2.5,tag_height=20,tag_line_width=1.0,cut_line_length=20.0,cut_line_thickness=2.5,arrow_size=10.0,arrow_radius=170.0,arrow_thickness=2.0,font=None,font_size=18.0,rotation_angle=-90.0):
    option.arguments.input_file=input_file
    option.arguments.input_json=input_json
    option.arguments.output_svg_file=output_svg_file
    option.arguments.output_png_file=output_png_file
    option.arguments.config_path=config_path
    option.arguments.picture_box=picture_box
    option.arguments.radius=radius
    option.arguments.plasmid_width=plasmid_width
    option.arguments.tag_height=tag_height
    option.arguments.tag_line_width=tag_line_width
    option.arguments.cut_line_length=cut_line_length
    option.arguments.cut_line_thickness=cut_line_thickness
    option.arguments.arrow_size=arrow_size
    option.arguments.arrow_radius=arrow_radius
    option.arguments.arrow_thickness=arrow_thickness
    option.arguments.font=font
    option.arguments.font_size=font_size
    option.arguments.rotation_angle=rotation_angle
    main()


def main():
    if option.arguments.config_path!=None:
        option.option_json(option.arguments.config_path)
    if option.arguments.input_file==None:
        gene_list=json.loads(option.arguments.input_json)
    else:
        with open( option.arguments.input_file,mode="r") as f:
            gene_list = json.load(f)
    svg_text=svg_draw.head_svg()
    angle=0.0
    id=0
    flag_before_item_is_tag=False
    
    #print(f"gene_list:{gene_list}")
    
    for gene_item in gene_list:
        if "font_color" in gene_item:
            font_color=gene_item["font_color"]
        else:
            font_color="black"
        if gene_item["type"]=="tag":
            if flag_before_item_is_tag:
                angle+=5
            svg_text+=svg_draw.annular_sector(angle,gene_item["angle"],gene_item["color"],gene_item["label"],id,font_color)
            
            flag_before_item_is_tag=True
            id+=1
            angle+=float(gene_item["angle"])
        elif gene_item["type"]=="line":
            svg_text+=svg_draw.point(angle,gene_item["color"],gene_item["label"],id,font_color)
            id+=1
            flag_before_item_is_tag=False
        elif gene_item["type"]=="arrow":
            svg_text+=svg_draw.arrow(angle,gene_item["angle"],gene_item["color"],gene_item["label"],id,font_color)
            id+=1
            flag_before_item_is_tag=False
        else:
            flag_before_item_is_tag=False
            angle+=float(gene_item["angle"])



    svg_text+='</g></svg>'
    

    if option.arguments.output_svg_file!=None:
        svg_draw.save_SVG(option.arguments.output_svg_file,svg_text)
        #print(f"svg_text:{svg_text}")

    if option.arguments.output_png_file!=None:
        svg_draw.save_png(option.arguments.output_png_file,svg_text)


def main_option():
    option.arguments = option.get_option()
    main()


if __name__ == '__main__':
    option.arguments = option.get_option()
    main()
