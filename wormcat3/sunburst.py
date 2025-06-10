"""
Create sunburst HTML file with Regulated Gene Set Data (RGS)
"""
import json
import os
import pandas as pd
from pathlib import Path
from wormcat3.sunburst_template import sunburst_template_front, sunburst_template_back


def _read_input_annotations(file_nm_in):
    """
    Read the RGS Data and create JSON based on Category 3 info
    """
    nodes_dict = {}
    try:
        df = pd.read_csv(file_nm_in)
        node1_list= []
        nodes_dict = {"name":"rgs", "children":node1_list}

        cat3 = df.groupby(["Category.3"]).count()

        cat3_dict={}
        for cat3_index, row in cat3.iterrows():
            count = int(row['Wormbase.ID'])
            cat3_dict[cat3_index] = count

        for key in cat3_dict:
            components = key.split(':')
            size = int(cat3_dict[key])
            if len(components) == 1:
                node1_list.append({"name": components[0].strip(), "size": size})
            elif len(components) == 2:
                node_list = __getChildrenFor(components[0].strip(),nodes_dict)
                node_list.append({"name": components[1].strip(), "size": size})
            else:
                node_list = __getChildrenFor2(components[0].strip(),components[1].strip(),nodes_dict)
                node_list.append({"name": components[2].strip(), "size": size})
    except Exception as e:
        print("Error in _read_input_annotations:", e)
        pass
    return nodes_dict

def __getChildrenFor(parent, nodes_dict):
    """
    Utility function getChildrenFor given parent
    """
    children = nodes_dict['children']
    node_list = None
    for key in children:
        if parent == key['name']:
            if 'children' in key:
                node_list = key['children']
            break

    if node_list is None:
        node_list = []
        children.append({"name":parent, "children":node_list})
    return node_list


def __getChildrenFor2(grand_parent, parent, nodes_dict):
    """
    Utility function getChildrenFor given grand parent
    """
    children = __getChildrenFor(grand_parent, nodes_dict)
    node_list = None
    for key in children:
        if parent == key['name']:
            if 'children' in key:
                node_list = key['children']
            break

    if node_list is None:
        node_list = []
        children.append({"name":parent, "children":node_list})

    return node_list


def create_sunburst(dir_path: str, run_number: str) -> None:
    """
    Create a sunburst HTML file from the data in 'input_annotated_###.csv' file.
    """
    dir_path = Path(dir_path)
    rgs_file = dir_path / f"input_annotated_{run_number}.csv"
    html_file = dir_path / f"sunburst_{run_number}.html"

    # Load data and convert to JSON
    data = _read_input_annotations(rgs_file)
    json_data = json.dumps(data)

    # Insert JSON data into the template and write to file
    var_json_data = f"var json_data = {json_data}"
    rendered_html = sunburst_template_front + var_json_data + sunburst_template_back
    html_file.write_text(rendered_html, encoding="utf-8")



if __name__ == "__main__":
    print(f"starting {os.getcwd()}")
    create_sunburst('RGS_Feb-14-2020-11_45_54')


