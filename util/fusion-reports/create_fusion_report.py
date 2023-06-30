#!/usr/bin/env python3
# encoding: utf-8

import os, sys, re
import json

from igv_reports import datauri


def create_fusion_report(template, fusions, output_filename, input_file_prefix):

    basedir = os.path.dirname(fusions)
    data_uris = {}

    ## read in template
    with open(template, "r") as f:
        data = f.readlines()

    ## determine where in the html template that the igv code should be inserted.
    for i, line in enumerate(data):
        j = line.find('<!-- start igv report here -->')
        if j >= 0:
            space = ' ' * j
            report_start = i + 1
            break

    else:
        print("file must contain the line \"<!-- start igv report here -->\"")
        return


    input_lines = data[report_start:]



    output_lines = []
    #insert json for fusion selection table
    with open(fusions, "r") as f:
        j = json.loads(f.read())
        flattend_json = json.dumps(j)

    output_lines.append('var tableJson = ' + flattend_json)


    ## make data uri's for each of the referenced input filenames in the html doc
    for line_index, line in enumerate(input_lines):
        
        is_index = line.find('indexURL:') > 0
        if is_index:
            # ignore the index files
            continue

        line = re.sub("__PREFIX__", input_file_prefix, line)
                
        m = re.search("url:\s*([\"\']([^\"\']+)[\"\'])", line, flags=re.IGNORECASE)
        if (m):
            core_match = m.group(1)
            filename = m.group(2)
            (line, count) = re.subn(core_match, "data[\"{}\"]".format(filename), line, count=1)
            if count != 1:
                raise(RuntimeError("couldnt perform replacement at: {}".format(line)))

                        
            if os.path.exists(os.path.join(basedir, filename)):
                data_uris[filename] = datauri.file_to_data_uri(os.path.join(basedir, filename))
            else:
                sys.stderr.write("Warning - not locating file: {}\n".format(os.path.join(basedir, filename)))
                
        output_lines.append(line)


    report_header =  data[:report_start]

    report_data_uris = create_data_var(data_uris, space)

    report_body = output_lines

    new_html_data = report_header + report_data_uris + report_body

    with open(output_filename, 'w') as f:
        f.writelines(new_html_data)


def create_data_var(data_uris, space=''):
    data = []
    for i, (key, value) in enumerate(data_uris.items()):
        data.append('{}"{}": "{}"{}\n'.format(space + ' ' * 4, key, value, ',' if i < len(data_uris) - 1 else ''))
    return [space + "var data = {\n"] + data + [space+"};\n"]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--html_template", help="the html file to be converted", required=True, type=str)
    parser.add_argument("--fusions_json", help="json file defining the fusions (fusion inspector output)", required=True, type=str)
    parser.add_argument("--html_output", help="filename for html output", required=True, type=str)
    parser.add_argument("--input_file_prefix", help="prefix for input files", required=True, type=str)
    
    args = parser.parse_args()
    
    create_fusion_report(args.html_template, args.fusions_json, args.html_output, args.input_file_prefix)
    
    sys.exit(0)
