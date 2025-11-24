#!/usr/bin/env python

import json
import glob
import pathlib

def main():
    with open('inst/extdata/galaxy/nidap/NIDAPBulkTemplate_parameterTo_MOSuiteMapping.json', 'r') as infile:
        mapping = json.load(infile)
    code_template_files = glob.glob('inst/extdata/galaxy/nidap/*.code-template.json')

    code_template_file = pathlib.Path('inst/extdata/galaxy/nidap/Clean_Raw_Counts_CCBR_.code-template.json')
    with open(code_template_file, 'r') as infile:
        code_template = json.load(infile)
        template_name = code_template_file.name
        new_template = {"title": code_template['title'], "description": code_template['description'], "codeTemplate": code_template['codeTemplate'], "columns": [], "inputDatasets": [], "parameters": []}
        for arg_type in ('columns', 'inputDatasets', 'parameters'):
            for param in code_template.get(arg_type, []):
                param_name = param.get('key')
                new_param  = mapping.get('template_mappings', {}).get(template_name, {}).get('parameter_mappings', {}).get(param_name, None)
                if new_param:
                    param['key'] = new_param
                    new_template[arg_type].append(param)
        with open(pathlib.Path(f'inst/extdata/galaxy/template-templates/{template_name}.json'), 'w') as outfile:
            json.dump(new_template, outfile, indent=4)


if __name__ == "__main__":
    main()