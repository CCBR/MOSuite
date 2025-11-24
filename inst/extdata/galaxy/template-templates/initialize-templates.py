#!/usr/bin/env python

import json
import glob
import pathlib


def main():
    with open(
        "inst/extdata/galaxy/nidap/NIDAPBulkTemplate_parameterTo_MOSuiteMapping.json",
        "r",
    ) as infile:
        mapping = json.load(infile)
    for filename in glob.glob("inst/extdata/galaxy/nidap/*.code-template.json"):
        code_template_file = pathlib.Path(filename)
        with open(code_template_file, "r") as infile:
            code_template = json.load(infile)
            template_name = code_template_file.name
            new_template = {
                "title": code_template["title"].replace(" [CCBR]", ""),
                "description": code_template["description"],
                "codeTemplate": code_template["codeTemplate"],
                "columns": [],
                "inputDatasets": [],
                "parameters": [],
            }
            for arg_type in ("columns", "inputDatasets", "parameters"):
                for param in code_template.get(arg_type, []):
                    param_name = param.get("key")
                    params_lst = next(
                        (
                            m
                            for k, m in mapping.get("template_mappings", {}).items()
                            if k == template_name
                        ),
                        {},
                    ).get("parameter_mappings", [])
                    new_param = next(
                        (p.get(param_name) for p in params_lst if param_name in p), None
                    )
                    if new_param:
                        param["key"] = new_param
                        new_template[arg_type].append(param)
                    else:
                        new_template[arg_type].append(param)
            with open(
                pathlib.Path(
                    f"inst/extdata/galaxy/template-templates/{template_name}".replace(
                        "_CCBR_", ""
                    )
                ),
                "w",
            ) as outfile:
                json.dump(new_template, outfile, indent=4)


if __name__ == "__main__":
    main()
