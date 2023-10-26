import yaml
import pathlib
import numpy as np
import json
import sys


#  check that config file is provided
def check_configfile():
    try:
        msg = f"""
----------- loading config file -----------
--- run config:
{config["run_config"]}

--- pileup
{config["pileup"]}

--- plots
{config["plots"]}
-------------------------------------------

"""
        print(msg)
    except:
        raise Exception("config file not specified. Please specify with --configfile flag.")
    
    

check_configfile()

# make create log folder, required for cluster execution
pathlib.Path("log").mkdir(exist_ok=True)


def extract_record_names(fasta_file):
    """Given a fasta file returns the list of record names.
    Each name is split at the first space occurrence.
    Checks that there are no duplicates."""
    records = []
    with open(fasta_file, "r") as f:
        for line in f.readlines():
            if line.startswith(">"):
                rec_name = line[1:].strip()
                assert (
                    " " not in rec_name
                ), f"space character not allowed in record name: {rec_name} from {fasta_file}"
                records.append(rec_name)
    assert len(records) == len(
        np.unique(records)
    ), f"warning, duplicated record name in {fasta_file}"
    return records


def parse_pileup_config(rc):
    """parser for pileup part of the config file."""
    # extract fasta record names for each reference
    records = {}
    for ref_name in rc["pileups"]:
        ref_file = rc["input"] + "/references/" + ref_name + ".fa"
        records[ref_name] = extract_record_names(ref_file)
    rc["ref_records"] = records

    return rc


# extract pileup config options
run_config = parse_pileup_config(config["run_config"])

in_fld = run_config["input"].removesuffix("/")
out_fld = run_config["output"].removesuffix("/")
pileups = run_config["pileups"]
ref_records = run_config["ref_records"]

# print run options
print("----------- run configuration ------------")
print("input folder:", in_fld)
print("\noutput folder:", out_fld)
print("\npileups:", json.dumps(pileups, indent=2))
print("\nreference records:", json.dumps(ref_records, indent=2))
print("------------------------------------------")


include: "rules/pileup.smk"
include: "rules/plots.smk"


rule all:
    input:
        rules.pileup_all.input,
        rules.plot_all.input,
