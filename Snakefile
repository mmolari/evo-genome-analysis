# configfile: "config/config.yaml"
import yaml
import re


def extract_record_names(fasta_file):
    records = []
    with open(fasta_file, "r") as f:
        for line in f.readlines():
            if line.startswith(">"):
                records.append(line[1:].strip())
    return records


def parse_run_config(yaml_file):
    with open(yaml_file, "r") as f:
        rc = yaml.safe_load(f)

    # remove extension from files

    # extract record name for each reference
    records = {}
    for ref_name in rc["pileups"]:
        ref_file = rc["input"] + "/references/" + ref_name + ".fa"
        records[ref_name] = extract_record_names(ref_file)

    rc["records"] = records

    return rc


# check that `run_config` is defined
if not "run_config" in config:
    raise RuntimeError(
        "the parameter run_config is not defined. Define it with --config run_config=myfile.yml"
    )


run_config = parse_run_config(config["run_config"])
print(run_config)

in_fld = run_config["input"]
out_fld = run_config["output"]
pileups = run_config["pileups"]
records = run_config["records"]


include: "rules/pileup.smk"


rule all:
    input:
        [
            expand(rules.map_reads.output, ref_id=ref, read_id=reads)
            for ref, reads in pileups.items()
        ],
