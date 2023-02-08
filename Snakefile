configfile: "config.yml"


import yaml
import pathlib
import numpy as np
import json

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


def parse_pileup_config(yaml_file):
    """parser for yaml pileup config file."""

    # read file
    with open(yaml_file, "r") as f:
        rc = yaml.safe_load(f)

    # extract fasta record names for each reference
    records = {}
    for ref_name in rc["pileups"]:
        ref_file = rc["input"] + "/references/" + ref_name + ".fa"
        records[ref_name] = extract_record_names(ref_file)
    rc["ref_records"] = records

    return rc


# check that `run_config` is defined
if not "run_config" in config:
    raise RuntimeError(
        "the parameter run_config is not defined. Define it with --config run_config=myfile.yml"
    )


run_config = parse_pileup_config(config["run_config"])


in_fld = run_config["input"]
out_fld = run_config["output"]
pileups = run_config["pileups"]
ref_records = run_config["ref_records"]

print("----------- run configuration ------------")
print("input folder:", in_fld)
print("\noutput folder:", out_fld)
print("\npileups:", json.dumps(pileups, indent=2))
print("\nreference records:", json.dumps(ref_records, indent=2))
print("------------------------------------------")


include: "rules/pileup.smk"


rule all:
    input:
        [
            expand(rules.map_reads.output, ref_id=ref, read_id=reads)
            for ref, reads in pileups.items()
        ],
        [
            expand(
                rules.build_pileup.output,
                ref_id=ref,
                rec_id=ref_records[ref],
                read_id=reads,
            )
            for ref, reads in pileups.items()
        ],
