import pysam
import argparse

import numpy as np
import pandas as pd


def read_info_dict(read):
    cigar = read.cigartuples
    Hstart = 0
    if cigar[0][0] == 5:
        Hstart = cigar[0][1]
    fwd = read.is_forward
    qs = read.query_alignment_start + Hstart
    qe = read.query_alignment_end + Hstart
    qL = read.infer_read_length()
    if not fwd:
        qs, qe = qL - qe, qL - qs
    res = {
        "qry": read.query_name,
        "ref": read.reference_name,
        "flag": read.flag,
        "fwd": fwd,
        "ref_len": read.reference_length,
        "qry_len": read.infer_read_length(),
        "n_matches": read.get_cigar_stats()[0][0],
        "sec": read.is_secondary,
        "suppl": read.is_supplementary,
        "rs": read.reference_start,
        "re": read.reference_end,
        "qs": qs,
        "qe": qe,
    }
    return res


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
    Script that given a sam/bam file, creates csv dataframe on info on secondary
    and supplementary mappings, with their corresponding primary one.
    """
    )
    parser.add_argument("--bam", type=str, help="sorted bam file")
    parser.add_argument("--csv", type=str, help="output csv dataframe")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    # list of unmapped reads
    non_primary_ids = []
    df = []

    with pysam.Samfile(args.bam) as sf:

        for read in sf.fetch():
            if read.is_secondary or read.is_supplementary:
                non_primary_ids.append(read.query_name)
        non_primary_ids = np.unique(non_primary_ids)

        for read in sf.fetch():
            if read.query_name in non_primary_ids:
                df.append(read_info_dict(read))

    if len(non_primary_ids) > 0:
        df = pd.DataFrame(df).sort_values(["qry", "flag", "ref"]).reset_index(drop=True)
        df.to_csv(args.csv, index=False)
