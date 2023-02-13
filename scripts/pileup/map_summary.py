import pysam
import argparse

import numpy as np
import pandas as pd
from collections import defaultdict


def read_info_dict(read):
    res = {
        "qry": read.query_name,
        "ref": read.reference_name,
        "flag": read.flag,
        "ref_len": read.reference_length,
        "qL_cigar": read.infer_read_length(),
        "qL_seq": read.query_length,
        "sec": read.is_secondary,
        "suppl": read.is_supplementary,
    }
    return res


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
    Script that given a sam/bam file, prints a summary of how 
    many reads are mapped to each reference. Also surveys secondary
    and supplementary reads.
    """
    )
    parser.add_argument("--bam", type=str, help="bam file")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()
    print("### mapping summary for file: ###")
    print(args.bam, "\n")

    with pysam.Samfile(args.bam) as sf:

        ref_L = {r: sf.get_reference_length(r) for r in sf.references}

        print("----- references -----")
        for r, L in ref_L.items():
            print(f"{r:<15} = {L} bp")

        df = []
        for read in sf.fetch(until_eof=True):
            df.append(read_info_dict(read))
        df = pd.DataFrame(df)

        mask_unmapped = df.flag & 4 > 0
        df.loc[mask_unmapped, "ref"] = "unmapped"
        mask_primary = (~df.suppl) & (~df.sec) & (~mask_unmapped)

        print("\n----- primary mappings -----")
        sdf = df[mask_primary | mask_unmapped]
        summary = {
            "n. reads": sdf["ref"].value_counts(),
            "avg. read len (bp)": sdf.groupby("ref")["qL_seq"].mean(),
            "tot mapped len (bp)": sdf.groupby("ref")["qL_seq"].sum(),
        }
        summary["avg. coverage"] = summary["tot mapped len (bp)"] / pd.Series(ref_L)
        print(pd.DataFrame(summary))

        # primary reference for each query
        qid_to_primary_ref = {
            row.qry: row.ref for idx, row in df[mask_primary][["qry", "ref"]].iterrows()
        }

        # select supplementary mappings
        for tp in ["supplementary", "secondary"]:
            print(f"\n----- {tp} mappings -----")

            mask = df.suppl if tp == "supplementary" else df.sec
            N = mask.sum()
            print(f"n. {tp} mappings: {N}")

            if N == 0:
                continue

            sdf = df[mask].copy()
            sdf["primary ref"] = sdf["qry"].apply(lambda x: qid_to_primary_ref[x])
            srf = f"{tp} ref"
            sdf = sdf.rename(columns={"ref": srf})

            print("\n--> n. reads")
            print(sdf[["primary ref", srf]].value_counts().unstack())

            print(f"\n--> {tp} read total length")
            print(sdf.groupby(["primary ref", srf])["qL_seq"].sum().unstack())
