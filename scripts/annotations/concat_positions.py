import numpy as np
import pandas as pd
import pathlib as pl
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gaps", type=str, required=True)
    parser.add_argument("--cons", type=str, required=True)
    parser.add_argument("--ins", type=str, required=True)
    parser.add_argument("--clips", type=str, required=True)
    parser.add_argument("--out", type=str, required=True)
    return parser.parse_args()


def try_read_csv(fld, fname):
    try:
        return pd.read_csv(fld / fname)
    except FileNotFoundError:
        return None


if __name__ == "__main__":
    args = parse_args()

    dfs = {}
    for label, plab, fld, fname in [
        ("consensus", "pos", pl.Path(args.cons), "consensus.csv"),
        ("insertion", "position", pl.Path(args.ins), "n_insertions.csv"),
        ("clip", "position", pl.Path(args.clips), "n_clips.csv"),
        ("gap", "pos", pl.Path(args.gaps), "gap.csv"),
    ]:
        # make sure path exists
        assert fld.exists(), f"Path {fld} does not exist"
        df = try_read_csv(fld, fname)
        if df is not None:
            assert plab in df.columns, f"Column {plab} not found in {fld / fname}"
            sdf = df[plab].unique()
            sdf = pd.DataFrame(sdf, columns=["position"])
            sdf["label"] = label
            dfs[label] = sdf

    # concatenate
    df = pd.concat(dfs.values())
    df = df.sort_values("position")

    # save
    df.to_csv(args.out, index=False)
