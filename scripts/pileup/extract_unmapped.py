import pysam
import argparse
import gzip

import numpy as np
import pandas as pd

from Bio import SeqIO, SeqRecord, Seq


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
    Script that given a sam/bam file, creates a .fastq.gz file with unmapped reads
    and a .csv dataframe with stats on each read (name, length, avg_quality, flag).
    """
    )
    parser.add_argument("--bam", type=str, help="sorted bam file")
    parser.add_argument("--csv", type=str, help="output csv dataframe")
    parser.add_argument("--fastq", type=str, help="output fastq.gz file")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()

    # list of unmapped reads
    unmapped_df = []
    unmapped_records = []

    with pysam.Samfile(args.bam) as sf:

        for read in sf.fetch():

            # capture only unmapped reads
            if not read.is_unmapped:
                continue

            unmapped_df.append(
                {
                    "read": read.query_name,
                    "len": len(read.seq),
                    "avg. qscore": np.mean(read.query_qualities),
                    "flag": read.flag,
                }
            )

            # create and append fastq-compatible record
            seq = Seq.Seq(read.seq)
            rec = SeqRecord.SeqRecord(seq, id=read.query_name, description="", name="")
            rec.letter_annotations["phred_quality"] = read.query_qualities
            unmapped_records.append(rec)

    # create and export dataframe
    df = pd.DataFrame(unmapped_df)
    if len(df) > 0:
        df = df.sort_values("len", ascending=False)
    df.to_csv(args.csv, index=False)

    # create and export fastq.gz file
    with gzip.open(args.fastq, "wt") as f:
        SeqIO.write(unmapped_records, f, "fastq")
