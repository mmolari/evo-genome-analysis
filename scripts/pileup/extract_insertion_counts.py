import argparse
import numpy as np
import pickle as pkl
import gzip
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ref_fasta", type=str, help="FASTA file with reference sequence"
    )
    parser.add_argument("--ref_record", type=str, help="Record name in FASTA file")
    parser.add_argument(
        "--ins_dicts", type=str, nargs="+", help="Insertion dictionaries"
    )
    parser.add_argument("--ins_ct", type=str, help="output insertion counts")
    return parser.parse_args()


def reference_seq_length(fname, record):
    """return the length of the reference sequence from the FASTA file"""
    # Load the FASTA file
    fasta = list(SeqIO.parse(fname, "fasta"))
    # Find the record with the same name as the pileup file
    for rec in fasta:
        if rec.id == record:
            # Return the consensus sequence
            return len(str(rec.seq))

    raise ValueError("Record not found in FASTA file")


def sample_name(fname):
    """Get the sample name from the insertion file"""
    return fname.split("/")[-2]


def ins_dict_to_ins_counts(L, ins):
    """Given an insertion dictionary, returns an array with shape (4,L).
    The first dimension stores:
    - fwd/rev count of insertions
    - fwd/rev total size of inserted sequence at position"""

    I = np.zeros((4, L), dtype=int)

    for p, d in ins.items():
        for s, v in d.items():
            l = len(s)
            I[:2, p] += v
            I[2:, p] += v * l
    return I


if __name__ == "__main__":
    args = parse_args()

    L = reference_seq_length(args.ref_fasta, args.ref_record)

    ins_files = args.ins_dicts

    I = {}
    for fname in ins_files:
        with gzip.open(fname, "rb") as f:
            ins = pkl.load(f)
        sample = sample_name(fname)
        I[sample] = ins_dict_to_ins_counts(L, ins)

    np.savez_compressed(args.ins_ct, **I)
