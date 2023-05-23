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
    parser.add_argument("--clip_dicts", type=str, nargs="+", help="clip dictionaries")
    parser.add_argument("--clip_ct", type=str, help="output clip counts")
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


def clip_dict_to_array(L, clips):
    """Given an insertion dictionary, returns an array with shape (6,L).
    The paris of dimension store:
    - fwd/rev count of clips
    - fwd/rev count of read ends
    - fwd/rev total size of clipped sequence at position"""

    C = np.zeros((6, L + 1), dtype=int)

    for p, v in clips["count"].items():
        C[:4, p] += v

    for p, d in clips["seqs"].items():
        fwd = d[0]
        rev = d[1]
        C[4, p] += sum(len(s) for s in fwd)
        C[5, p] += sum(len(s) for s in rev)
    return C


if __name__ == "__main__":
    args = parse_args()

    L = reference_seq_length(args.ref_fasta, args.ref_record)

    Cs = {}
    for cf in args.clip_dicts:
        with gzip.open(cf, "rb") as f:
            clip = pkl.load(f)

        sample = sample_name(cf)
        Cs[sample] = clip_dict_to_array(L, clip)

    np.savez_compressed(args.clip_ct, **Cs)
