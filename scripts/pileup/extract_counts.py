import argparse
import numpy as np
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pileups", type=str, nargs="+")
    parser.add_argument("--ref_seq", type=str)
    parser.add_argument("--ref_record", type=str)
    parser.add_argument("--cov_ct", type=str)
    parser.add_argument("--cons_ct", type=str)
    parser.add_argument("--gap_ct", type=str)
    return parser.parse_args()


def load_pileup(fname):
    """Load npz pileup file into a numpy array"""
    with np.load(fname) as data:
        # (2, 6, L)
        return data["pileup"]


def sample_name(fname):
    """Get the sample name from the pileup file"""
    return fname.split("/")[-2]


def pileup_to_coverage(pp):
    """Convert pileup to counts"""
    # (2, 6, L) -> (2, L)
    # exclude gap and ambiguous nucleotides
    return pp[:, :4, :].sum(axis=1)


def pileup_to_gap_count(pp):
    """Convert pileup to gap count"""
    # (2, 6, L) -> (2, L)
    # only consider gaps
    return pp[:, 4, :]


def pileup_to_consensus_count(pp, ref_seq):
    """Convert pileup to consensus count"""
    nucl_id = {n: i for i, n in enumerate("ACGT")}
    ref_vec = np.array([nucl_id[b] for b in ref_seq], dtype=int)

    L = len(ref_vec)
    cons_ct = pp[:, ref_vec, np.arange(L)]

    return cons_ct


def load_consensus_seq(fname, record):
    """Load the consensus sequence from the FASTA file"""
    # Load the FASTA file
    fasta = list(SeqIO.parse(fname, "fasta"))
    # Find the record with the same name as the pileup file
    for rec in fasta:
        if rec.id == record:
            # Return the consensus sequence
            return str(rec.seq)

    raise ValueError("Record not found in FASTA file")


if __name__ == "__main__":

    args = parse_args()

    ref_seq = load_consensus_seq(args.ref_seq, args.ref_record)

    # load pileups
    pps = {sample_name(s): load_pileup(s) for s in args.pileups}

    # evaluate fwd/rev coverage (2, L)
    cov_ct = {s: pileup_to_coverage(pp) for s, pp in pps.items()}

    # evaluate consensus count (2, L)
    cons_ct = {s: pileup_to_consensus_count(pp, ref_seq) for s, pp in pps.items()}

    # evaluate gap count (2, L)
    gap_ct = {s: pileup_to_gap_count(pp) for s, pp in pps.items()}

    # save coverage, consensus and gap counts
    np.savez_compressed(args.cov_ct, **cov_ct)
    np.savez_compressed(args.cons_ct, **cons_ct)
    np.savez_compressed(args.gap_ct, **gap_ct)
