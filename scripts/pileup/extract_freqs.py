import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Evaluate consensus counts given a pileup and reference genome."""
    )
    parser.add_argument("--ref_genome", help="reference genome")
    parser.add_argument("--record", help="record to consider in the reference file")
    parser.add_argument("--pileup", help="reads pileup")
    parser.add_argument("--npz", help="output npz array")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def load_record_from_fasta(fasta_file, record_name):
    """Loads desired record from the fasta file"""
    ref_genome = None
    for rg in SeqIO.parse(fasta_file, format="fasta"):
        if rg.name == record_name:
            ref_genome = rg
    assert ref_genome is not None, f"record {record_name} not found in fasta file."
    return ref_genome


def ref_genome_to_nucl_id(ref_genome):
    """return an array of indices (0 to 5), one per nucleotide of a
    DNA sequence, according to the order (A,C,G,T,-,N)"""
    nuc_alpha = ["A", "C", "G", "T", "-", "N"]
    nuc_idx = {n: i for i, n in enumerate(nuc_alpha)}
    nseq = np.array(ref_genome)
    return np.array([nuc_idx[x] for x in nseq])


def consensus_count(pileup, ref_genome):
    """Returns the number of consensus nucleotides (n. times a nucleotide is
    the one on the reference genome) and the total number of nucleotides. Gaps are discarded"""
    ref_genome_idxs = ref_genome_to_nucl_id(ref_genome)
    nucl_pileup = pileup[:, :4, :]  # keep only nucleotides, no N or gap
    assert np.all(ref_genome_idxs < 4), "reference sequence includes gaps or Ns"
    # evaluate number of consensus
    n_consensus_f = np.choose(ref_genome_idxs, nucl_pileup[0])
    n_consensus_r = np.choose(ref_genome_idxs, nucl_pileup[1])
    n_tot = np.sum(nucl_pileup, axis=1)
    return {
        "ncons_fwd": n_consensus_f,
        "ncons_rev": n_consensus_r,
        "ntot_fwd": n_tot[0],
        "ntot_rev": n_tot[1],
    }

def gap_count(pileup):
    """Given the pileup retuns the count of fwd/rev gaps per position,
    along with the total number of reads per position including gaps."""
    ngaps = pileup[:, 4, :]
    ntot = pileup.sum(axis=1)
    return {
        "ngap_fwd" : ngaps[0],
        "ngap_rev" : ngaps[1],
        "ntot_fwd" : ntot[0],
        "ntot_rev" : ntot[1],
    }

if __name__ == "__main__":

    args = parse_args()

    if args.verbose:
        vprint = print
    else:
        vprint = lambda x: None

    vprint(f"extracting consensus frequency for {args.ref_genome} / {args.pileup}")

    vprint(f"load pileup and reference genome")
    pileup = np.load(args.pileup)["pileup"]
    ref_genome = load_record_from_fasta(args.ref_genome, args.record)
    assert len(ref_genome) == pileup.shape[2], "expected record length does not match."
    

    vprint(f"evaluate consensus count")
    cons_ct = consensus_count(pileup, ref_genome)

    vprint(f"evaluate gap count")
    gap_ct = gap_count(pileup)

    vprint(f"save results")
    cons_arr = pd.DataFrame(cons_ct)[["ncons_fwd", "ncons_rev", "ntot_fwd", "ntot_rev"]].to_numpy()
    gap_arr = pd.DataFrame(gap_ct)[["ngap_fwd", "ngap_rev", "ntot_fwd", "ntot_rev"]].to_numpy()
    np.savez_compressed(args.npz, cons_ct=cons_arr, gap_ct=gap_arr)

