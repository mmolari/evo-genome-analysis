# %%
import numpy as np
from Bio import SeqIO


def load_pileup(fname):
    """Load npz pileup file into a numpy array"""
    with np.load(fname) as data:
        # (2, 6, L)
        pp = data["pileup"]
        # exclude gaps and ambiguous bases
        # (2, 4, L)
        return pp[:, :4, :]


def sample_name(fname):
    """Get the sample name from the pileup file"""
    return fname.split("/")[-2]


def pileup_to_counts(pp):
    """Convert pileup to counts"""
    # (2, 4, L) -> (2, L)
    return pp.sum(axis=1)


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


# %%
rec_id = "R4_51_kbp"
samples = [
    # f"../test_data/output-test-leo/pileup/reference_ST131I/R1_4963_kbp/sample{i}/allele_counts.npz"
    f"../test_data/output-test-leo/pileup/reference_ST131I/{rec_id}/sample{i}/allele_counts.npz"
    for i in range(1, 5)
]

ref = "../test_data/input-test-leo/references/reference_ST131I.fa"
ref_seq = load_consensus_seq(ref, rec_id)

pps = {sample_name(s): load_pileup(s) for s in samples}
cts = {s: pileup_to_counts(pp) for s, pp in pps.items()}
cons_ct = {s: pileup_to_consensus_count(pp, ref_seq) for s, pp in pps.items()}


# %%
