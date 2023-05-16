# %%

import pickle as pkl
import gzip
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


def load_insertion_dict(fname):
    with gzip.open(fname, "rb") as f:
        insertions = pkl.load(f)
    return insertions


def load_coverage(fname):
    all_count = {}
    with np.load(fname) as f:
        all_count = f["pileup"]
    return np.sum(all_count, axis=1)


def safe_division(a, b, extra):
    "return the division between the two vectors. Return extra where b is zero"
    mask = b == 0
    result = np.divide(a, b, out=np.zeros_like(a, dtype=float), where=~mask)
    result[mask] = extra
    return result


# %%

# rec_id = "R1_4963_kbp"
# rec_id = "R2_127_kbp"
# rec_id = "R3_89_kbp"
rec_id = "R4_51_kbp"
sample_id = "sample3"

prefix = f"../test_data/output-test-leo/pileup/reference_ST131I/{rec_id}/{sample_id}"

ins_file = f"{prefix}/insertions.pkl.gz"
all_ct_file = f"{prefix}/allele_counts.npz"
# %%
insertions = load_insertion_dict(ins_file)
ac = load_coverage(all_ct_file)
L = ac.shape[1]

# %%
def ins_numbers(insertions, L):
    """Return (2, L) matrices of insertion number and total length"""
    In = np.zeros((2, L), dtype=int)
    Il = np.zeros((2, L), dtype=int)
    for k, ins_dict in insertions.items():
        for seq, n_ins in ins_dict.items():
            In[:, k] += n_ins
            Il[:, k] += n_ins * len(seq)
    return In, Il


In, Il = ins_numbers(insertions, L)
Cn = (ac + np.roll(ac, -1)) / 2
# Cnmin = np.minimum(ac, np.roll(ac, 1))
# %%
If = In.sum(axis=0) / Cn.sum(axis=0)
plt.hist(If, bins=100)
plt.yscale("log")
plt.show()

# %%
plt.scatter(In.mean(axis=0), Il.mean(axis=0), marker=".", alpha=0.1)
plt.show()

plt.scatter(Cn.sum(axis=0), In.sum(axis=0), c=Il.sum(axis=0), marker=".", alpha=0.1)
plt.show()
# %%

# problem with quality filter -> coverage of flanking sites as a proxy
# test mapping and make sure it works (but I think it does)
# select a step and parallelize with snakemake?
# what quantities?
# - define better dimentsions (L-2) for insertions and average coverage,
# - make sure that everything matches
# - then per window:
#   - n. insertions vs coverage
#   - avg. inserted length vs coverage.

Int = In.sum(axis=0)
Ilt = Il.sum(axis=0) / Int
Cnt = Cn.sum(axis=0)
# %%
step = 1
Inc = np.convolve(np.ones(step), Int, mode="same")
Cnc = np.convolve(np.ones(step) / step, Cnt, mode="same")
# %%
plt.hist(
    Inc / Cnc,
    bins=100,
)
plt.yscale("log")
plt.show()

freq_threshold = 0.6
x = np.arange(L)
y = Inc / Cnc
mask = y > freq_threshold

plt.scatter(x, y, color="gray", marker=".")
plt.axhline(freq_threshold, color="k", ls="--")
plt.scatter(x[mask], y[mask], color="red", marker=".")
plt.tight_layout()
plt.show()
# %%
