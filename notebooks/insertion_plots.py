# %%

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


def load_npz(fname):
    A = {}
    with np.load(fname) as f:
        for s in f:
            A[s] = f[s]
    return A


def safe_division(a, b, extra):
    "return the division between the two vectors. Return extra where b is zero"
    mask = b == 0
    result = np.divide(a, b, out=np.zeros_like(a, dtype=float), where=~mask)
    result[mask] = extra
    return result


def insertion_arrays(Is, Cs):
    """Given the insertion and coverage dictionaries,
    returns the tensor with insertion number, frequencies and coverage,
    together with the order of samples"""
    samples = sorted(list(Is.keys()))
    I_num, I_freq, C_num = [], [], []
    for s in samples:
        I = Is[s][:2]
        I = np.vstack([I, I.sum(axis=0)])
        I_num.append(I)

        C = Cs[s]
        C = np.vstack([C, C.sum(axis=0)])
        C = (C + np.roll(C, -1)) / 2
        # C = np.maximum(C, np.roll(C, -1))
        I_freq.append(safe_division(I, C, np.nan))
        C_num.append(C)
    # stack along new axis
    I_num = np.stack(I_num, axis=0)
    I_freq = np.stack(I_freq, axis=0)
    C_num = np.stack(C_num, axis=0)

    # shape: (S, 3, L)
    return I_num, I_freq, C_num, samples


def relevant_deltaF(F, C, freq_thr, cov_thr):
    """
    Evaluate a per-position delta-frequency (Fmax - Fmin) vector.
    Only sites that have above thresholds fwd and rev coverage are considered
    for Fmax and Fmin.
    Moreover, only sites with above frequency-threshold Fmax are considered.
    All the other sites are assigned a deltaF of -inf.
    """

    # mask relevant coverages (sites that have fwd & rev coverage above threshold)
    # shape (n_samples, L)
    mask_fwd = C[:, 1, :] >= cov_thr
    mask_rev = C[:, 2, :] >= cov_thr
    mask_cov = mask_fwd & mask_rev

    # evaluate maximum/minimum non-masked frequency
    allnan = np.all(~mask_cov, axis=0)

    F_max = np.copy(F[:, 0, :])
    F_max[~mask_cov] = np.nan
    F_max = np.nanmax(F_max, axis=0, initial=0)
    F_max[allnan] = np.nan

    F_min = np.copy(F[:, 0, :])
    F_min[~mask_cov] = np.nan
    F_min = np.nanmin(F_min, axis=0, initial=1)
    F_min[allnan] = np.nan

    # evaluate distance between covered max and min frequencies
    # shape (L)
    dF = F_max - F_min

    # at least one above-frequency-threshold and above-coverage point
    # shape (L)
    mask = F_max > freq_thr

    # set masked positions to -inf
    dF[~mask] = -np.inf

    return dF


# %%

prefix = "../test_data/output-test-leo/pileup/reference_ST131I"
ref_id = "R2_127_kbp"

ins_file = f"{prefix}/{ref_id}/insertions.npz"
cov_file = f"{prefix}/{ref_id}/coverage.npz"

freq_thr = 0.8
cov_thr = 10

# %%


Is = load_npz(ins_file)
Cs = load_npz(cov_file)
In, If, Cn, samples = insertion_arrays(Is, Cs)
dF = relevant_deltaF(If, Cn, freq_thr, cov_thr)
S = len(samples)
# %%

Ifmax = np.nanmax(If[:, 0, :])
bins = np.linspace(0, Ifmax, 25)

fig, axs = plt.subplots(2, 2, figsize=(10, 6), gridspec_kw={"width_ratios": [1, 3]})

ax = axs[0, 0]
for ns, sample in enumerate(samples):
    ax.hist(If[ns, 0, :], bins=bins, histtype="step")
ax.set_xlabel("n. insertions")
ax.set_ylabel("n. positions")
ax.set_yscale("log")

plt.tight_layout()
plt.show()
# %%

# 4-panel plot for distribution of insertions
S = len(Is)
fig, axs = plt.subplots(S, 1, figsize=(10, S * 4), sharex=True)
axs = axs.flatten()
for i, (s, I) in enumerate(Is.items()):
    axs[i].plot(I[0, :], marker=".", ls="")
    axs[i].plot(-I[1, :], marker=".", ls="")
    axs[i].set_title(s)
plt.show()

# %%

for s, I in Is.items():
    plt.plot(I[0, :], marker=".", ls="")
    plt.plot(-I[1, :], marker=".", ls="")
plt.show()


# %%

Max = max([I[:2].max() for I in Is.values()])
for s, I in Is.items():
    bins = np.arange(0, Max + 2, 1) - 0.5
    plt.hist(I[:2].sum(axis=0), bins=bins, histtype="step", label=s)
plt.yscale("log")
plt.show()
# %%
