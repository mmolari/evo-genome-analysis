# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd


def load_npz_dict(fname):
    npz_dict = {}
    with np.load(fname) as f:
        for sample, data in f.items():
            npz_dict[sample] = data
    return npz_dict


def safe_division(a, b, extra):
    "return the division between the two vectors. Return extra where b is zero"
    mask = b == 0
    result = np.divide(a, b, out=np.zeros_like(a, dtype=float), where=~mask)
    result[mask] = extra
    return result


def load_tensors(gap_file, cov_file):
    """Load gap and coverage files and generate tensors of gap frequencies and coverages"""
    gaps = load_npz_dict(gap_file)
    cov = load_npz_dict(cov_file)

    # always use this sample order
    samples = list(gaps.keys())

    # evaluate frequencies
    freqs = []
    coverages = []
    for s in samples:
        g = gaps[s]
        c = cov[s]
        Ff = safe_division(g[0, :], c[0, :] + g[0, :], extra=0)
        Fr = safe_division(g[1, :], c[1, :] + g[1, :], extra=0)
        gtot = g.sum(axis=0)
        ctot = c.sum(axis=0)
        Ft = safe_division(gtot, ctot + gtot, extra=0)
        freqs.append(np.vstack([Ft, Ff, Fr]))
        gapped_cov = [ctot + gtot, c[0, :] + g[0, :], c[1, :] + g[1, :]]
        coverages.append(np.vstack(gapped_cov))
    # shape: (n_samples, 3, n_positions)
    F = np.stack(freqs, axis=0)
    C = np.stack(coverages, axis=0)
    return F, C, samples


def relevant_deltaF(F, C, freq_thr, cov_thr):
    """
    Evaluate a per-position delta-frequency (Fmax - Fmin).
    Only sites that have above thresholds fwd and rev coverage are considered
    for Fmax and Fmin.
    Moreover, only sites with above threshold Fmax are returned.
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

    # exclude masked positions
    dF[~mask] = -np.inf
    return dF


def plot_gap_overview(F, samples, freq_thr, ax1, ax2):

    S = len(samples)
    cmap = mpl.colormaps.get("tab10")
    c_dict = {s: cmap(s) for s in range(S)}

    # plot histogram
    ax = ax1
    bins = np.linspace(0, 1, 50)
    for s in range(S):
        Ft = F[s, 0, :]
        ax.hist(Ft, histtype="step", bins=bins, label=samples[s], color=c_dict[s])
    ax.axvline(freq_thr, color="k", linestyle="--")
    ax.set_xlabel("gap frequency")
    ax.set_ylabel("n. sites")
    ax.set_yscale("log")
    ax.legend(loc="upper right")

    # plot high frequency positions
    ax = ax2
    for s in range(S):
        Ft = F[s, 0, :]
        mask = Ft > freq_thr
        x = np.argwhere(mask).flatten() + 1
        y = np.ones_like(x) * s
        ax.plot(x, y, label=s, color=c_dict[s], marker="|", linestyle="None")
    ax.set_xlabel("position (bp)")
    ax.set_yticks(np.arange(S))
    ax.set_yticklabels(samples)
    ax.grid(alpha=0.2)


def plot_dF(dF, ax1, ax2):

    mask = dF >= 0

    # plot histogram
    ax = ax1
    ax.hist(dF[mask], bins=50)
    ax.set_xlabel("dF")
    ax.set_ylabel("n. sites")
    ax.set_yscale("log")
    ax.set_xlim(left=0)

    # plot high frequency positions
    ax = ax2
    x = np.arange(len(dF))
    x = x[mask]
    y = dF[mask]
    ax.scatter(x, y, marker="|", alpha=0.3)
    ax.set_xlabel("position (bp)")
    ax.set_ylabel("dF")
    axt = ax.twinx()
    axt.hist(x, bins=100, histtype="step", color="k", alpha=0.2)
    axt.set_ylabel("n. sites")
    ax.grid(alpha=0.2)


def plot_summary(F, samples, freq_thr, dF, savename):
    fig, axs = plt.subplots(
        2, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1, 4]}, sharex="col"
    )
    plot_gap_overview(F, samples, freq_thr, axs[0, 0], axs[0, 1])
    plot_dF(dF, axs[1, 0], axs[1, 1])
    plt.tight_layout()
    plt.savefig(savename)
    plt.show()


def plot_traj(F, C, dF, samples, idxs, freq_thr, cov_thr, savename):
    "Plot selected frequency trajectories"

    I = len(idxs)
    # number of columns and rows
    nC = 3
    nR = I // nC + 1 * (I % nC > 0)
    if I == 0:
        nR = 1
    fig, axs = plt.subplots(
        nR, nC, figsize=(4 * nC, 2 * nR), sharex=True, sharey=True, squeeze=False
    )

    for ni, i in enumerate(idxs):
        ax = axs[ni // nC, ni % nC]

        A = 0.3 + 0.7 * (C[:, :, i] >= cov_thr)
        x = np.arange(len(samples))
        ax.plot(
            F[:, 0, i],
            label="tot",
            color="k",
            marker="o",
            linestyle=":",
            zorder=1,
        )
        ax.scatter(
            x,
            F[:, 1, i],
            label="fwd",
            alpha=A[:, 1],
            color="C0",
            marker=">",
        )
        ax.scatter(
            x,
            F[:, 2, i],
            label="rev",
            alpha=A[:, 2],
            color="C1",
            marker="<",
        )
        ax.set_title(f"pos {i+1}  |  dF={dF[i]:.2f}")
        ax.grid(alpha=0.2)
        ax.axhline(freq_thr, color="lightgray", linestyle="--")

    for ax in axs[-1, :]:
        # ax.set_xlabel("sample")
        ax.set_xticks(np.arange(len(samples)))
        ax.set_xticklabels(samples, rotation=90)
        ax.set_ylim(0, 1)
    for ax in axs[:, 0]:
        ax.set_ylabel("gap frequency")
    plt.tight_layout()
    plt.savefig(savename)
    plt.show()


# %%

# define parameters

cov_thr = 5
freq_thr = 0.75

# rec_id = "R1_4963_kbp"
# rec_id = "R2_127_kbp"
rec_id = "R3_89_kbp"
# rec_id = "R4_51_kbp"
gap_file = f"../test_data/output-test-leo/pileup/reference_ST131I/{rec_id}/gaps.npz"
cov_file = f"../test_data/output-test-leo/pileup/reference_ST131I/{rec_id}/coverage.npz"

# (n_samples, 3, L) tensors of gap frequencies and coverages, plust ordered list of samples
# the second index of the tensor is: 0: tot, 1: fwd, 2: rev
F, C, samples = load_tensors(gap_file, cov_file)


# %%


deltaF = relevant_deltaF(F, C, freq_thr=freq_thr, cov_thr=cov_thr)

# %%

# to dataframe
df = pd.Series(deltaF, name="deltaF").to_frame()
df["pos"] = np.arange(len(deltaF)) + 1
df = df[df["deltaF"] >= 0]
df.sort_values("deltaF", ascending=False, inplace=True)
idxs = df["pos"].values - 1
for ns, s in enumerate(samples):
    df[f"F_{s}"] = F[ns, 0, idxs]

# %%
n_top = 33
idxs = np.array(df.iloc[:n_top]["pos"].values - 1)

# %%

svfig = f"test_fig/gap_summary_{rec_id}.png"
plot_summary(F, samples, freq_thr, deltaF, savename=svfig)

svfig = f"test_fig/gap_traj_{rec_id}.png"
plot_traj(F, C, deltaF, samples, idxs, freq_thr, cov_thr, savename=svfig)
# %%
# export
df.to_csv(f"test_fig/gap_summary_{rec_id}.csv", index=False)
# %%
df
# %%
