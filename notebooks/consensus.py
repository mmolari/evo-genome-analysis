# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import plotly.figure_factory as ff

# %%


def load_npz_dict(fname):
    """Load a npz file and return a dictionary of arrays"""
    npz_dict = {}
    with np.load(fname) as f:
        for sample, data in f.items():
            npz_dict[sample] = data
    return npz_dict


def safe_division(a, b, extra):
    "return the division between the two vectors. Return `extra` where b is zero"
    mask = b == 0
    result = np.divide(a, b, out=np.zeros_like(a, dtype=float), where=~mask)
    result[mask] = extra
    return result


def load_tensors(cons_file, cov_file):
    """Load consensus and coverage files. Generate tensors of
    non-consensus frequencies and coverages"""

    # load files
    cons = load_npz_dict(cons_file)
    cov = load_npz_dict(cov_file)

    # always use this sample order
    samples = list(cons.keys())

    # evaluate frequencies
    freqs = []
    coverages = []
    for s in samples:
        # (2, L) tensors
        c = cons[s]
        n = cov[s]

        # 2 (L) tensor
        Ff, Fr = safe_division(c, n, extra=np.nan)
        Ft = safe_division(c.sum(axis=0), n.sum(axis=0), extra=np.nan)
        freqs.append(np.vstack([Ft, Ff, Fr]))

        coverages.append(np.vstack([n.sum(axis=0), n[0, :], n[1, :]]))
    # shape: (n_samples, 3, n_positions)
    F = np.stack(freqs, axis=0)
    C = np.stack(coverages, axis=0)

    # return non-consensus frequencies
    F = 1 - F
    return F, C, samples


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


def deltaF_to_dataframe(deltaF, F, samples):
    """Given the deltaF vector, return a dataframe with only positions that have deltaF > 0.
    The dataframe contains the following columns:
    - deltaF
    - pos
    - F_{sample} for each sample
    """

    # create dataframe and assign genomic positions
    df = pd.DataFrame(deltaF, columns=["deltaF"])
    df["pos"] = df.index + 1

    # only consider positive frequencies
    df = df[df["deltaF"] >= 0]

    # sort by deltaF
    df.sort_values("deltaF", ascending=False, inplace=True)

    # assign frequencies for each sample
    for ns, s in enumerate(samples):
        df[f"F_{s}"] = F[ns, 0, df.index]

    return df


def __plot_freq(F, samples, freq_thr, ax1, ax2):
    """Plot non-consensus frequency histogram and high frequency positions"""

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
    ax.set_xlabel("non-consensus frequency")
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


def __plot_dF(dF, ax1, ax2):
    """Plot deltaF histogram and positions"""

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
    """Plot summary of non-consensus frequencies and deltaF"""
    fig, axs = plt.subplots(
        2, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1, 4]}, sharex="col"
    )
    __plot_freq(F, samples, freq_thr, axs[0, 0], axs[0, 1])
    __plot_dF(dF, axs[1, 0], axs[1, 1])
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


def plot_traj(F, C, dF, samples, idxs, freq_thr, cov_thr, savename):
    """Plot selected frequency trajectories"""

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
        ax.set_ylabel("non-consensus frequency")
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


def plotly_noncons(F, samples, freq_thr, savename):
    """Plot interactive non-consensus distribution with plotly"""
    Ft = F[:, 0, :]

    keep_samples = []
    sample_nums = []
    for ns, s in enumerate(samples):
        mask = Ft[ns, :] > freq_thr
        if np.any(mask):
            keep_samples.append(s)
            sample_nums.append(ns)

    cmap = mpl.colormaps.get("tab10")
    colors = [mpl.colors.to_hex(cmap(s)) for s in sample_nums]
    hist_data = []
    for s in sample_nums:
        mask = Ft[s, :] > freq_thr
        pos = np.argwhere(mask).flatten() + 1
        hist_data.append(pos)
    fig = ff.create_distplot(
        hist_data, keep_samples, bin_size=1000, show_curve=False, colors=colors
    )
    fig.update_layout(
        title=f"distribution of non-consensus sites - frequency threshold = {freq_thr}",
        xaxis_title="position (bp)",
        yaxis_title="prob. density",
        legend_title="samples",
    )
    fig.write_html(savename)


# %%
cov_thr = 5
freq_thr = 0.75
n_top = 33

# rec_id = "R1_4963_kbp"
# rec_id = "R2_127_kbp"
# rec_id = "R3_89_kbp"
rec_id = "R4_51_kbp"
cons_file = (
    f"../test_data/output-test-leo/pileup/reference_ST131I/{rec_id}/consensus.npz"
)
cov_file = f"../test_data/output-test-leo/pileup/reference_ST131I/{rec_id}/coverage.npz"

# %%


# (n_samples, 3, L) tensors of consensus frequencies and coverages, plust ordered list of samples
# the second index of the tensor is: 0: tot, 1: fwd, 2: rev
F, C, samples = load_tensors(cons_file, cov_file)


# compute deltaF for each position
deltaF = relevant_deltaF(F, C, freq_thr=freq_thr, cov_thr=cov_thr)

out_summary = f"test_fig/consensus_summary.png"
plot_summary(F, samples, freq_thr, deltaF, savename=out_summary)

# create a dataframe with only relevant positions (deltaF >= 0)
df = deltaF_to_dataframe(deltaF, F, samples)


# %%

if len(df) > 0:
    # order clusters by max frequency, and select groups of trajectories
    idxs = df.index.values[:n_top]

    out_trajs = f"test_fig/consensus_trajs.png"
    plot_traj(F, C, deltaF, samples, idxs, freq_thr, cov_thr, savename=out_trajs)

    out_html = f"test_fig/consensus.html"
    plotly_noncons(F, samples, freq_thr, savename=out_html)

    # export dataframe
    out_csv = f"test_fig/consensus.csv"
    df.to_csv(out_csv, index=False)

# %%
