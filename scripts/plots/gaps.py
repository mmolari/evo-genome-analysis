import argparse
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import plotly.figure_factory as ff


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freq_thr", type=float)
    parser.add_argument("--cov_thr", type=int)
    parser.add_argument("--noise_thr", type=float)
    parser.add_argument("--noise_tol", type=float)
    parser.add_argument("--max_initial_freq", type=float)
    parser.add_argument("--n_top_trajs", type=int)
    parser.add_argument("--gap_npz", type=str)
    parser.add_argument("--cov_npz", type=str)
    parser.add_argument("--plot_fld", type=str)
    return parser.parse_args()


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


def relevant_deltaF(F, C, freq_thr, cov_thr, noise_thr, noise_tol, max_fi):
    """
    Evaluate a per-position delta-frequency (Fmax - Fmin) vector.
    Only sites that have above thresholds fwd and rev coverage are considered
    for Fmax and Fmin.
    Only sites with above frequency-threshold Fmax are considered.
    All the other sites are assigned a deltaF of -inf.
    Also, sites in which |freq_fwd - freq_rev| > noise_thr * avg_freq + noise_tol are discarded.
    Finally, trajectories whose initial frequency is above max_fi are discarded.
    """

    # true where site should not be discarded
    mask_keep = np.ones_like(F[:, 0, :], dtype=bool)

    # mask relevant coverages (sites that have fwd & rev coverage above threshold)
    # shape (n_samples, L)
    mask_fwd = C[:, 1, :] >= cov_thr
    mask_rev = C[:, 2, :] >= cov_thr
    mask_keep &= mask_fwd & mask_rev

    # remove sites that have |fwd_freq-rev_freq| < avg_freq * noise_thr + noise_tol
    noise_mask = np.abs(F[:, 1, :] - F[:, 2, :]) <= F[:, 0, :] * noise_thr + noise_tol
    mask_keep &= noise_mask

    # evaluate maximum/minimum non-masked frequency
    allnan = np.all(~mask_keep, axis=0)

    F_max = np.copy(F[:, 0, :])
    F_max[~mask_keep] = np.nan
    F_max = np.nanmax(F_max, axis=0, initial=0)
    F_max[allnan] = np.nan

    F_min = np.copy(F[:, 0, :])
    F_min[~mask_keep] = np.nan
    F_min = np.nanmin(F_min, axis=0, initial=1)
    F_min[allnan] = np.nan

    # evaluate distance between covered max and min frequencies
    # shape (L)
    dF = F_max - F_min

    # exclude sites that have initial frequency above max_fi at the first timepoint
    mask = F[0, 0, :] <= max_fi
    dF[~mask] = -np.inf

    # at least one above-frequency-threshold and above-coverage point
    # shape (L)
    mask = F_max > freq_thr

    # set masked positions to -inf
    dF[~mask] = -np.inf

    return dF


def cluster_pos_dict(pos):
    """assigns cluster ids to positions. Sets of adjacent positions
    are assigned to the same cluster. Returns a dictionary mapping."""
    p = np.sort(pos)
    dp = np.diff(p)
    clusters = np.zeros_like(p)
    clusters[0] = 1
    clusters[1:] = dp > 1
    clusters = np.cumsum(clusters)

    pos_to_cluster = dict(zip(p, clusters))
    return pos_to_cluster


def deltaF_to_dataframe(deltaF, F, samples):
    """Given the deltaF vector, return a dataframe with only positions that have deltaF > 0.
    The dataframe contains the following columns:
    - deltaF
    - pos
    - F_{sample} for each sample
    - gap_cluster: cluster id for each position
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

    if len(df) == 0:
        return df

    cluster_dict = cluster_pos_dict(df["pos"])

    df["gap_cluster"] = df["pos"].map(cluster_dict)

    return df


def plot_gap_overview(F, samples, freq_thr, ax1, ax2):
    """Plot gap frequency histogram and high frequency positions"""

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
    """Plot summary of gap frequencies and deltaF"""
    fig, axs = plt.subplots(
        2, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1, 4]}, sharex="col"
    )
    plot_gap_overview(F, samples, freq_thr, axs[0, 0], axs[0, 1])
    plot_dF(dF, axs[1, 0], axs[1, 1])
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


def plot_traj(F, C, dF, samples, idxs, freq_thr, f0_thr, cov_thr, savename):
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

    for ni, idx in enumerate(idxs):
        ax = axs[ni // nC, ni % nC]

        for i in idx:
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
        if len(idx) == 1:
            i = idx[0]
            ax.set_title(f"pos {i+1}  |  dF={dF[i]:.2f}")
        else:
            avg_dF = np.mean(dF[idx])
            ax.set_title(
                f"pos [{np.min(idx)+1}-{np.max(idx)+1}]  |  <dF> = {avg_dF:.2f}"
            )

        ax.grid(alpha=0.2)
        ax.axhline(freq_thr, color="lightblue", linestyle="--")
        ax.axhline(f0_thr, color="gold", linestyle="--")

    for ax in axs[-1, :]:
        # ax.set_xlabel("sample")
        ax.set_xticks(np.arange(len(samples)))
        ax.set_xticklabels(samples, rotation=90)
        ax.set_ylim(0, 1)
    for ax in axs[:, 0]:
        ax.set_ylabel("gap frequency")
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


def plotly_gaps(F, samples, freq_thr, savename):
    """Plot interactive gap distribution with plotly"""
    Ft = F[:, 0, :]
    cmap = mpl.colormaps.get("tab10")
    hist_data, keep_samples, colors = [], [], []
    for s, sample in enumerate(samples):
        mask = Ft[s, :] > freq_thr
        if mask.sum() > 0:
            pos = np.argwhere(mask).flatten() + 1
            hist_data.append(pos)
            colors.append(mpl.colors.to_hex(cmap(s)))
            keep_samples.append(sample)

    fig = ff.create_distplot(
        hist_data, keep_samples, bin_size=1000, show_curve=False, colors=colors
    )
    fig.update_layout(
        title=f"Gap distribution - frequency threshold = {freq_thr}",
        xaxis_title="position (bp)",
        yaxis_title="prob. density",
        legend_title="samples",
    )
    fig.write_html(savename)


if __name__ == "__main__":
    args = parse_args()

    # define parameters
    cov_thr = args.cov_thr
    freq_thr = args.freq_thr
    noise_thr = args.noise_thr
    noise_tol = args.noise_tol
    max_fi = args.max_initial_freq
    n_top_trajs = args.n_top_trajs

    out_fld = pathlib.Path(args.plot_fld)
    out_fld.mkdir()
    out_summary = out_fld / "gap_summary.png"
    out_trajs = out_fld / "gap_trajs.png"
    out_html = out_fld / "gap_distribution.html"
    out_csv = out_fld / "gap.csv"

    # (n_samples, 3, L) tensors of gap frequencies and coverages, plust ordered list of samples
    # the second index of the tensor is: 0: tot, 1: fwd, 2: rev
    F, C, samples = load_tensors(args.gap_npz, args.cov_npz)

    # compute deltaF for each position
    deltaF = relevant_deltaF(
        F,
        C,
        freq_thr=freq_thr,
        cov_thr=cov_thr,
        noise_thr=noise_thr,
        noise_tol=noise_tol,
        max_fi=max_fi,
    )

    plot_summary(F, samples, freq_thr, deltaF, savename=out_summary)

    # create a dataframe with only relevant positions (deltaF >= 0)
    df = deltaF_to_dataframe(deltaF, F, samples)

    # do not plot if no relevant positions are found
    if len(df) > 0:
        # order clusters by max frequency, and select groups of trajectories
        cl_freq = df.groupby("gap_cluster")["deltaF"].max().sort_values(ascending=False)
        cl_idxs = cl_freq.index.values[:n_top_trajs]
        idxs_groups = [df["pos"][df["gap_cluster"] == i].values - 1 for i in cl_idxs]

        plot_traj(
            F,
            C,
            deltaF,
            samples,
            idxs_groups,
            freq_thr,
            max_fi,
            cov_thr,
            savename=out_trajs,
        )
        plotly_gaps(F, samples, freq_thr, savename=out_html)

        # export dataframe
        df.to_csv(out_csv, index=False)
