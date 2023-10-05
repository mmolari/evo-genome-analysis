import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import pathlib


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--freq_thr",
        type=float,
        help="parameter: frequency threshold for trajectory selection",
    )
    parser.add_argument(
        "--cov_thr",
        type=int,
        help="parameter: threshold coverage for trajectory selection",
    )
    parser.add_argument(
        "--noise_thr",
        type=float,
        help="parameter: noise threshold for fwd/rev frequencies",
    )
    parser.add_argument(
        "--noise_tol",
        type=float,
        help="parameter: noise tolerance for fwd/rev frequencies",
    )
    parser.add_argument(
        "--max_initial_freq",
        type=float,
        help="parameter: max frequency at t0 to accept trajectory",
    )
    parser.add_argument(
        "--cov_window",
        type=int,
        help="parameter: convolution window for coverage smoothing n. insertion selection",
    )
    parser.add_argument(
        "--cov_fraction",
        type=float,
        help="parameter: fraction of the maximum coverage for n. insertion selection",
    )
    parser.add_argument(
        "--n_top_trajs",
        type=int,
        help="parameter: max. number of selected trajectories",
    )
    parser.add_argument("--ins_npz", type=str, help="input: npz file with insertions")
    parser.add_argument("--cov_npz", type=str, help="input: npz file with coverages")
    parser.add_argument(
        "--plot_fld", type=str, help="output: folder where to save plots"
    )
    return parser.parse_args()


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
        I = np.vstack([I.sum(axis=0), I])
        I_num.append(I)

        C = Cs[s]
        C = np.vstack([C.sum(axis=0), C])
        C = (C + np.roll(C, -1, axis=1)) / 2
        # C = np.minimum(C, np.roll(C, -1, axis=1))
        I_freq.append(safe_division(I, C, np.nan))
        C_num.append(C)
    # stack along new axis
    I_num = np.stack(I_num, axis=0)
    I_freq = np.stack(I_freq, axis=0)
    C_num = np.stack(C_num, axis=0)

    # shape: (S, 3, L)
    return I_num, I_freq, C_num, samples


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
    Fmax = np.nanmax(F[:, 0, :])
    bins = np.linspace(0, Fmax, 50)
    for s in range(S):
        Ft = F[s, 0, :]
        ax.hist(Ft, histtype="step", bins=bins, label=samples[s], color=c_dict[s])
    ax.axvline(freq_thr, color="k", linestyle="--")
    ax.set_xlabel("insertion frequency")
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
    Fmax = max(np.nanmax(dF), 1)
    bins = np.linspace(0, Fmax, 50)
    ax.hist(dF[mask], bins=bins)
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
    """Plot summary of insertion frequencies and deltaF"""
    fig, axs = plt.subplots(
        2, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1, 4]}, sharex="col"
    )
    __plot_freq(F, samples, freq_thr, axs[0, 0], axs[0, 1])
    __plot_dF(dF, axs[1, 0], axs[1, 1])
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
        ax.axhline(freq_thr, color="lightblue", linestyle="--")
        ax.axhline(f0_thr, color="gold", linestyle="--")

    for ax in axs[-1, :]:
        # ax.set_xlabel("sample")
        ax.set_xticks(np.arange(len(samples)))
        ax.set_xticklabels(samples, rotation=90)
        ax.set_ylim(bottom=0)
    for ax in axs[:, 0]:
        ax.set_ylabel("insertion frequency")
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


def plotly_insertions(F, samples, freq_thr, savename):
    """Plot interactive insertions distribution with plotly"""
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
        title=f"distribution of insertions - frequency threshold = {freq_thr}",
        xaxis_title="position (bp)",
        yaxis_title="prob. density",
        legend_title="samples",
    )
    fig.write_html(savename)


def n_insertion_scatterplot(
    In, Cn, samples, coverage_window, coverage_fraction, savename
):
    """Scatter-plot of number of insertions per position in the genome.
    This is compared with the rolling window average of the coverage."""
    positions = []

    S = len(samples)
    L = In.shape[2]
    fig, axs = plt.subplots(S, 1, figsize=(10, S * 4), sharex=True)
    pos = np.arange(L, dtype=int)
    kernel = np.ones(coverage_window) / coverage_window
    for k, s in enumerate(samples):
        ax = axs[k]

        for strand in [0, 1]:
            factor = 1 if strand == 0 else -1
            # rolling window average of Cn
            Cn_avg = np.convolve(Cn[k, strand + 1, :], kernel, mode="same")
            ax.scatter(
                pos,
                factor * Cn_avg,
                color="lightgray",
                marker=".",
                rasterized=True,
                zorder=-1,
            )
            I = In[k, strand + 1, :]
            mask = I > Cn_avg * coverage_fraction
            ax.scatter(pos, I * factor, marker=".", c=mask, cmap="bwr", rasterized=True)

            for i in pos[mask]:
                positions.append(
                    {
                        "sample": s,
                        "strand": "fwd" if strand == 0 else "rev",
                        "position": i,
                        "n_insertions": I[i],
                    }
                )

        ax.axhline(0, color="lightgray", linestyle="--")
        ax.set_title(s)
        ax.grid(alpha=0.2)
        ax.set_ylabel("fwd / rev n. insertions")
    axs[-1].set_xlabel("position (bp)")
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)

    return pd.DataFrame(positions)


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
    out_summary = out_fld / "summary.png"
    out_trajs = out_fld / "traj.png"
    out_html = out_fld / "insertions.html"
    out_csv_n = out_fld / "n_insertions.csv"
    out_csv_traj = out_fld / "insertion_traj.csv"
    out_n_insertions = out_fld / "n_insertions.png"

    print("loading data...")
    Is = load_npz(args.ins_npz)
    Cs = load_npz(args.cov_npz)

    print("creating tensors...")
    In, If, Cn, samples = insertion_arrays(Is, Cs)
    del Is, Cs

    print("computing delta-frequency...")
    dF = relevant_deltaF(
        If,
        Cn,
        freq_thr=freq_thr,
        cov_thr=cov_thr,
        noise_thr=noise_thr,
        noise_tol=noise_tol,
        max_fi=max_fi,
    )

    print("summary plot...")
    plot_summary(If, samples, freq_thr, dF, savename=out_summary)

    print("n. of insertions scatterplot...")
    df = n_insertion_scatterplot(
        In,
        Cn,
        samples,
        coverage_window=args.cov_window,
        coverage_fraction=args.cov_fraction,
        savename=out_n_insertions,
    )
    df.to_csv(out_csv_n, index=False)

    # create a dataframe with only relevant positions (deltaF >= 0)
    df = deltaF_to_dataframe(dF, If, samples)

    if len(df) > 0:
        print("plotting relevant trajectories...")
        # order clusters by max frequency, and select groups of trajectories
        idxs = df.index.values[:n_top_trajs]

        plot_traj(
            If, Cn, dF, samples, idxs, freq_thr, max_fi, cov_thr, savename=out_trajs
        )

        plotly_insertions(If, samples, freq_thr, savename=out_html)

        # export dataframe
        df.to_csv(out_csv_traj, index=False)
