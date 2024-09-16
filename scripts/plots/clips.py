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
        help="parameter: convolution window for coverage smoothing n. clips selection",
    )
    parser.add_argument(
        "--cov_fraction",
        type=float,
        help="parameter: fraction of the maximum coverage for n. clips selection",
    )
    parser.add_argument(
        "--n_top_trajs",
        type=int,
        help="parameter: max. number of selected trajectories",
    )
    parser.add_argument("--clips_npz", type=str, help="input: npz file with clips")
    parser.add_argument("--cov_npz", type=str, help="input: npz file with coverages")
    parser.add_argument(
        "--plot_fld", type=str, help="output: folder where to save plots"
    )
    parser.add_argument("--samples_order", nargs="+", type=str, help="order of samples")
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


def load_data(clip_file, coverage_file):
    Cl_dict = load_npz(clip_file)
    Cov_dict = load_npz(coverage_file)

    samples = sorted(list(Cl_dict.keys()))
    Cls, Covs = [], []
    for s in samples:
        Cl = Cl_dict[s][:2]
        Cl = np.vstack([Cl.sum(axis=0), Cl])
        Cls.append(Cl)

        Cov = Cov_dict[s]
        Cov = np.vstack([Cov.sum(axis=0), Cov])
        Covs.append(Cov)

    Cls = np.stack(Cls, axis=0)  # shape: (S, 3, L)
    Covs = np.stack(Covs, axis=0)  # shape: (S, 3, L)

    return Cls, Covs, samples


def coverage_comparison(Covs):
    col_1 = Covs[:, :, 0]
    middle = np.maximum(Covs[:, :, 1:], Covs[:, :, :-1])
    col_n = Covs[:, :, -1]

    # concatenate
    comparison_coverages = np.concatenate(
        [col_1[:, :, None], middle, col_n[:, :, None]], axis=-1
    )

    return comparison_coverages


def n_clip_scatterplots(
    Cls, Covs, samples, convolution_window, freq_threshold, savename
):
    S = len(samples)
    L = Cls.shape[-1]
    pos = np.arange(L, dtype=int)
    fig, axs = plt.subplots(S, 1, figsize=(12, 3 * S), sharex=True)

    kernel = np.ones(convolution_window) / convolution_window

    positions = []

    for i, s in enumerate(samples):
        ax = axs[i]
        Cl = Cls[i]
        Cov = Covs[i]

        # plot coverage
        for strand in [1, 2]:
            factor = 1 if strand == 1 else -1
            avg_cov = np.convolve(Cov[strand], kernel, mode="same")
            ax.plot(avg_cov * factor, color="lightgray", zorder=-1)

            # relevant positions
            relevant = Cl[strand] > Cov[strand] * freq_threshold

            # plot clips
            mask = Cl[strand] > 0
            y = Cl[strand][mask] * factor
            x = pos[mask]
            c = relevant[mask]
            ax.scatter(x, y, marker=".", c=c, cmap="bwr")

            for p in pos[relevant]:
                positions.append(
                    {
                        "sample": s,
                        "strand": "fwd" if strand == 1 else "rev",
                        "position": p,
                        "coverage": Cov[strand][p],
                        "clip": Cl[strand][p],
                    }
                )

        ax.set_title(s)
        ax.axhline(0, color="black", linewidth=0.5, ls="--")
        ax.grid(alpha=0.2)
    axs[-1].set_xlabel("position (bp)")

    # custom legend
    handles = [
        mpl.lines.Line2D(
            [],
            [],
            color="lightgray",
            linestyle="-",
            label=f"coverage ({convolution_window} bp window)",
        ),
    ]
    fig.legend(handles=handles, loc="upper right")

    plt.tight_layout()
    plt.savefig(savename, dpi=300)
    plt.close(fig)

    return pd.DataFrame(positions)


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
    ax.set_xlabel("clip frequency")
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
    """Plot summary of clip frequencies and deltaF"""
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
        ax.set_ylabel("clip frequency")
    plt.tight_layout()
    plt.savefig(savename)
    plt.close(fig)


def plotly_clips(F, samples, freq_thr, savename):
    """Plot interactive non-consensus distribution with plotly"""
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
        title=f"distribution of clip sites - frequency threshold = {freq_thr}",
        xaxis_title="position (bp)",
        yaxis_title="prob. density",
        legend_title="samples",
    )
    fig.write_html(savename)


if __name__ == "__main__":
    args = parse_args()

    out_fld = pathlib.Path(args.plot_fld)
    out_fld.mkdir()
    out_summary = out_fld / "clips_summary.png"
    out_html = out_fld / "clips.html"
    out_n_clips_csv = out_fld / "n_clips.csv"
    out_n_clips = out_fld / "n_clips.png"
    out_trajs = out_fld / "clip_trajs.png"
    out_trajs_csv = out_fld / "clip_trajs.csv"

    Cls, Covs, samples = load_data(args.clips_npz, args.cov_npz)
    assert set(samples) == set(args.sample_order), "samples != sample_order"
    samples = args.sample_order

    Covs = coverage_comparison(Covs)
    Fs = safe_division(Cls, Covs, extra=np.nan)

    df = n_clip_scatterplots(
        Cls, Covs, samples, args.cov_window, args.cov_fraction, savename=out_n_clips
    )
    df.to_csv(out_n_clips_csv, index=False)

    dF = relevant_deltaF(
        Fs,
        Covs,
        freq_thr=args.freq_thr,
        cov_thr=args.cov_thr,
        noise_thr=args.noise_thr,
        noise_tol=args.noise_tol,
        max_fi=args.max_initial_freq,
    )

    plot_summary(Fs, samples, args.freq_thr, dF, savename=out_summary)

    # create a dataframe with only relevant positions (deltaF >= 0)
    df = deltaF_to_dataframe(dF, Fs, samples)

    if len(df) > 0:
        # order clusters by max frequency, and select groups of trajectories
        idxs = df.index.values[: args.n_top_trajs]

        plot_traj(
            Fs,
            Covs,
            dF,
            samples,
            idxs,
            freq_thr=args.freq_thr,
            f0_thr=args.max_initial_freq,
            cov_thr=args.cov_thr,
            savename=out_trajs,
        )

        plotly_clips(Fs, samples, args.freq_thr, savename=out_html)

        # export dataframe
        df.to_csv(out_trajs_csv, index=False)
