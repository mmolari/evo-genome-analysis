import argparse
import numpy as np
from itertools import cycle
import matplotlib as mpl
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage", type=str)
    parser.add_argument("--maxbins", type=int)
    parser.add_argument("--out_png", type=str)
    parser.add_argument("--out_html", type=str)
    return parser.parse_args()


def step_size(L, maxbin):
    step = 1
    factor = cycle([2, 5])
    while L // step > maxbin:
        step *= next(factor)
    return step


def load_coverage_dict(fname):
    "Load the coverage from all samples in a dictionary"
    cov_dict = {}
    with np.load(fname) as f:
        for sample, cov in f.items():
            # (2, L) -> L
            cov_dict[sample] = cov.sum(axis=0)

    return cov_dict


def average_coverage(cov, step):
    return np.array(
        [
            np.mean(cov[i * step : (i + 1) * step])
            for i in range(cov.size // step + 1 * bool(cov.size % step))
        ]
    )


def plot_png(cov_dict, step, fig_name):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1, 4]})

    cov_max = max([cov.max() for cov in cov_dict.values()])
    bins = np.linspace(0, cov_max, 1000)
    colors = {sname: f"C{i}" for i, sname in enumerate(cov_dict.keys())}

    # plot coverage cumulative distribution
    ax = axs[0]
    for sname, cov in cov_dict.items():
        ax.hist(
            cov,
            bins=bins,
            cumulative=True,
            histtype="step",
            color=colors[sname],
            density=True,
        )
        ax.axvline(np.median(cov), linestyle="--", color=colors[sname])
    ax.grid(alpha=0.3)
    ax.set_xlim(0, cov_max)
    ax.set_xlabel("coverage")
    ax.set_ylabel("c.d.f.")
    ax.set_title("coverage cumulative distribution")

    # plot coverage over genome
    ax = axs[1]
    for sname, cov in cov_dict.items():
        # average coverage over windows
        avg_cov = average_coverage(cov, step)
        x = np.arange(avg_cov.size) * step
        ax.plot(x, avg_cov / avg_cov.mean(), label=sname, color=colors[sname])
    ax.grid(alpha=0.3)

    ax.set_xlabel("reference position (bp)")
    ax.set_ylabel("coverage / avg. coverage")
    ax.set_title("coverage / avg. coverage distribution")
    ax.legend()
    plt.tight_layout()
    plt.savefig(fig_name)
    plt.close(fig)


def plot_html(cov_dict, step, fig_name):
    fig = make_subplots(
        rows=1,
        cols=2,
        specs=[[{"type": "xy"}, {"type": "xy"}]],
        column_widths=[0.2, 0.8],
        subplot_titles=(
            "coverage cum. distr.",
            f"coverage / avg cov ({step} bp window)",
        ),
    )

    cmap = mpl.colormaps["tab10"]
    colors = {
        sname: mpl.colors.to_hex(cmap(i)) for i, sname in enumerate(cov_dict.keys())
    }

    cov_max = max([cov.max() for cov in cov_dict.values()])
    bins = np.linspace(0, cov_max, 1000)

    for sample, cov in cov_dict.items():

        # plot coverage cumulative distribution
        y = np.histogram(cov, bins=bins)[0].cumsum() / cov.size
        x = bins[:-1]
        fig.add_trace(
            go.Scatter(
                x=x, y=y, name=sample, legendgroup=sample, line_color=colors[sample]
            ),
            row=1,
            col=1,
        )

        # plot coverage over genome
        avg_cov = average_coverage(cov, step)
        x = np.arange(avg_cov.size) * step
        y = avg_cov / avg_cov.mean()
        fig.add_trace(
            go.Scatter(
                x=x, y=y, name="", legendgroup=sample, line_color=colors[sample]
            ),
            row=1,
            col=2,
        )

    fig["layout"]["xaxis"]["title"] = "coverage"
    fig["layout"]["yaxis"]["title"] = "cdf"
    fig["layout"]["xaxis2"]["title"] = "reference position (bp)"
    fig["layout"]["yaxis2"]["title"] = "coverage / avg. coverage"

    fig.write_html(fig_name)


if __name__ == "__main__":

    args = parse_args()

    # load coverage
    cov_dict = load_coverage_dict(fname=args.coverage)
    ref_L = cov_dict[list(cov_dict.keys())[0]].size
    step = step_size(ref_L, args.maxbins)
    print(f"{ref_L=}, {args.maxbins=}, {step=}")

    # plot
    plot_png(cov_dict, step, args.out_png)
    plot_html(cov_dict, step, args.out_html)
