# %%
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import plotly.subplots as subp
import plotly.graph_objects as go


# %%
step = 5000

cov_dict = {}

for i in range(1, 5):
    cov_file = f"../test_data/output-test-leo/pileup/reference_ST131I/R1_4963_kbp/sample{i}/allele_counts.npz"

    # load numpy array
    with np.load(cov_file) as f:
        # shape: (2, 6, L)
        pu = f["pileup"]

    # nucl order: ["A", "C", "G", "T", "-", "N"]
    cov = np.sum(pu, axis=0)[:4].sum(axis=0)
    cov_dict[f"sample{i}"] = cov

# %%
def average_coverage(cov, step):
    return np.array(
        [
            np.mean(cov[i * step : (i + 1) * step])
            for i in range(cov.size // step + 1 * bool(cov.size % step))
        ]
    )


def plot_png(cov_dict, step, fig_name):
    fig, axs = plt.subplots(1, 2, figsize=(15, 4), gridspec_kw={"width_ratios": [1, 4]})

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

    # plot coverage over genome
    ax = axs[1]
    for sname, cov in cov_dict.items():
        # average coverage over windows
        avg_cov = average_coverage(cov, step)
        x = np.arange(avg_cov.size) * step
        ax.plot(x, avg_cov / avg_cov.mean(), label=sname, color=colors[sname])
    ax.grid(alpha=0.3)

    ax.set_xlabel("position (bp)")
    ax.set_ylabel("coverage / avg. coverage")
    plt.legend()
    plt.savefig(fig_name)
    plt.show(fig)


plot_png(cov_dict, step, "test_fig/coverage.png")

# %%


fig = go.Figure()
for i in range(1, 5):
    x = np.arange(cov_dict[i].size) * step
    y = cov_dict[i] / cov_dict[i].mean()
    fig.add_trace(go.Scatter(x=x, y=y, name=f"sample {i}", visible=True))

fig.update_layout(
    title=f"Coverage ({step} bp window)",
    xaxis_title="position (bp)",
    yaxis_title="coverage / avg. coverage",
    legend_title="sample",
)
fig.write_html("test_fig/coverage.html")

# %%


def plot_html(cov_dict, step, fig_name):
    fig = subp.make_subplots(
        rows=1,
        cols=2,
        specs=[[{"type": "xy"}, {"type": "xy"}]],
        column_widths=[0.2, 0.8],
        subplot_titles=(
            "coverage cum. distr.",
            f"(coverage / avg cov) vs genome ({step} bp window)",
        ),
    )

    cmap = mpl.colormaps["tab10"]
    colors = {
        sname: mpl.colors.to_hex(cmap(i)) for i, sname in enumerate(cov_dict.keys())
    }

    cov_max = max([cov.max() for cov in cov_dict.values()])
    bins = np.linspace(0, cov_max, 1000)

    for sample, cov in cov_dict.items():

        # cumulative distr
        y = np.histogram(cov, bins=bins)[0].cumsum() / cov.size
        x = bins[:-1]

        # plot coverage cumulative distribution
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

    fig['layout']['xaxis']['title']='coverage'
    fig['layout']['yaxis']['title']='cdf'
    fig['layout']['xaxis2']['title']='genome position (bp)'
    fig['layout']['yaxis2']['title']='coverage / avg. coverage'


    fig.write_html(fig_name)


fig = plot_html(cov_dict, step, "test_fig/coverage.html")

# %%

# 191
