import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=str, nargs="+")
    parser.add_argument("--step", type=int)
    parser.add_argument("--out_png",type=str)
    parser.add_argument("--out_html",type=str)
    return parser.parse_args()

def sample_name(fname):
    return fname.split("/")[-1]

def load_coverage(fname):
    with np.load(fname) as f:
        # shape: (2, 6, L)
        pu = f["pileup"]

    # nucl order: ["A", "C", "G", "T", "-", "N"]
    cov = np.sum(pu, axis=0)[:4].sum(axis=0)
    return cov

def average_coverage(cov, step):
    return np.array([
        np.mean(cov[i*step:(i+1)*step]) for i in range(cov.size//step + 1*bool(cov.size%step))
    ])

def plot_png(cov_dict, avg_cov_dict, fig_name):
    fig, axs = plt.subplots(1,2, figsize=(12,5), gridspec_kw={"width_ratios": [1, 4]})

    cov_max = max([cov.max() for cov in cov_dict.values()])
    bins = np.linspace(0, cov_max, 1000)
    colors = {sname: f"C{i}" for i, sname in enumerate(cov_dict.keys())}

    # plot coverage cumulative distribution
    ax = axs[0]
    for sname, cov in cov_dict.items():
        ax.hist(cov, bins=bins, cumulative=True, histtype="step", color=colors[sname])
        ax.axvline(cov.median(), linestyle="--", color=colors[sname])
    ax.grid()

    # plot coverage over genome
    ax = axs[1]
    for sname, cov in avg_cov_dict.items():
        x = np.arange(cov.size)*args.step
        ax.plot(x, cov/cov.mean(), label=sname, color=colors[sname])
    ax.grid()

    ax.set_xlabel("position (bp)")
    ax.set_ylabel("coverage / avg. coverage")
    plt.legend()
    plt.savefig(fig_name)
    plt.close(fig)

if __name__ == "__main__":
    
    args = parse_args()

    # load coverage
    cov_dict = { sample_name(fname): load_coverage(fname) for fname in args.samples}
    avg_cov_dict = {sname: average_coverage(cov, args.step) for sname, cov in cov_dict.items()}

    # plot
    plot_png(cov_dict, avg_cov_dict, args.out_png)




