import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns
import pathlib


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--in_csv", type=str, help="CSV file with non-primary alignments"
    )
    parser.add_argument("--plot_fld", type=str, help="output plot folder")
    parser.add_argument(
        "--L_thr",
        type=int,
        help="minimum alignment length on primary/supplementary reference",
    )
    return parser.parse_args()


def links(prim, suppl):
    """Finds link between positions of primary and secondary mappings."""

    # primary before suppl on query
    prim_first = np.mean([prim["qs"], prim["qe"]]) < np.mean([suppl["qs"], suppl["qe"]])

    # same strandedness
    strand = prim["fwd"] == suppl["fwd"]

    # link point
    if prim_first:
        x = prim["re"] if prim["fwd"] else prim["rs"]
        y = suppl["rs"] if suppl["fwd"] else suppl["re"]
    else:
        x = prim["rs"] if prim["fwd"] else prim["re"]
        y = suppl["re"] if suppl["fwd"] else suppl["rs"]
    return x, y, strand


def build_plot_df(df_suppl, df_prim, L_thr):
    plot_df = []
    for idx, suppl in df_suppl.iterrows():
        qry = suppl["qry"]  # query name
        prim = df_prim.loc[qry]  # corresponding primary mapping
        ref_suppl = suppl["ref"]  # reference for supplementary mapping
        ref_prim = prim["ref"]  # reference for primary mapping

        # skip if primary or supplementary mappings are too short
        if (prim["ref_len"] < L_thr) or (suppl["ref_len"] < L_thr):
            continue

        x, y, s = links(prim, suppl)

        plot_df.append(
            {
                "primary": ref_prim,
                "supplementary": ref_suppl,
                "pos. primary": x,
                "pos. suppl.": y,
                "strand": "same" if s else "opposite",
                "qry": qry,
            }
        )
    return pd.DataFrame(plot_df)


def seaborn_plot(plot_df, refs, savename):
    grid = sns.FacetGrid(
        plot_df,
        col="primary",
        row="supplementary",
        hue="strand",
        sharex="col",
        sharey="row",
        row_order=refs,
        col_order=refs,
        height=4,
        legend_out=True,
        hue_order=["same", "opposite"],
    )
    grid.map(plt.scatter, "pos. primary", "pos. suppl.", alpha=0.1)
    for (row_val, col_val), ax in grid.axes_dict.items():
        ax.grid(True, alpha=0.2)
    grid.add_legend()
    plt.savefig(savename)


def plolty_plot(plot_df, refs, savename):
    fig = px.scatter(
        plot_df,
        x="pos. primary",
        y="pos. suppl.",
        color="strand",
        facet_col="primary",
        facet_row="supplementary",
        category_orders={"primary": refs, "supplementary": refs},
        opacity=0.1,
    )
    for i in range(len(refs)):
        for j in range(len(refs)):
            fig.update_xaxes(matches=f"x{j+1}", row=i + 1, col=j + 1)
            fig.update_yaxes(matches=f"y{len(refs)*i+1}", row=i + 1, col=j + 1)
    fig.write_html(savename)


if __name__ == "__main__":
    args = parse_args()

    # output files
    out_fld = pathlib.Path(args.plot_fld)
    out_fld.mkdir()
    out_sns = out_fld / "supplementary_alignments.png"
    out_html = out_fld / "supplementary_alignments.html"

    # import read datafrane
    df = pd.read_csv(args.in_csv)

    # separate in primary and supplementary alignments
    df_suppl = df[df["suppl"]]
    df_prim = df[(~df["sec"]) & (~df["suppl"])].set_index("qry")

    # list of references
    refs = np.sort(df["ref"].unique())

    # build plot dataframe
    plot_df = build_plot_df(df_suppl, df_prim, args.L_thr)

    if len(plot_df) > 0:
        # seaborn plot
        seaborn_plot(plot_df, refs, out_sns)

        # plotly plot
        plolty_plot(plot_df, refs, out_html)
