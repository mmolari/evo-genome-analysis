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
        help="minimum alignment length on primary/secondary reference",
    )
    return parser.parse_args()


def links(prim, sec):
    """Finds link between positions of primary and secondary mappings."""

    # same strandedness
    strand = prim["fwd"] == sec["fwd"]

    # link point
    x = np.mean([prim["rs"], prim["re"]])
    y = np.mean([sec["rs"], sec["re"]])
    return x, y, strand


def build_plot_df(df_sec, df_prim, L_thr):
    plot_df = []
    for idx, sec in df_sec.iterrows():
        qry = sec["qry"]  # query name
        prim = df_prim.loc[qry]  # corresponding primary mapping
        ref_sec = sec["ref"]  # reference for secondary mapping
        ref_prim = prim["ref"]  # reference for primary mapping

        # skip if primary or secondary mappings are too short
        if (prim["ref_len"] < L_thr) or (sec["ref_len"] < L_thr):
            continue

        x, y, s = links(prim, sec)

        plot_df.append(
            {
                "primary": ref_prim,
                "secondary": ref_sec,
                "pos. primary": x,
                "pos. secondary": y,
                "strand": "same" if s else "opposite",
                "qry": qry,
            }
        )
    return pd.DataFrame(plot_df)


def seaborn_plot(plot_df, refs, savename):
    grid = sns.FacetGrid(
        plot_df,
        col="primary",
        row="secondary",
        hue="strand",
        sharex="col",
        sharey="row",
        row_order=refs,
        col_order=refs,
        height=4,
        legend_out=True,
        hue_order=["same", "opposite"],
    )
    grid.map(plt.scatter, "pos. primary", "pos. secondary", alpha=0.1)
    for (row_val, col_val), ax in grid.axes_dict.items():
        ax.grid(True, alpha=0.2)
    grid.add_legend()
    plt.savefig(savename)


def plolty_plot(plot_df, refs, savename):
    fig = px.scatter(
        plot_df,
        x="pos. primary",
        y="pos. secondary",
        color="strand",
        facet_col="primary",
        facet_row="secondary",
        category_orders={"primary": refs, "secondary": refs},
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
    out_sns = out_fld / "secondary_alignments.png"
    out_html = out_fld / "secondary_alignments.html"

    # import read datafrane
    df = pd.read_csv(args.in_csv)

    # separate in primary and secondary alignments
    df_sec = df[df["sec"]]
    df_prim = df[(~df["sec"]) & (~df["suppl"])].set_index("qry")

    # list of references
    refs = np.sort(df["ref"].unique())

    # build plot dataframe
    plot_df = build_plot_df(df_sec, df_prim, args.L_thr)

    if len(plot_df) > 0:
        # seaborn plot
        seaborn_plot(plot_df, refs, out_sns)

        # plotly plot
        plolty_plot(plot_df, refs, out_html)
