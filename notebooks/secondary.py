# %%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns


# %%
non_prim_file = "../test_data/output-test-leo/pileup/reference_ST131I/non_primary/sample1/non_primary.csv"

df = pd.read_csv(non_prim_file)
df.head()

# %%

df_suppl = df[df["suppl"]]
df_prim = df[(~df["sec"]) & (~df["suppl"])].set_index("qry")
refs = np.sort(df["ref"].unique())

# %%


def links(prim, suppl):
    pqs, pqe = prim["qs"], prim["qe"]
    sqs, sqe = suppl["qs"], suppl["qe"]

    # primary -> suppl on query
    prim_first = np.mean([pqs, pqe]) < np.mean([sqs, sqe])

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


# %%
plot_df = []
for idx, row in df_suppl.iterrows():
    qry = row["qry"]
    prim = df_prim.loc[qry]
    ref = row["ref"]
    ref_prim = prim["ref"]

    if (prim["ref_len"] < 100) or (row["ref_len"] < 100):
        continue

    x, y, s = links(prim, row)

    plot_df.append(
        {
            "primary": ref_prim,
            "supplementary": ref,
            "pos. primary": x,
            "pos. suppl.": y,
            "strand": "same" if s else "opposite",
        }
    )
plot_df = pd.DataFrame(plot_df)

# %%

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
fig.write_html("supplementary_alignments.html")
# %%

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
)
grid.map(plt.scatter, "pos. primary", "pos. suppl.", alpha=0.1)
for (row_val, col_val), ax in grid.axes_dict.items():
    ax.grid(True, alpha=0.2)
grid.add_legend()
# grid.fig.tight_layout()
plt.show()

# %%
