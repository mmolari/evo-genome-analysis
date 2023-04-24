# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# %%
step = 5000

cov_dict = {}

for i in range(1,5):
    cov_file = f"../test_data/output-test-leo/pileup/reference_ST131I/R1_4963_kbp/sample{i}/allele_counts.npz"

    # load numpy array
    with np.load(cov_file) as f:
        # shape: (2, 6, L)
        pu = f["pileup"]

    # nucl order: ["A", "C", "G", "T", "-", "N"]
    cov = np.sum(pu, axis=0)[:4].sum(axis=0)
    cov_dict[i] = np.array([
        np.mean(cov[i*step:(i+1)*step]) for i in range(cov.size//step + 1*bool(cov.size%step))
    ])

# %%

fig, ax = plt.subplots(1,1, figsize=(12,5))
for i in range(1,5):
    x = np.arange(cov_dict[i].size)*step
    ax.plot(x, cov_dict[i]/cov_dict[i].mean(), label=f"sample {i}")
ax.grid()
ax.set_xlabel("position (bp)")
plt.show()

# %%
import plotly.express as px
import plotly.graph_objects as go

fig = go.Figure()
for i in range(1,5):
    x = np.arange(cov_dict[i].size)*step
    y = cov_dict[i]/cov_dict[i].mean()
    fig.add_trace(go.Scatter(x=x, y=y, name=f"sample {i}", visible=True))

fig.update_layout(
    title=f"Coverage ({step} bp window)",
    xaxis_title="position (bp)",
    yaxis_title="coverage / avg. coverage",
    legend_title="sample",
)
fig.write_html("test_fig/coverage.html")
# %%
