# %%

import argparse
import pickle as pkl
import gzip
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# %%

ins_file = (
    "../test_data/output-test-leo/pileup/reference_ST131I/R2_127_kbp/insertions.npz"
)

# load insertions
Is = {}
with np.load(ins_file) as f:
    for s in f:
        Is[s] = f[s]
# %%

for s, I in Is.items():
    plt.plot(I[0, :], marker=".", ls="")
    plt.plot(-I[1, :], marker=".", ls="")
plt.show()



# %%

Max = max([I[:2].max() for I in Is.values()])
for s, I in Is.items():
    bins = np.arange(0, Max+2, 1) - 0.5
    plt.hist(I[:2].sum(axis=0), bins=bins, histtype="step", label=s)
plt.yscale("log")
plt.show()
# %%
