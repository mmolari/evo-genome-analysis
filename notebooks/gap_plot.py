# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def load_npz_dict(fname):
    npz_dict = {}
    with np.load(fname) as f:
        for sample, data in f.items():
            npz_dict[sample] = data
    return npz_dict


def safe_division(a, b, extra):
    "return the division between the two vectors. Return nan if b is zero"
    # create a mask to identify zero values in the denominator
    mask = b == 0

    # perform the division, replacing zero denominators with NaN
    result = np.divide(a, b, out=np.zeros_like(a, dtype=float), where=~mask)
    result[mask] = extra
    return result


# %%

gap_file = "../test_data/output-test-leo/pileup/reference_ST131I/R1_4963_kbp/gaps.npz"
cons_file = (
    "../test_data/output-test-leo/pileup/reference_ST131I/R1_4963_kbp/coverage.npz"
)

gaps = load_npz_dict(gap_file)
cons = load_npz_dict(cons_file)
samples = list(gaps.keys())
# %%

# gap frequencies
bins = np.linspace(0, 1, 100)
for s in samples:
    g = gaps[s].sum(axis=0)
    c = cons[s].sum(axis=0)
    freq = safe_division(g, c + g, extra=np.nan)
    plt.hist(freq, histtype="step", bins=bins, label=s)
plt.legend()
plt.yscale("log")
plt.show()
# %%
S = len(samples)
L = gaps[samples[0]].shape[1]
bins = np.arange(0, L, 5000)
fig, axs = plt.subplots(S, 1, figsize=(10, S * 3), sharex=True)
for ns, s in enumerate(samples):
    ax = axs[ns]
    for i in range(2):
        g = gaps[s][i, :]
        c = cons[s][i, :]
        freq = safe_division(g, c + g, extra=0)
        print(freq.shape)
        ax.hist(np.arange(L), weights=freq, histtype="step", bins=bins)
plt.show()
# %%
thr = 0.85
S = len(samples)
L = gaps[samples[0]].shape[1]
bins = np.arange(0, L, 5000)
fig, ax = plt.subplots(1, 1, figsize=(10, 3))
for ns, s in enumerate(samples):
    g = gaps[s]
    c = cons[s]
    freq = safe_division(g, c + g, extra=0)
    freq = freq.min(axis=0)
    ticks = np.argwhere(freq > thr).flatten()
    ax.plot(ticks, np.ones_like(ticks) * ns, marker="|", linestyle="None", alpha=0.3)
plt.tight_layout()
plt.show()

#  select and plot relevant trajectories (positions plus variation)

# %%
