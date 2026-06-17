import os, glob, time, h5py, warnings
import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from scipy.optimize import curve_fit

from liblibra_core import *
import util.libutil as comn

import matplotlib.ticker as ticker

import libra_py
from libra_py import units, data_conv
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.workflows.nbra.decoherence_times as decoherence_times
import libra_py.data_visualize

warnings.filterwarnings('ignore')

# =========================
# Define functions
# =========================

def exp_functFSSH(t, _x):
    return 1 - np.exp(-np.power(t / _x, 1.3))

def exp_functMSDM(t, _x):
    return 1 - np.exp(-np.power(t / _x, 1.1))


ICONDS = list(range(1, 3000, 100))

# =========================
# Figure
# =========================

fig, ax = plt.subplots(figsize=(10, 6.5))

methods = ['FSSH', 'MSDM']
colors = ['red', 'blue']

timescale_labels = []
all_times = []

# =========================
# Loop over methods
# =========================

for idx, method in enumerate(methods):

    taus = []

    if method == 'FSSH':
        exp_func = exp_functFSSH
    elif method == 'MSDM':
        exp_func = exp_functMSDM
    else:
        raise ValueError(f"Unknown method: {method}")

    for icond in ICONDS:

        try:
            with h5py.File(f'{method}_icond_{icond}/mem_data.hdf', 'r') as F:

                sh_pop = np.array(F['sh_pop_adi/data'][:, 0])
                md_time = np.array(F['time/data'][:]) * units.au2fs

            all_times.append(md_time)

            # Background trajectories
            ax.plot(
                md_time,
                sh_pop,
                alpha=0.10,
                color=colors[idx],
                linewidth=1.0
            )

            # Fit
            popt, pcov = curve_fit(
                exp_func,
                md_time,
                sh_pop,
                bounds=([0.0], [np.inf])
            )

            _tau = popt[0]

            # R-squared
            residuals = sh_pop - exp_func(md_time, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((sh_pop - np.mean(sh_pop))**2)

            if ss_tot > 0:
                r_squared = 1.0 - (ss_res / ss_tot)
            else:
                r_squared = np.nan

            print(f"Method: {method}, IC {icond}, R2 = {r_squared}")

            if np.isfinite(r_squared) and r_squared > 0.0:
                taus.append(_tau)

        except Exception as e:
            print(f"Error processing {method}, IC {icond}: {e}")
            continue

    # =========================
    # Statistics
    # =========================

    taus = np.array(taus)

    if len(taus) == 0:
        print(f"No valid fits for {method}")
        continue

    ave_tau = np.average(taus)
    s = np.std(taus)

    Z = 1.96
    N = taus.shape[0]

    error_bar = Z * s / np.sqrt(N)

    ave_tau_ps = ave_tau / 1000.0
    error_bar_ps = error_bar / 1000.0

    print(f'Timescales for {method}: {ave_tau_ps:.3f} ± {error_bar_ps:.3f} ps')

    timescale_labels.append(
        f"{method}: {ave_tau_ps:.1f} ± {error_bar_ps:.1f} ps"
    )

    # =========================
    # Plot fitted curve
    # =========================

    if len(all_times) > 0:
        t_fit = np.linspace(0.0, max([np.max(t) for t in all_times]), 1000)
    else:
        t_fit = md_time

    ax.plot(
        t_fit,
        exp_func(t_fit, ave_tau),
        linewidth=3,
        color=colors[idx],
        label=method
    )

# =========================
# Axis formatting
# =========================

if len(all_times) > 0:
    xmax = max([np.max(t) for t in all_times])
else:
    xmax = 6000.0

ax.set_xlim(0.0, xmax)

ax.set_xlabel("Time (fs)", fontsize=30)
ax.set_ylabel("Ground State Population", fontsize=30)

ax.tick_params(axis='both', which='major', labelsize=28)

# Remove automatic scientific notation/offset if present
ax.yaxis.get_offset_text().set_visible(False)

# Spine thickness
for spine in ax.spines.values():
    spine.set_linewidth(1.5)

# =========================
# Legend
# =========================

ax.legend(fontsize=28)

# =========================
# Optional timescale box
# =========================
# Uncomment if you want the fitted times written inside the plot.

# timescale_text = "\n".join(timescale_labels)
# props = dict(boxstyle='round', facecolor='white', alpha=0.8)
# ax.text(
#     0.05, 0.95,
#     timescale_text,
#     transform=ax.transAxes,
#     fontsize=20,
#     verticalalignment='top',
#     bbox=props
# )

# =========================
# Save & show
# =========================

plt.tight_layout()

plt.savefig(
    "newTiO2_5_100K_GS_pop_evolution_FSSH_MSDM.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    "newTiO2_5_100K_GS_pop_evolution_FSSH_MSDM.pdf",
    bbox_inches="tight"
)

plt.show()

