#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from liblibra_core import *
from libra_py import units, data_stat
import libra_py.packages.cp2k.methods as CP2K_methods

get_ipython().run_line_magic('matplotlib', 'notebook')

systems = {
    "100 K": {
        "path": "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/0500_100K/step3/res-sd-30-occ",
        "color": "tab:blue"
    },
    "300 K": {
        "path": "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/TiO2_5/H2O_0/step3_NACs_CAM_B3LYP/best/300K/AS_30_40",
        "color": "tab:orange"
    }
}

def analyze_energy_gap_0_1(systems, temperature):

    plt.figure(figsize=(8,5))

    for label, info in systems.items():

        params = {
            "path_to_energy_files": info["path"],
            "dt": 1.0,
            "prefix": "Hvib_ci_",
            "suffix": "_re",
            "istep": 1000,
            "fstep": 3990
        }

        try:
            md_time, energies = CP2K_methods.extract_energies_sparse(params)
        except Exception as e:
            print(f"[Error] Couldn't extract energies for {label}: {e}")
            continue

        energies = energies * units.au2ev

        if energies.shape[1] < 2:
            print(f"[Skip] Not enough states for {label}")
            continue

        gap = np.abs(energies[:, 1] - energies[:, 0])

        energy_gaps = []
        for val in gap:
            m = MATRIX(1, 1)
            m.set(0, 0, val)
            energy_gaps.append(m)

        bin_supp, dens, cum = data_stat.cmat_distrib(energy_gaps, 0, 0, 0, 0, 50, 0.01)

        plt.plot(
            bin_supp,
            dens,
            label=label,
            color=info["color"],
            linewidth=2,
            alpha=0.6
        )

    plt.xlabel('Energy gap (S$_1$ − S$_0$) (eV)', fontsize=26)
    plt.ylabel('PD (1/eV)', fontsize=26)
    plt.yscale('log')
    plt.xlim(2, 5)
    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.legend(fontsize=22)
    plt.tight_layout()
    plt.savefig(f'energy_gap_0_1_{temperature}K.jpg', dpi=600)
    plt.show()

analyze_energy_gap_0_1(systems, 300)


# In[ ]:


import os
import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

from liblibra_core import *
from libra_py import units, data_stat

get_ipython().run_line_magic('matplotlib', 'inline')

systems = {
    "100 K": {
        "path": "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/0500_100K/step3/res-sd-30-occ",
        "color": "tab:blue"
    },
    "300 K": {
        "path": "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/TiO2_5/H2O_0/step3_NACs_CAM_B3LYP/best/300K/AS_30_40",
        "color": "tab:orange"
    }
}

def analyze_nac_s0_s1s3(systems, states=[1], nbins=80, dx=0.05):

    plt.figure(figsize=(8, 5))

    for label, info in systems.items():

        folder = info["path"]
        color = info["color"]

        nac_values = []

        # Imaginary part of Hvib contains NACs
        nac_files = sorted(glob.glob(os.path.join(folder, "Hvib_ci_*_im.npz")))

        if not nac_files:
            print(f"[Warning] No NAC files found for {label} in {folder}")
            continue

        for nac_file in nac_files:

            hvib = sp.load_npz(nac_file)
            hvib_dense = np.asarray(hvib.todense()).real

            nstates = hvib_dense.shape[0]

            for j in states:
                if j < nstates:
                    # Hvib_im is in Hartree; convert to meV
                    nac_0j = abs(hvib_dense[0, j]) * units.au2ev * 1000.0

                    m = MATRIX(1, 1)
                    m.set(0, 0, nac_0j)
                    nac_values.append(m)

        if len(nac_values) == 0:
            print(f"[Skip] No valid NAC values for {label}")
            continue

        bin_supp, dens, cum = data_stat.cmat_distrib(
            nac_values,
            0, 0,
            0, 0,
            nbins,
            dx
        )

        plt.plot(
            bin_supp,
            dens,
            label=label,
            color=color,
            linewidth=2,
            alpha=0.8
        )

        print(f"{label}: collected {len(nac_values)} NAC values")

    plt.xlabel(r'|NAC(S$_0$, S$_1$)| (meV)', fontsize=26)
    plt.ylabel(r'PD (1/meV)', fontsize=26)

    plt.xscale('log')
    plt.yscale('log')

    # log scale cannot start at 0
    plt.xlim(0.05, 5.0)

    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
  

    plt.tight_layout()
    plt.savefig("nac_s0_s1_100K_300K.jpg", dpi=600, bbox_inches="tight")
    plt.show()


analyze_nac_s0_s1s3(systems)


# In[ ]:


import os
import glob
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt

from liblibra_core import *
from libra_py import units, data_stat

# =========================
# Input paths (one isomer)
# =========================
systems = {
    "100 K": {
        "path": "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/0500_100K/step3/res-sd-30-occ",
        "color": "tab:blue"
    },
    "300 K": {
        "path": "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/TiO2_5/H2O_0/step3_NACs_CAM_B3LYP/best/300K/AS_30_40",
        "color": "tab:orange"
    }
}

# =========================
# Main analysis function
# =========================
def analyze_nac_adjacent_states(systems, max_state=10, nbins=80, dx=0.05):

    plt.figure(figsize=(8, 5))

    for label, info in systems.items():

        folder = info["path"]
        color = info["color"]

        nac_values = []

        # Load imaginary Hvib matrices (NACs)
        nac_files = sorted(glob.glob(os.path.join(folder, "Hvib_ci_*_im.npz")))

        if not nac_files:
            print(f"[Warning] No NAC files found for {label} in {folder}")
            continue

        for nac_file in nac_files:

            try:
                hvib = sp.load_npz(nac_file)
                hvib_dense = np.asarray(hvib.todense()).real
            except Exception as e:
                print(f"[Skip] Failed loading {nac_file}: {e}")
                continue

            nstates = hvib_dense.shape[0]
            upper_state = min(max_state, nstates - 1)

            # Adjacent pairs: (i, i+1)
            for i in range(upper_state):
                j = i + 1

                nac_ij = abs(hvib_dense[i, j]) * units.au2ev * 1000.0  # meV

                m = MATRIX(1, 1)
                m.set(0, 0, nac_ij)
                nac_values.append(m)

        if len(nac_values) == 0:
            print(f"[Skip] No valid NAC values for {label}")
            continue

        print(f"{label}: collected {len(nac_values)} NAC values")

        # Histogram via Libra
        bin_supp, dens, cum = data_stat.cmat_distrib(
            nac_values,
            0, 0,
            0, 0,
            nbins,
            dx
        )

        plt.plot(
            bin_supp,
            dens,
            label=label,
            color=color,
            linewidth=2,
            alpha=0.85
        )

    # =========================
    # Plot styling
    # =========================
    plt.xlabel(r'|NAC(S$_i$, S$_{i+1}$)| (meV)', fontsize=26)
    plt.ylabel(r'PD (1/meV)', fontsize=26)

    plt.xscale('log')
    plt.yscale('log')

    # Avoid zero for log scale
    plt.xlim(0.05, 5.0)

    plt.xticks(fontsize=26)
    plt.yticks(fontsize=26)
    plt.legend(fontsize=22)

    plt.tight_layout()
    plt.savefig("nac_adjacent_states_100K_300K.jpg", dpi=600, bbox_inches="tight")
    plt.show()


# =========================
# Run
# =========================

# Low-energy manifold (recommended)
analyze_nac_adjacent_states(systems, max_state=3)

# If you want more states:
# analyze_nac_adjacent_states(systems, max_state=10)

