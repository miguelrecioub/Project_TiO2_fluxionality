#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.optimize import curve_fit
from scipy.ndimage import uniform_filter1d
from libra_py.units import au2ev, au2fs
get_ipython().run_line_magic('matplotlib', 'inline')

# ----------------------------
# Stretched exponential decay with fixed E0 and E_inf
# ----------------------------
def stretched_exponential_fixed_bounds(t, tau, beta, E0, Einf):
    return Einf + (E0 - Einf) * np.exp(-(t / tau) ** beta)

# ----------------------------
# Load Adiabatic Energies
# ----------------------------
def load_adiabatic_energies(start_step, end_step, path_template, scale=au2ev, time_scale=au2fs):
    energies = []
    time = None
    for step in range(start_step, end_step):
        path = path_template.format(step=step)
        energy = sp.load_npz(path).todense().real
        if time is None:
            time = np.arange(energy.shape[0]) * time_scale
        energies.append(np.diag(energy))
    return np.array(energies) * scale, time

# ----------------------------
# Main Analysis Function
# ----------------------------
def calculate_total_energy_decay_to_S1(
    adiabatic_energies, methods, icond_range, method_config,
    smoothing_window=100, fit_time_cutoff=3000, verbose=True
):
    energy_results = {}
    average_energies = {}
    time_results = {}
    decay_fit_results = {}

    # --- Compute population-weighted average S1 energy ---
    all_pops = []
    all_E1 = []

    for method in methods:
        path_template = method_config[method]["path_template"]
        for icond in icond_range:
            with h5py.File(path_template.format(icond=icond), 'r') as f:
                sh_pop = np.array(f['sh_pop_adi/data'])
                n_steps = min(sh_pop.shape[0], adiabatic_energies.shape[0])
                pops = sh_pop[:n_steps, 1]  # S1 population
                E1 = adiabatic_energies[:n_steps, 2]  # S1 energy
                all_pops.append(pops)
                all_E1.append(E1)

    all_pops = np.array(all_pops)
    all_E1 = np.array(all_E1)
    E1_avg = np.sum(all_pops * all_E1) / np.sum(all_pops)

    for method in methods:
        print(f"\nProcessing method: {method}")
        try:
            method_results = []
            method_times = None
            path_template = method_config[method]["path_template"]

            initial_total_energies = []

            for icond in icond_range:
                path = path_template.format(icond=icond)
                with h5py.File(path, 'r') as f:
                    sh_pop = np.array(f['sh_pop_adi/data'])
                    time = np.array(f['time/data']) * au2fs

                    n_steps = min(sh_pop.shape[0], adiabatic_energies.shape[0])
                    sh_pop = sh_pop[:n_steps, :]
                    time = time[:n_steps]
                    Hvib_rolled = np.roll(adiabatic_energies[:n_steps, :], -icond, axis=0)

                    total_energy = np.sum(sh_pop * Hvib_rolled, axis=1)
                    method_results.append(total_energy)
                    initial_total_energies.append(total_energy[0])

                    if method_times is None:
                        method_times = time

            total_energy_avg = np.mean(method_results, axis=0)
            total_energy_std = np.std(method_results, axis=0)
            smoothed_total = uniform_filter1d(total_energy_avg, size=smoothing_window)

            E0_fixed = np.mean(initial_total_energies)
            Einf_fixed = E1_avg

            fit_mask = method_times < fit_time_cutoff
            fit_times = method_times[fit_mask]
            fit_values = smoothed_total[fit_mask]

            def fit_func_stretched(t, tau, beta):
                return stretched_exponential_fixed_bounds(t, tau, beta, E0_fixed, Einf_fixed)

            popt, pcov = curve_fit(
                fit_func_stretched,
                fit_times,
                fit_values,
                p0=[300, 0.8],
                bounds=([1e-2, 0.01], [1e4, 1.0])
            )

            tau_fit, beta_fit = popt
            tau_std, beta_std = np.sqrt(np.diag(pcov))
            fit_curve = fit_func_stretched(method_times, tau_fit, beta_fit)

            if verbose:
                plt.figure()
                plt.plot(method_times, total_energy_avg, label='Raw Total Energy', alpha=0.4)
                plt.fill_between(
                    method_times,
                    total_energy_avg - total_energy_std,
                    total_energy_avg + total_energy_std,
                    color='gray', alpha=0.3, label='±1σ'
                )
                plt.plot(method_times, smoothed_total, label='Smoothed Total Energy', color='red')
                plt.plot(method_times, fit_curve, '--',
                         label=f'Stretched Fit (t<{fit_time_cutoff} fs): τ = {tau_fit:.2f} ± {tau_std:.2f} fs, β = {beta_fit:.2f} ± {beta_std:.2f}',
                         color='black')
                plt.axhline(E0_fixed, color='green', linestyle=':', label=f'Initial ≈ {E0_fixed:.3f} eV')
                plt.axhline(Einf_fixed, color='blue', linestyle=':', label=f'S₁ Energy ≈ {Einf_fixed:.3f} eV (weighted)')

                plt.xlabel("Time (fs)", fontsize=22)
                plt.ylabel("Average Excitation Energy (eV)", fontsize=22)
                plt.tick_params(axis='both', which='major', labelsize=22)
                plt.gcf().set_size_inches(12, 6)
                plt.tight_layout(rect=[0, 0, 0.75, 1])
                plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))

                plt.show()

            print(f"Average total energy at t=0: {E0_fixed:.6f} eV")
            print(f"Fitted stretched exponential (t < {fit_time_cutoff} fs):")
            print(f"  τ = {tau_fit:.2f} ± {tau_std:.2f} fs")
            print(f"  β = {beta_fit:.3f} ± {beta_std:.3f}")
            print(f"  Einf (fixed) = {Einf_fixed:.6f} eV")

            energy_results[method] = (method_results, [tau_fit, beta_fit])
            average_energies[method] = smoothed_total
            time_results[method] = method_times
            decay_fit_results[method] = {
                "fit_curve": fit_curve,
                "params": (tau_fit, beta_fit, Einf_fixed),
                "errors": (tau_std, beta_std, 0.0)
            }

        except Exception as e:
            print(f"Error in method {method}: {e}")

    return energy_results, average_energies, time_results, decay_fit_results


# In[2]:


# ----------------------------
# Energy file info
# ----------------------------
start_step = 1000
end_step = 3990
energy_path_template = "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/0500_300K/isomer2/step3/AS_30_39/Hvib_ci_{step}_re.npz"
get_ipython().run_line_magic('matplotlib', 'inline')

# ----------------------------
# Load adiabatic energies
# ----------------------------
adiabatic_energies, _ = load_adiabatic_energies(start_step, end_step, energy_path_template)

# ----------------------------
# Method configuration
# ----------------------------
method_config = {
    "DISH": {
        "path_template": "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/0500_300K/isomer2/step4/DISH_icond_{icond}/mem_data.hdf"
    }
}

# ----------------------------
# Initial condition indices
# ----------------------------
icond_range = range(1, 3000, 100)

# ----------------------------
# Run energy decay analysis with fixed decay to S₁
# ----------------------------
methods = ["DISH"]
total_results, total_avg, total_times, total_fits = calculate_total_energy_decay_to_S1(
    adiabatic_energies=adiabatic_energies,
    methods=methods,
    icond_range=icond_range,
    method_config=method_config,
    smoothing_window=100,
    verbose=True
)


# In[3]:


import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.optimize import curve_fit
from scipy.ndimage import uniform_filter1d
from libra_py.units import au2ev, au2fs
get_ipython().run_line_magic('matplotlib', 'inline')

# ----------------------------
# Bi-exponential decay function with fixed E_inf
# ----------------------------
def bi_exponential_decay(t, A1, tau1, A2, tau2, Einf):
    return Einf + A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)

# ----------------------------
# Load Adiabatic Energies
# ----------------------------
def load_adiabatic_energies(start_step, end_step, path_template, scale=au2ev, time_scale=au2fs):
    energies = []
    time = None
    for step in range(start_step, end_step):
        path = path_template.format(step=step)
        energy = sp.load_npz(path).todense().real
        if time is None:
            time = np.arange(energy.shape[0]) * time_scale
        energies.append(np.diag(energy))
    return np.array(energies) * scale, time

# ----------------------------
# Main Analysis Function
# ----------------------------
def calculate_total_energy_decay_to_S1(
    adiabatic_energies, methods, icond_range, method_config,
    smoothing_window=100, fit_time_cutoff=1000, verbose=True
):
    energy_results = {}
    average_energies = {}
    time_results = {}
    decay_fit_results = {}

    all_pops = []
    all_E1 = []

    # --- Compute weighted average S1 energy ---
    for method in methods:
        path_template = method_config[method]["path_template"]
        for icond in icond_range:
            with h5py.File(path_template.format(icond=icond), 'r') as f:
                sh_pop = np.array(f['sh_pop_adi/data'])
                n_steps = min(sh_pop.shape[0], adiabatic_energies.shape[0])
                pops = sh_pop[:n_steps, 1]
                E1 = adiabatic_energies[:n_steps, 1]
                all_pops.append(pops)
                all_E1.append(E1)

    all_pops = np.array(all_pops)
    all_E1 = np.array(all_E1)
    E1_avg = np.sum(all_pops * all_E1) / np.sum(all_pops)

    for method in methods:
        print(f"\nProcessing method: {method}")
        try:
            method_results = []
            method_times = None
            path_template = method_config[method]["path_template"]

            tau1_list = []
            tau2_list = []
            initial_total_energies = []

            for icond in icond_range:
                path = path_template.format(icond=icond)
                with h5py.File(path, 'r') as f:
                    sh_pop = np.array(f['sh_pop_adi/data'])
                    time = np.array(f['time/data']) * au2fs

                    n_steps = min(sh_pop.shape[0], adiabatic_energies.shape[0])
                    sh_pop = sh_pop[:n_steps, :]
                    time = time[:n_steps]
                    Hvib_rolled = np.roll(adiabatic_energies[:n_steps, :], -icond, axis=0)

                    total_energy = np.sum(sh_pop * Hvib_rolled, axis=1)
                    initial_total_energies.append(total_energy[0])

                    # Smooth individual energy
                    smoothed = uniform_filter1d(total_energy, size=smoothing_window)

                    # Fit to bi-exponential for this trajectory
                    fit_mask = time < fit_time_cutoff
                    fit_times = time[fit_mask]
                    fit_values = smoothed[fit_mask]

                    def fit_func_biexp(t, A1, tau1, A2, tau2):
                        return bi_exponential_decay(t, A1, tau1, A2, tau2, E1_avg)

                    E0_i = total_energy[0]
                    p0 = [E0_i - E1_avg, 300, (E0_i - E1_avg) * 0.5, 1000]
                    bounds = ([-np.inf, 1e-2, -np.inf, 1e-2], [np.inf, 1e4, np.inf, 1e4])

                    try:
                        popt, _ = curve_fit(fit_func_biexp, fit_times, fit_values, p0=p0, bounds=bounds)
                        A1_i, tau1_i, A2_i, tau2_i = popt
                        tau1_list.append(tau1_i)
                        tau2_list.append(tau2_i)
                    except Exception as e_fit:
                        print(f"  Fit failed for icond {icond}: {e_fit}")

                    method_results.append(total_energy)
                    if method_times is None:
                        method_times = time

            # Aggregate energy and fit summary
            total_energy_avg = np.mean(method_results, axis=0)
            total_energy_std = np.std(method_results, axis=0)
            smoothed_total = uniform_filter1d(total_energy_avg, size=smoothing_window)

            avg_tau1 = np.mean(tau1_list)
            std_tau1 = np.std(tau1_list)
            avg_tau2 = np.mean(tau2_list)
            std_tau2 = np.std(tau2_list)

            # Fit average energy for plotting (dashed black line)
            fit_mask = method_times < fit_time_cutoff
            fit_times = method_times[fit_mask]
            fit_values = smoothed_total[fit_mask]

            def fit_func_biexp(t, A1, tau1, A2, tau2):
                return bi_exponential_decay(t, A1, tau1, A2, tau2, E1_avg)

            E0_fixed = np.mean(initial_total_energies)
            p0 = [E0_fixed - E1_avg, 300, (E0_fixed - E1_avg) * 0.5, 1000]

            try:
                popt, _ = curve_fit(fit_func_biexp, fit_times, fit_values, p0=p0,
                                    bounds=([-np.inf, 1e-2, -np.inf, 1e-2], [np.inf, 1e4, np.inf, 1e4]))
                fit_curve = fit_func_biexp(method_times, *popt)
                tau1_fit, tau2_fit = popt[1], popt[3]
            except Exception as e_avg_fit:
                print(f"Average fit failed: {e_avg_fit}")
                fit_curve = None

            # Plot
            if verbose:
                plt.figure()
                plt.plot(method_times, total_energy_avg, label='Raw Total Energy', alpha=0.4)
                plt.fill_between(
                    method_times,
                    total_energy_avg - total_energy_std,
                    total_energy_avg + total_energy_std,
                    color='gray', alpha=0.3, label='±1σ'
                )
                plt.plot(method_times, smoothed_total, label='Smoothed Avg Energy', color='red')
                if fit_curve is not None:
                    plt.plot(method_times, fit_curve, '--', color='black',
                             label=f'Avg Fit: τ₁={tau1_fit:.1f} fs, τ₂={tau2_fit:.1f} fs')

                plt.axhline(E0_fixed, color='green', linestyle=':', label=f'Initial ≈ {E0_fixed:.3f} eV')
                plt.axhline(E1_avg, color='blue', linestyle=':', label=f'S₁ Energy ≈ {E1_avg:.3f} eV')

                plt.xlabel("Time (fs)", fontsize=22)
                plt.ylabel("Average Excitation Energy (eV)", fontsize=22)
                plt.tick_params(axis='both', which='major', labelsize=22)
                plt.gcf().set_size_inches(12, 6)
                plt.tight_layout(rect=[0, 0, 0.75, 1])
              #  plt.legend(fontsize=14)
                plt.show()

            print(f"τ₁ = {avg_tau1:.2f} ± {std_tau1:.2f} fs")
            print(f"τ₂ = {avg_tau2:.2f} ± {std_tau2:.2f} fs")

            energy_results[method] = (method_results, (avg_tau1, std_tau1, avg_tau2, std_tau2))
            average_energies[method] = smoothed_total
            time_results[method] = method_times
            decay_fit_results[method] = {
                "tau1_avg": avg_tau1,
                "tau1_std": std_tau1,
                "tau2_avg": avg_tau2,
                "tau2_std": std_tau2,
                "Einf": E1_avg,
                "individual_tau1": tau1_list,
                "individual_tau2": tau2_list
            }

        except Exception as e:
            print(f"Error in method {method}: {e}")

    return energy_results, average_energies, time_results, decay_fit_results


# In[4]:


# ----------------------------
# Energy file info
# ----------------------------
start_step = 1000
end_step = 3990
energy_path_template = "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/0500_300K/isomer2/step3/AS_30_39/Hvib_ci_{step}_re.npz"
get_ipython().run_line_magic('matplotlib', 'inline')

# ----------------------------
# Load adiabatic energies
# ----------------------------
adiabatic_energies, _ = load_adiabatic_energies(start_step, end_step, energy_path_template)

# ----------------------------
# Method configuration
# ----------------------------
method_config = {
    "DISH": {
        "path_template": "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/0500_300K/isomer2/step4/DISH_icond_{icond}/mem_data.hdf"
    }
}

# ----------------------------
# Initial condition indices
# ----------------------------
icond_range = range(1, 3000, 100)

# ----------------------------
# Run energy decay analysis with fixed decay to S₁
# ----------------------------
methods = ["DISH"]
total_results, total_avg, total_times, total_fits = calculate_total_energy_decay_to_S1(
    adiabatic_energies=adiabatic_energies,
    methods=methods,
    icond_range=icond_range,
    method_config=method_config,
    smoothing_window=100,
    verbose=True
)


# In[5]:


import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.optimize import curve_fit
from scipy.ndimage import uniform_filter1d
from libra_py.units import au2ev, au2fs
get_ipython().run_line_magic('matplotlib', 'inline')

# ----------------------------
# Single exponential decay function with fixed E_inf
# ----------------------------
def single_exponential_decay(t, A, tau, Einf):
    return Einf + A * np.exp(-t / tau)

# ----------------------------
# Load Adiabatic Energies
# ----------------------------
def load_adiabatic_energies(start_step, end_step, path_template, scale=au2ev, time_scale=au2fs):
    energies = []
    time = None
    for step in range(start_step, end_step):
        path = path_template.format(step=step)
        energy = sp.load_npz(path).todense().real
        if time is None:
            time = np.arange(energy.shape[0]) * time_scale
        energies.append(np.diag(energy))
    return np.array(energies) * scale, time

# ----------------------------
# Main Analysis Function
# ----------------------------
def calculate_total_energy_decay_to_S1(
    adiabatic_energies, methods, icond_range, method_config,
    smoothing_window=100, fit_time_cutoff=3000, verbose=True
):
    energy_results = {}
    average_energies = {}
    time_results = {}
    decay_fit_results = {}

    # --- Compute population-weighted average S1 energy ---
    all_pops = []
    all_E1 = []

    for method in methods:
        path_template = method_config[method]["path_template"]
        for icond in icond_range:
            with h5py.File(path_template.format(icond=icond), 'r') as f:
                sh_pop = np.array(f['sh_pop_adi/data'])
                n_steps = min(sh_pop.shape[0], adiabatic_energies.shape[0])
                pops = sh_pop[:n_steps, 1]  # S1 population
                E1 = adiabatic_energies[:n_steps, 1]  # S1 energy
                all_pops.append(pops)
                all_E1.append(E1)

    all_pops = np.array(all_pops)
    all_E1 = np.array(all_E1)
    E1_avg = np.sum(all_pops * all_E1) / np.sum(all_pops)

    for method in methods:
        print(f"\nProcessing method: {method}")
        try:
            method_results = []
            method_times = None
            path_template = method_config[method]["path_template"]

            initial_total_energies = []

            for icond in icond_range:
                path = path_template.format(icond=icond)
                with h5py.File(path, 'r') as f:
                    sh_pop = np.array(f['sh_pop_adi/data'])
                    time = np.array(f['time/data']) * au2fs

                    n_steps = min(sh_pop.shape[0], adiabatic_energies.shape[0])
                    sh_pop = sh_pop[:n_steps, :]
                    time = time[:n_steps]
                    Hvib_rolled = np.roll(adiabatic_energies[:n_steps, :], -icond, axis=0)

                    total_energy = np.sum(sh_pop * Hvib_rolled, axis=1)
                    method_results.append(total_energy)
                    initial_total_energies.append(total_energy[0])

                    if method_times is None:
                        method_times = time

            total_energy_avg = np.mean(method_results, axis=0)
            total_energy_std = np.std(method_results, axis=0)
            smoothed_total = uniform_filter1d(total_energy_avg, size=smoothing_window)

            E0_fixed = np.mean(initial_total_energies)
            Einf_fixed = E1_avg

            fit_mask = method_times < fit_time_cutoff
            fit_times = method_times[fit_mask]
            fit_values = smoothed_total[fit_mask]

            # --- Fit to single exponential decay ---
            def fit_func_singleexp(t, A, tau):
                return single_exponential_decay(t, A, tau, Einf_fixed)

            p0 = [E0_fixed - Einf_fixed, 500]
            bounds = (
                [-np.inf, 1e-2],
                [np.inf, 1e4]
            )

            popt, pcov = curve_fit(
                fit_func_singleexp,
                fit_times,
                fit_values,
                p0=p0,
                bounds=bounds
            )

            A, tau = popt
            perr = np.sqrt(np.diag(pcov))
            fit_curve = fit_func_singleexp(method_times, A, tau)

            if verbose:
                plt.figure()
                plt.plot(method_times, total_energy_avg, label='Raw Total Energy', alpha=0.4)
                plt.fill_between(
                    method_times,
                    total_energy_avg - total_energy_std,
                    total_energy_avg + total_energy_std,
                    color='gray', alpha=0.3, label='±1σ'
                )
                plt.plot(method_times, smoothed_total, label='Smoothed Total Energy', color='red')
                plt.plot(method_times, fit_curve, '--',
                         label=f'Single-exp Fit: A={A:.2f}, τ={tau:.1f} ± {perr[1]:.1f} fs',
                         color='black')
                plt.axhline(E0_fixed, color='green', linestyle=':', label=f'Initial ≈ {E0_fixed:.3f} eV')
                plt.axhline(Einf_fixed, color='blue', linestyle=':', label=f'S₁ Energy ≈ {Einf_fixed:.3f} eV (weighted)')

                plt.xlabel("Time (fs)", fontsize=20)
                plt.ylabel("Average Excitation Energy (eV)", fontsize=22)
             #   plt.title(f"{method} – Single Exponential Fit to Energy Relaxation", fontsize=22)
                plt.tick_params(axis='both', which='major', labelsize=22)
                #plt.legend(fontsize=14, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=True, ncol=1)
                plt.gcf().set_size_inches(12, 6)
                plt.tight_layout(rect=[0, 0, 0.75, 1])
                #plt.grid(True)
                plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))

                plt.show()

            print(f"Average total energy at t=0: {E0_fixed:.6f} eV")
            print(f"Fitted single exponential (t < {fit_time_cutoff} fs):")
            print(f"  A = {A:.2f}, τ = {tau:.2f} ± {perr[1]:.2f} fs")
            print(f"  Einf (fixed) = {Einf_fixed:.6f} eV")

            energy_results[method] = (method_results, [A, tau])
            average_energies[method] = smoothed_total
            time_results[method] = method_times
            decay_fit_results[method] = {
                "fit_curve": fit_curve,
                "params": (A, tau, Einf_fixed),
                "errors": (*perr, 0.0)
            }

        except Exception as e:
            print(f"Error in method {method}: {e}")

    return energy_results, average_energies, time_results, decay_fit_results


# In[6]:


# ----------------------------
# Energy file info
# ----------------------------
start_step = 1000
end_step = 3990
energy_path_template = "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/0500_300K/isomer2/step3/AS_30_39/Hvib_ci_{step}_re.npz"
get_ipython().run_line_magic('matplotlib', 'inline')

# ----------------------------
# Load adiabatic energies
# ----------------------------
adiabatic_energies, _ = load_adiabatic_energies(start_step, end_step, energy_path_template)

# ----------------------------
# Method configuration
# ----------------------------
method_config = {
    "DISH": {
        "path_template": "/projects/academic/alexeyak/miguelrecio/TiO2/second_research/0500_300K/isomer2/step4/DISH_icond_{icond}/mem_data.hdf"
    }
}

# ----------------------------
# Initial condition indices
# ----------------------------
icond_range = range(1, 3000, 100)

# ----------------------------
# Run energy decay analysis with fixed decay to S₁
# ----------------------------
methods = ["DISH"]
total_results, total_avg, total_times, total_fits = calculate_total_energy_decay_to_S1(
    adiabatic_energies=adiabatic_energies,
    methods=methods,
    icond_range=icond_range,
    method_config=method_config,
    smoothing_window=100,
    verbose=True
)


# In[ ]:




