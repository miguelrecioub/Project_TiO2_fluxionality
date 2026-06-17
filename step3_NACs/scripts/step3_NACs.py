import os, glob, time, h5py, warnings
import numpy as np
import scipy.sparse as sp
from libra_py import units, data_stat, influence_spectrum, data_conv
import matplotlib.pyplot as plt
from liblibra_core import *
from libra_py.workflows.nbra import step3
import libra_py.packages.cp2k.methods as CP2K_methods
import multiprocessing as mp

import util.libutil as comn
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.workflows.nbra.decoherence_times as decoherence_times

params_active_space = {
    'lowest_orbital': 60-45, 'highest_orbital': 61+45, 'num_occ_orbitals': 30, 'num_unocc_orbitals': 30,
    'path_to_npz_files': os.getcwd()+'/../res', 'logfile_directory': os.getcwd()+'/../all_logfiles',
    'path_to_save_npz_files': os.getcwd()+'/../new_res'
}
new_lowest_orbital, new_highest_orbital = step3.limit_active_space(params_active_space)

params_mb_sd = {
          'lowest_orbital': new_lowest_orbital, 'highest_orbital': new_highest_orbital, 'num_occ_states': 30, 'num_unocc_states': 30,
          'isUKS': 0, 'number_of_states': 20, 'tolerance': 0.0, 'verbosity': 0, 'use_multiprocessing': True, 'nprocs': 4,
          'is_many_body': True, 'time_step': 1.0, 'es_software': 'cp2k',
          'path_to_npz_files': params_active_space['path_to_save_npz_files'],
          'logfile_directory': os.getcwd()+'/../all_logfiles',
          'path_to_save_sd_Hvibs': os.getcwd()+'/res-sd-30-occ',
          'outdir': os.getcwd()+'/res-sd-30-occ', 'start_time': 1000, 'finish_time': 3999, 'sorting_type': 'energy',
         }

step3.run_step3_sd_nacs_libint(params_mb_sd)
clear_output()
