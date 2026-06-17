# Project_TiO2_fluxionality

**scripts and selected outputs from performing the nonadiabatic molecular dynamics (NA-MD) runs**

**scripts**:

The inputs for this step are obtained from the previous steps. Therefore, there is no separate `inputs` folder for step 4. The required files are loaded directly by the NA-MD scripts from the corresponding step 3 output directories.

The `NAMD_FSSH.py` and `NAMD_MSDM.py` scripts are used to run the NA-MD simulations with different trajectory surface hopping (TSH) schemes. In this project, we have used FSSH and mSDM approaches to simulate the relaxation and recombination dynamics.

Before running these scripts, the paths to the required input files from step 3 must be modified according to the system, temperature, trajectory range, and type of dynamics to be simulated.

The `recipes` folder contains the predefined Libra recipes for the different TSH schemes. This folder should not be modified. It should simply be kept in the same directory where the NA-MD calculations are conducted.

The `submit_FSSH.slm` and `submit_MSDM.slm` files are provided to submit the corresponding NA-MD calculations to a HPC system.

To prepare this step, make sure that the previous steps have been properly completed, check that the paths to the step 3 output files are correctly defined in the corresponding `NAMD_*.py` script, and run:

**sbatch `submit_FSSH.slm`**

or

**sbatch `submit_MSDM.slm`**

depending on the selected method.

After finishing the NA-MD runs, the `newGS_populations.py` script can be used to plot the ground-state population evolution associated with electron-hole recombination.

The `average_excitation_energy_decay_3functions.py` script can be used to analyze and fit the average excitation-energy decay associated with relaxation from a higher-energy bright state to S1. This analysis requires the corresponding NA-MD output files and should be adapted to the location of the selected input data.

**outputs**:

The `outputs` directory contains the selected NA-MD results used in this project. The results are organized into two subdirectories, one for each temperature: `100K` and `300K`.

Each temperature folder contains two types of compressed `.bz2` files:

1. Recombination dynamics, corresponding to simulations starting from S1 and following the population transfer to the ground state S0.

2. Relaxation dynamics, corresponding to simulations starting from a higher-energy bright excited state and following the relaxation toward S1.

For the 100 K case, the provided compressed files are:

`isomer1_100K_recomb_FSSH.bz2`

`isomer1_100K_recomb_mSDM.bz2`

`isomer1_100K_relax_FSSH.bz2`

`isomer1_100K_relax_mSDM.bz2`

For the 300 K case, the provided compressed files are:

`isomer1_300K_recomb_FSSH.bz2`

`isomer1_300K_recomb_MSDM.bz2`

`isomer1_300K_relax_FSSH.bz2`

`isomer1_300K_relax_MSDM.bz2`

To unpack the compressed output files, decompress each `.bz2` file separately. For example:

**tar -xvjf isomer1_100K_recomb_FSSH.bz2**

**tar -xvjf isomer1_100K_relax_FSSH.bz2**

**tar -xvjf isomer1_300K_recomb_FSSH.bz2**

**tar -xvjf isomer1_300K_relax_FSSH.bz2**

Repeat the same command for the remaining recombination and relaxation files.

In addition to the compressed NA-MD output files, the `outputs` directory also contains representative figures generated from the step 4 results, including the relaxation and recombination dynamics reported in the manuscript and Supporting Information.

-------------

Take into account that the provided NA-MD scripts are templates. The paths to the input files generated in step 3 must be adapted depending on the system, temperature, trajectory range, selected initial state, and type of process under study. In particular, recombination calculations start from S1, whereas relaxation calculations start from a higher-energy bright excited state and follow the population decay toward S1.