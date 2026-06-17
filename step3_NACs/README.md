# Project_TiO2_fluxionality

**scripts and selected outputs from computing the nonadiabatic couplings (NACs) between pairs of TD-DFT states**

**scripts**:

The inputs for this step are obtained from the step 2 TD-DFT calculations. Therefore, there is no separate `inputs` folder for step 3. The required files are loaded directly by the `step3_NACs.py` script from the corresponding step 2 output directories.

The provided `step3_NACs.py` script is written as a case example for isomer 1 at 100 K. Before running the calculation for a different temperature, system, or trajectory segment, the paths to the CP2K `.log` files and to the corresponding `res` folders generated in step 2 must be modified accordingly inside the script.

The `submit_template.slm` file is provided to submit the NAC calculation to a HPC system.

To prepare this step, make sure that the previous step has been properly completed, check that the paths to the step 2 output files are correctly defined in `step3_NACs.py`, and run:

**sbatch `submit_template.slm`**

The `representation_PDs.py` script can be used, after completing the step 3 calculations, to generate probability distribution (PD) plots for the NACs and energy gaps. The script should be modified depending on the desired property, temperature, and input files to be represented.

The `representation_PDs.py` script was used to generate, for both 100 K and 300 K, the probability distributions of the NACs between the ground state and the first excited states, NACs between adjacent excited states, the S0-S1 energy gap, and the energy gaps between adjacent excited states.

**outputs**:

The complete step 3 NAC output files are large. Therefore, as in step 2, only the results corresponding to the last 500 geometries of the AIMD trajectories are provided. These correspond to geometries 2500-3000.

The `key_outputs` directory contains two subdirectories, one for each temperature: `100K` and `300K`. Each subdirectory contains the compressed `.bz2` files with the NAC results for isomer 1 at the corresponding temperature.

For the 100 K case, the NAC results are provided as compressed files named:

`isomer1_100K_step3_part*.bz2`

For the 300 K case, the NAC results are provided as compressed files named:

`isomer1_300K_step3_part*.bz2`

To unpack the compressed output files, decompress each `.bz2` file separately. For example:

**tar -xvjf isomer1_100K_step3_part1.bz2**

**tar -xvjf isomer1_300K_step3_part1.bz2**

Repeat the same command for the remaining `step3_part*.bz2` files.

In addition to the compressed NAC files, the `key_outputs` folder also contains representative figures generated from the step 3 results. These include the probability distributions of the S0-S1 energy gap, the adjacent-state NACs, and selected S0-S1 NAC distributions at both temperatures.

-------------

Take into account that the provided `step3_NACs.py` script is written for the 100 K case of isomer 1. The paths to the step 2 `.log` files and `res` folders must be adapted depending on the system, temperature, and trajectory range under study. Similarly, the `representation_PDs.py` script must be edited to generate the desired probability distribution plots from the selected NAC or energy-gap data.