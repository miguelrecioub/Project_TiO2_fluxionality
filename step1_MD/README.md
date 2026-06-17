# Project_TiO2_fluxionality

**inputs and scripts should be placed in the same directory to start the calculation**

**inputs**:

AIMD simulations are performed with CP2K. This folder contains the input files required to run the molecular dynamics calculations. As a case example, the provided `md.inp` file is written for a 100 K AIMD calculation of isomer 1. The corresponding starting geometries for the three isomers are also included as `.xyz` files.

**scripts**:

The `submit.slm` file is provided to submit the step 1 calculation to a HPC system. The `RMSD.py` script can be used to compute and analyze the structural deviations along the AIMD trajectories.

**outputs**:

This folder contains the relevant AIMD results used in the project. These include the AIMD trajectories for isomers 1 and 2 at both temperatures, the evolution of the Ti_flux–O_flux distance, the evolution of representative 4NN Ti–O distances, and the RMSD analysis figure.

-------------

To prepare this step, combine the required input and script files in the same directory and run:

**sbatch `submit.slm`**