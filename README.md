# Project\_titania\_nanocluster\_fluxionality



This project methodology is composed of four sequential steps:

(i) generation of the nuclear guiding trajectory through ab initio molecular dynamics (AIMD);

(ii) computation of the TD-DFT excitations, yielding excitation energies, molecular orbital (MO) overlaps, and MO time-overlaps;

(iii) computation of the nonadiabatic couplings (NACs); and

(iv) nonadiabatic molecular dynamics (NA-MD) simulations.



Accordingly, the repository is organized into four main directories, each containing the inputs and scripts required to run the corresponding step. Additional scripts for visualizing the results discussed in the manuscript and Supporting Information are also provided. These scripts are intended as templates and should be adapted as needed for the specific system under study.



We recommend consulting [this tutorial](https://github.com/compchem-cybertraining/Tutorials_Libra/tree/master/6_dynamics/2_nbra_workflows) for detailed guidance on how to successfully perform each of the methodological steps.



Along with the inputs and scripts for the four sequential steps, each methodological-step directory also contains a folder with representative results for the fluxional and rigid isomers. These results include optimized structures of the three isomers; AIMD trajectories for isomers 1 and 2 at both temperatures; the evolution of the fluxional-pair interatomic distance; excited-state energy evolutions; nonadiabatic couplings; excited-state probability distribution functions; and relaxation and recombination dynamics. Results reported in the Supporting Information, such as the minimum-energy path between isomers 1 and 3 and the dynamical-mode analysis, are also included.



Due to the large size of the complete output files, the results' files provided for steps 2, 3, and 4 are limited to the last 500 geometries of the corresponding trajectories. However, the input files and README files included within each step directory describe the workflow used for the full production calculations, which were performed over longer trajectories containing up to 3000 geometries.

