import os, glob, time, h5py, warnings
import multiprocessing as mp
import matplotlib.pyplot as plt   # plots
import numpy as np
import scipy.sparse as sp
from scipy.optimize import curve_fit
from liblibra_core import *
import util.libutil as comn
import libra_py
from libra_py import units, data_conv #, dynamics_plotting
import libra_py.dynamics.tsh.compute as tsh_dynamics
#import libra_py.dynamics.tsh.plot as tsh_dynamics_plot
#import libra_py.data_savers as data_savers
import libra_py.workflows.nbra.decoherence_times as decoherence_times
import libra_py.data_visualize
from recipes import dish_rev2023_nbra, fssh_nbra, msdm_nbra
#from matplotlib.mlab import griddata
#%matplotlib inline 
warnings.filterwarnings('ignore')
path_to_save_sd_Hvibs = '../step3/res-sd-30-occ'
istep = 1000    # the first timestep to read
fstep = 3990 # the last timestep to read
nsteps = fstep - istep
NSTEPS = nsteps
print(F"Number of steps = {nsteps}")
# ================= Reading the data
#================== Read energies =====================
E = []
for step in range(istep,fstep):
    energy_filename = F"{path_to_save_sd_Hvibs}/Hvib_ci_{step}_re.npz"
    energy_mat = sp.load_npz(energy_filename)
    # For data conversion we need to turn np.ndarray to np.array so that 
    # we can use data_conv.nparray2CMATRIX
    E.append( np.array( np.diag( energy_mat.todense() ) ) )
E = np.array(E)
NSTATES = E[0].shape[0]
#================== Read time-overlap =====================
St = []
for step in range(istep,fstep):        
    St_filename = F"{path_to_save_sd_Hvibs}/St_ci_{step}_re.npz"
    St_mat = sp.load_npz(St_filename)
    St.append( np.array( St_mat.todense() ) )
St = np.array(St)
#================ Compute NACs and vibronic Hamiltonians along the trajectory ============    
NAC = []
Hvib = [] 
for c, step in enumerate(range(istep,fstep)):
    nac_filename = F"{path_to_save_sd_Hvibs}/Hvib_ci_{step}_im.npz"
    nac_mat = sp.load_npz(nac_filename)
    NAC.append( np.array( nac_mat.todense() ) )
    Hvib.append( np.diag(E[c, :])*(1.0+1j*0.0)  - (0.0+1j)*nac_mat[:, :] )
NAC = np.array(NAC)
Hvib = np.array(Hvib)
class abstr_class:
    pass
def compute_model(q, params, full_id):
    timestep = params["timestep"]
    nst = params["nstates"]
    obj = abstr_class()
    obj.ham_adi = data_conv.nparray2CMATRIX( np.diag(E[timestep, : ]) )
    obj.nac_adi = data_conv.nparray2CMATRIX( NAC[timestep, :, :] )
    obj.hvib_adi = data_conv.nparray2CMATRIX( Hvib[timestep, :, :] )
    obj.basis_transform = CMATRIX(nst,nst); obj.basis_transform.identity()  #basis_transform
    obj.time_overlap_adi = data_conv.nparray2CMATRIX( St[timestep, :, :] )
    
    return obj
#print('Number of steps:', NSTEPS)
#print('Number of states:', NSTATES)
# ================= Computing the energy gaps and decoherence times
HAM_RE = []
for step in range(E.shape[0]):
    HAM_RE.append( data_conv.nparray2CMATRIX( np.diag(E[step, : ]) ) )
# Average decoherence times and rates
tau, rates = decoherence_times.decoherence_times_ave([HAM_RE], [0], NSTEPS, 0)
# Computes the energy gaps between all states for all steps
dE = decoherence_times.energy_gaps_ave([HAM_RE], [0], NSTEPS)
# Decoherence times in fs
deco_times = data_conv.MATRIX2nparray(tau) * units.au2fs
# Zero all the diagonal elements of the decoherence matrix
np.fill_diagonal(deco_times, 0)
# Saving the average decoherence times
np.savetxt('decoherence_times.txt',deco_times.real)
# Computing the average energy gaps
gaps = MATRIX(NSTATES, NSTATES)
for step in range(NSTEPS):
    gaps += dE[step]
gaps /= NSTEPS
rates.show_matrix("decoherence_rates.txt")
gaps.show_matrix("average_gaps.txt")
#sys.exit(0)
#================== Model parameters ====================
model_params = { "timestep":0, "icond":0,  "model0":0, "nstates":NSTATES }
#=============== Some automatic variables, related to the settings above ===================
#=============== Some automatic variables, related to the settings above ===================
#############
NSTEPS = 3000
dyn_general = { "nsteps":NSTEPS, "ntraj":200, "nstates":NSTATES, "dt":1.0*units.fs2au,                                                 
                "decoherence_rates":rates, "ave_gaps":gaps,                
                "progress_frequency":0.1, "which_adi_states":range(NSTATES), "which_dia_states":range(NSTATES),
                "mem_output_level":2,
                "properties_to_save":[ "timestep", "time","se_pop_adi", "sh_pop_adi" ],
                "prefix":F"NBRA", "prefix2":F"NBRA", "isNBRA":0, "nfiles": nsteps-1 
              }
#dyn_general.update({"ham_update_method":2}) 
#
#dyn_general.update( {"ham_transform_method":0 }) 
#
#dyn_general.update( {"time_overlap_method":0 }) 
#
#dyn_general.update({"nac_update_method":0 }) 
#
#dyn_general.update( {"hvib_update_method":0 }) 
#
#dyn_general.update( {"force_method":0, "rep_force":1} ) 
#
#dyn_general.update({"hop_acceptance_algo":32, "momenta_rescaling_algo":0 }) 
##########################################################
#============== Select the method =====================
#dish_nbra.load(dyn_general); prf = "DISH"  # DISH
#fssh_nbra.load(dyn_general); prf = "FSSH"  # FSSH
#fssh2_nbra.load(dyn_general); prf = "FSSH2"  # FSSH2
#gfsh_nbra.load(dyn_general); prf = "GFSH"  # GFSH
#ida_nbra.load(dyn_general); prf = "IDA"  # IDA
#mash_nbra.load(dyn_general); prf = "MASH"  # MASH
msdm_nbra.load(dyn_general); prf = "MSDM"  # MSDM
##########################################################
#=================== Initial conditions =======================
#============== Nuclear DOF: these parameters don't matter much in the NBRA calculations ===============
nucl_params = {"ndof":1, "init_type":3, "q":[-10.0], "p":[0.0], "mass":[2000.0], "force_constant":[0.01], "verbosity":-1 }
#============== Electronic DOF: Amplitudes are sampled ========
elec_params = {"ndia":NSTATES, "nadi":NSTATES, "verbosity":-1, "init_dm_type":0}
###########
istate = 10 
###########
elec_params.update( {"init_type":1,  "rep":1,  "istate":istate } )  # how to initialize: random phase, adiabatic representation
if prf=="MASH":
    istates = list(np.zeros(NSTATES))
    istates[istate] = 1.0
    elec_params.update( {"init_type":4,  "rep":1,  "istate":3, "istates":istates } )  # different initialization for MASH
def function1(icond):
#   time.sleep(rnd.uniform(0.0, 1.0) * 20 )
    time.sleep(icond * 0.01 )
    rnd = Random()
    mdl = dict(model_params)
    mdl.update({"icond": icond})  #create separate cop
    dyn_gen = dict(dyn_general)
    dyn_gen.update({"prefix":F"{prf}_icond_{icond}", "prefix2":F"{prf}_icond_{icond}" })
    res = tsh_dynamics.generic_recipe(dyn_gen, compute_model, mdl, elec_params, nucl_params, rnd)
################################
nthreads = 4
ICONDS = list(range(1,3000,100))
################################
pool = mp.Pool(nthreads)
pool.map(function1, ICONDS)
pool.close()                            
pool.join()

