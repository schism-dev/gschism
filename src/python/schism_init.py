import sys
import os
import python.cuda_functions as cufunc

current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.dirname(current_dir)  # go up two levels to reach src/
sys.path.append(src_dir)

from os import path
import cupy as cp
from iop.logger import *
from iop.readNamelist import readNamelist, readParam
from python.schism_consts import *
from python.schism_consts import *


def readNamelist(in_dir):
    # Reads the namelist file params.nml
    nml = readParam(in_dir)

    return nml

def schism_init(core, opt):



    # Count the number of tracer models and tracers
    # Each trace module has a pre-defined ID as follows:
    #  1: T
    #  2: S
    #  3: GEN
    #  4: AGE
    #  5: SED3D
    #  6: EcoSim
    #  7: ICM
    #  8: CoSINE
    #  9: Feco
    # 10: TIMOR
    # 11: FABM
    # 12: DVD numerical mixing analysis of Klingbeit

    # Init. # of tracers in models 1:natrm
    # A tracer model is activated iff ntrs()>0
    ntrs = cp.zeros(natrm, dtype=int)
    ntrs[0:2] = cp.array([1, 1])  # T,S counted
    tr_mname = ['TEM', 'SAL']

    # Total number of tracers (including T,S)
    # The big tracer arrays are: tr_el[ntracers,nvrt,nea2],tr_nd0[ntracers,nvrt,npa]
    ntracers = ntrs.sum().item() # Including T,S
    if (mntracers < ntracers):
        error('INIT: mntracers<ntracers')


    # XXXXX Initialize time history variables
    #  uthnd=-99.d0; vthnd=-99.d0; eta_mean=-99.d0; !uth=-99.d0; vth=-99.d0; !flags


    # Compute total number of time steps 
    ntime = cp.round(core.rnday*86400.0/core.dt)
    nrec = min(ntime, core.ihfskip)/core.nspool


    # Calculate the index ranges for each tracer module
    irange_tr = cp.zeros((2, natrm), dtype=int)
    irange_tr[0, 1:] = cp.cumsum(ntrs[:-1])    
    irange_tr[1,  :] = cp.cumsum(ntrs[:])
    
    log(f'# of tracers in each module: {ntrs}')
    log(f'Total # of tracers: {ntracers}')
    log(f'Index ranges of each module: {irange_tr}')

    # True b-tropic model cannot have tracer models (to ensure T,S always come first in arrays)
    if(core.ibc == 1 and core.ibtp == 0 and ntracers > 2):
        error('INIT: true barotropic model cannot have tracer models')


    if(opt.inunfl==0):

    
        cufunc.ul0UpdateZCoords()
        
        cufunc.ul0ComputeDryElements()

        cufunc.ul0ComputeDryNodes()
        
        cufunc.ul0UpdateRewettedNodes(0, 0)
        
        cufunc.ul0UpdateSides(0, 0)

        cufunc.ul0UpdateSides2()
        

        

    else:
        print("inunfl = 1 not implemented")
        sys.exit()
        # levels1(0, 0)

    return int(ntime.item()), int(nrec)


    

