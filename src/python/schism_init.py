import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.dirname(current_dir)  # go up two levels to reach src/
sys.path.append(src_dir)

import f90nml  # XXX check if we need a more standard package
from os import path
import cupy as cp
from iop.logger import *
from python.schism_consts import *

in_dir = '/home/jo.gonzalez/src/gschism_old/test/schism_verification_tests/Test_CORIE'


def schism_init(in_dir):

    # Reads the namelist file params.nml
    nml = f90nml.read(path.join(in_dir, 'param.nml'))
    
    
    core = nml['core']

    # Check core parameters

    # Check model type parameters
    if core['ibc'] not in {0, 1}:
        error('Unknown ibcc')
    if core['ibtp'] not in {0, 1}:
        error('Unknown ibtp')

    # Print model configuration
    if core['ibc'] == 0:
        log('You are using baroclinic model')
    else:
        message = ('Barotropic model without ST calculation' 
                   if core['ibtp'] == 0 
                   else 'Barotropic model with ST calculation')
        log(message)

    # Validate time parameters
    if core['dt'] <= 0 or core['rnday'] <= 0:
        error('Illegal dt,rnday')

    # Validate pool settings
    if core['nspool'] == 0 or core['ihfskip'] == 0:
        error('Zero nspool')
    if core['ihfskip'] % core['nspool'] != 0:
        error('ihfskip/nspool /= integer')

    # Validate vegetation parameters
    if core['nbins_veg_vert'] <= 0:
        error('INIT: nbins_veg_vert<=0')


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
    ntrs[0:2] = 1  # T,S counted as separate model
    tr_mname = ['TEM', 'SAL']

    # Total number of tracers (including T,S)
    # The big tracer arrays are: tr_el[ntracers,nvrt,nea2],tr_nd0[ntracers,nvrt,npa]
    ntracers = ntrs.sum() # Including T,S
    if (mntracers < ntracers):
        error('INIT: mntracers<ntracers')


    # Calculate the index ranges for each tracer module
    irange_tr = cp.zeros((2, natrm), dtype=int)
    irange_tr[0, 1:] = cp.cumsum(ntrs[:-1])    
    irange_tr[1,  :] = cp.cumsum(ntrs[:])
    
    log(f'# of tracers in each module: {ntrs}')
    log(f'Total # of tracers: {ntracers}')
    log(f'Index ranges of each module: {irange_tr}')

    # True b-tropic model cannot have tracer models (to ensure T,S always come first in arrays)
    if(core['ibc'] == 1 and core['ibtp'] == 0 and ntracers > 2):
        error('INIT: true barotropic model cannot have tracer models')





schism_init(in_dir)