import f90nml
from os import path
from types import SimpleNamespace
from iop.logger import log, error

def readNamelist(pathNamelist):
    """
    Reads a Fortran namelist file and returns its contents as a dictionary.

    Parameters:
        namelist_path (str): Path to the namelist file.

    Returns:
        dict: Dictionary containing the namelist parameters.
    """
    try:
        nml = f90nml.read(pathNamelist)
        return nml
    
    except Exception as e:
        print(f"Error reading namelist: {e}")
        return None
    

class ParamConfig:
    """Class to hold SCHISM parameter configuration with easy access to sections."""
    def __init__(self, param_dict):
        # Create sections as attributes
        for section, values in param_dict.items():
            # Convert section name to lowercase for cleaner access
            section_name = section.lower()
            setattr(self, section_name, SimpleNamespace(**values))
    
    def __repr__(self):
        sections = [attr for attr in dir(self) if not attr.startswith('_')]
        return f"ParamConfig(sections={sections})"

def readParam(pathNamelist):
    """
    Specialized function for reading param.nml. Verifies all mandatory CORE parameters are present
    and sets default values for optional parameters.
    
    Parameters:
        pathNamelist (str): Path to the directory containing param.nml
        
    Returns:
        ParamConfig: Object containing the namelist parameters with attribute access
        
    Raises:
        ValueError: If any mandatory parameters are missing
    """
    # Mandatory parameters in CORE section
    mandatory_params = {
        'ipre',      # Pre-processor flag
        'ibc',       # Baroclinic option
        'ibtp',      # Barotropic option
        'rnday',     # Total run time in days
        'dt',        # Time step in seconds
        'msc2',      # Grid parameter for WWM
        'mdc2',      # Grid parameter for WWM
        'ntracer_gen', # Number of generic tracers
        'ntracer_age', # Number of age tracers
        'sed_class',   # Number of sediment classes
        'eco_class',   # Number of ecosystem classes
        'nspool',      # Output step spool
        'ihfskip',     # Stack spool
        'nbins_veg_vert' # Number of vertical bins in vegetation
    }

    # Complete default values for OPT section
    opt_defaults = {
        # Basic model configuration
        'gen_wsett': 0.0,
        'flag_fib': 1,
        'ics': 1,
        'rearth_pole': 6378206.4,
        'rearth_eq': 6378206.4,
        'indvel': 0,
        'imm': 0,
        'ibdef': 10,
        'ihot': 0,
        'ihydraulics': 0,
        'izonal5': 0,
        'slam0': -124.0,
        'sfea0': 45.0,

        # Numerical parameters
        'iupwind_mom': 0,
        'ihorcon': 0,
        'hvis_coef0': 0.025,
        'ishapiro': 1,
        'shapiro0': 0.5,
        'niter_shap': 1,
        'ihdif': 0,
        'thetai': 0.6,
        'drampbc': 0.0,
        'dramp': 1.0,
        
        # Advection parameters
        'nadv': 1,
        'dtb_min': 10.0,
        'dtb_max': 30.0,
        'h0': 0.01,
        'nchi': 0,
        'dzb_min': 0.5,
        'hmin_man': 1.0,
        
        # Physical parameters
        'ncor': 0,
        'rlatitude': 46.0,
        'coricoef': 0.0,
        'nws': 0,
        'wtiminc': None,  # Will be set to dt
        'iwind_form': 1,
        'iwindoff': 0,
        'drampwind': 1.0,
        
        # Conservation parameters
        'ihconsv': 0,
        'isconsv': 0,
        'i_hmin_airsea_ex': 2,
        'i_hmin_salt_ex': 2,
        
        # Turbulence parameters
        'itur': 0,
        'dfv0': 0.01,
        'dfh0': 1.0e-4,
        'h1_pp': 20.0,
        'h2_pp': 50.0,
        'vdmax_pp1': 0.01,
        'vdmax_pp2': 0.01,
        'vdmin_pp1': 1.0e-5,
        'vdmin_pp2': 1.0e-5,
        'tdmin_pp1': 1.0e-5,
        'tdmin_pp2': 1.0e-5,
        'mid': 'KL',
        'stab': 'KC',
        'xlsc0': 0.1,

        # Time and initialization
        'ibcc_mean': 0,
        'flag_ic': [1] * 11,  # Array of 11 ones
        'start_year': 2000,
        'start_month': 1,
        'start_day': 1,
        'start_hour': 0.0,
        'utc_start': 0.0,

        # Transport parameters
        'itr_met': 3,
        'h_tvd': 5.0,
        'eps1_tvd_imp': 1.0e-4,
        'eps2_tvd_imp': 1.0e-14,
        'ip_weno': 2,
        'courant_weno': 0.5,
        'ntd_weno': 1,
        'nquad': 2,
        'epsilon1': 1.0e-15,
        'i_epsilon2': 1,
        'epsilon2': 1.0e-10,
        'epsilon3': 1.0e-25,
        'ielad_weno': 0,
        'small_elad': 1.0e-4,
        'i_prtnftl_weno': 0,

        # Nudging parameters
        'inu_tr': [0] * 11,  # Array of 11 zeros
        'step_nu_tr': 86400.0,
        'vnh1': 400.0,
        'vnh2': 500.0,
        'vnf1': 0.0,
        'vnf2': 0.0,

        # Solver parameters
        'moitn0': 50,
        'mxitn0': 1500,
        'rtol0': 1.0e-12,
        'iflux': 0,
        'inter_mom': 0,
        'h_bcc1': 100.0,
        'inu_elev': 0,
        'inu_uv': 0,

        # Additional parameters
        'ihhat': 1,
        'kr_co': 1,
        'rmaxvel': 5.0,
        'velmin_btrack': 1.0e-4,
        'btrack_nudge': 9.013e-3,
        'ibtrack_test': 0,
        'irouse_test': 0,
        'inunfl': 0,
        'shorewafo': 0,
        'ic_elev': 0,
        'nramp_elev': 0,
        'inv_atm_bnd': 0,
        'prmsl_ref': 101325.0,
        's1_mxnbt': 0.5,
        's2_mxnbt': 3.5,

        # Wave and coupling parameters
        'iharind': 0,
        'icou_elfe_wwm': 0,
        'drampwafo': 0.0,
        'nstep_wwm': 1,
        'hmin_radstress': 1.0,
        'turbinj': 0.15,
        'turbinjds': None,  # Not specified in defaults
        'alphaw': 1.0,

        # Wave-related flags
        'fwvor_advxy_stokes': 1,
        'fwvor_advz_stokes': 1,
        'fwvor_gradpress': 1,
        'fwvor_breaking': 1,
        'fwvor_streaming': 1,
        'fwvor_wveg': 0,
        'fwvor_wveg_NL': 0,
        'wafo_obcramp': 0,

        # Physical properties
        'iwbl': 0,
        'cur_wwm': 0,
        'if_source': 0,
        'dramp_ss': 2.0,
        'ieos_type': 0,
        'ieos_pres': 0,
        'eos_a': -0.1,
        'eos_b': 1001.0,
        'slr_rate': 120.0,
        'rho0': 1000.0,
        'shw': 4184.0,
        'iveg': 0,
        'nstep_ice': 1,
        'h1_bcc': 50.0,
        'h2_bcc': 100.0,
        'hw_depth': 1.0e6,
        'hw_ratio': 0.5,
        'iunder_deep': 0,
        'level_age': -999,
        'vclose_surf_frac': 1.0,
        'iadjust_mass_consv0': 0,
        'ipre2': 0,
        'ielm_transport': 0,
        'max_subcyc': 10,
        'hmin_airsea_ex': 0.2,
        'hmin_salt_ex': 0.2,
        'itransport_only': 0,
        'iloadtide': 0,
        'loadtide_coef': 0.1,
        'nu_sum_mult': 1,
        'h_massconsv': 2.0,
        'rinflation_icm': 1.0e-3,
        'lev_tr_source': -9,
        'iprecip_off_bnd': 0,
        'model_type_pahm': 10,
        'stemp_stc': 0,
        'stemp_dz': 1.0,

        # Vegetation parameters
        'veg_lai': 1.0,  # Leaf Area Index [-]; used if iveg=2. Range: 0-10
        'veg_cw': 1.5,   # Calibration coefficient in diameter of bent leaf [-]; used if iveg=2. Range: 0-25
    }
    
    # Read the namelist file
    param_dict = readNamelist(path.join(pathNamelist, 'param.nml'))
    if param_dict is None:
        raise ValueError("Failed to read param.nml")
        
    # Verify all mandatory parameters are present in CORE section
    if 'CORE' not in param_dict:
        raise ValueError("CORE section missing from param.nml")
        
    missing_params = mandatory_params - set(param_dict['CORE'].keys())
    if missing_params:
        raise ValueError(f"Missing mandatory parameters in CORE section: {missing_params}")
    
    # Ensure OPT section exists
    if 'OPT' not in param_dict:
        param_dict['OPT'] = {}
    
    # Set default values for optional parameters if not present
    for key, default_value in opt_defaults.items():
        if key not in param_dict['OPT']:
            param_dict['OPT'][key] = default_value


    # Check core parameters
    core = param_dict['CORE']
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

    # Convert to ParamConfig object
    return ParamConfig(param_dict)
    

