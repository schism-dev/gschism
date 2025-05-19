from python.schism_consts import *

def processConstants(namelist, tidalData, ivcor: int, nvrt: int, kz: int, nsiglev: int,
                     nope: int, mnond: int,
                     templatePath="src/kernels/constants_template.h", 
                     outputPath="src/kernels/constants.h"):
    """
    Process constants_template.h using values from namelist, tidal data and grid parameters
    
    Parameters:
    -----------
    namelist : ParamConfig
        Configuration object from readNamelist
    tidalData : TidalData
        Object containing tidal forcing information
    ivcor : int
        Vertical coordinate type
    nvrt : int
        Number of vertical levels
    kz : int
        Number of z-levels
    nsiglev : int
        Number of sigma levels
    templatePath : str
        Path to constants_template.h
    outputPath : str
        Path where constants.h will be generated
    """
   
    # Define mapping of template placeholders to values
    # Format: 'PLACEHOLDER': (value, source_description)
    constantsMapping = {
        # From namelist CORE section
        '__BLOCK_SIZE__': (blockSize, 'cuda'),
        '__N_TRACERS__': (namelist.core.ntracer_gen, 'namelist.core'),
        '__IBC__': (namelist.core.ibc, 'namelist.core'),
        '__DT__': (namelist.core.dt, 'namelist.core'),
        
        # From namelist OPT section
        '__ICS__': (namelist.opt.ics, 'namelist.opt'),
        '__IMM__': (namelist.opt.imm, 'namelist.opt'),
        '__IVEG__': (namelist.opt.iveg, 'namelist.opt'),
        '__IHHAT__': (namelist.opt.ihhat, 'namelist.opt'),
        '__IHYDRAULICS__': (namelist.opt.ihydraulics, 'namelist.opt'),
        '__NCHI__': (namelist.opt.nchi, 'namelist.opt'),
        '__NRAMP_ELEV__': (namelist.opt.nramp_elev, 'namelist.opt'),
        '__ILOADTIDE__': (namelist.opt.iloadtide, 'namelist.opt'),
        '__LOADTIDE_COEF__': (namelist.opt.loadtide_coef, 'namelist.opt'),
        '__INV_ATM_BND__': (namelist.opt.inv_atm_bnd, 'namelist.opt'),
        '__PRMSL_REF__': (namelist.opt.prmsl_ref, 'namelist.opt'),
        '__VEG_LAI__': (namelist.opt.veg_lai, 'namelist.opt'),
        '__VEG_CW__': (namelist.opt.veg_cw, 'namelist.opt'),
        '__IUNDER_DEEP__': (namelist.opt.iunder_deep, 'namelist.opt'),
        '__H1_BCC__': (namelist.opt.h1_bcc, 'namelist.opt'),
        '__H2_BCC__': (namelist.opt.h2_bcc, 'namelist.opt'),
        '__HW_DEPTH__': (namelist.opt.hw_depth, 'namelist.opt'),
        '__HW_RATIO__': (namelist.opt.hw_ratio, 'namelist.opt'),
        '__REARTH_EQ__': (namelist.opt.rearth_eq, 'namelist.opt'),
        '__REARTH_POLE__': (namelist.opt.rearth_pole, 'namelist.opt'),
        '__THETAI__': (namelist.opt.thetai, 'namelist.opt'),
        '__RHO0__': (namelist.opt.rho0, 'namelist.opt'),
        '__DTB_MIN__': (namelist.opt.dtb_min, 'namelist.opt'),
        '__DTB_MAX__': (namelist.opt.dtb_max, 'namelist.opt'),
        '__VELMIN_BTRACK__': (namelist.opt.velmin_btrack, 'namelist.opt'),
        '__DZB_MIN__': (namelist.opt.dzb_min, 'namelist.opt'),
        '__H0__': (namelist.opt.h0, 'namelist.opt'),
        
        # From grid parameters
        '__IVCOR__': (ivcor, 'grid'),
        '__NVRT__': (nvrt, 'grid'),
        '__NZLEV__': (kz, 'grid'),  # kz is number of z-levels
        '__NSIGLEV__': (nsiglev, 'grid'),
        '__NOPEN_BOUNDARY_SEGMENTS__': (nope, 'grid'),
        '__MAX_NODES_PER_OPEN_BOUNDARY_SEGMENT__': (mnond, 'grid'),
        
        # From tidal data
        '__NTIP__': (tidalData.ntip, 'tidal'),
        '__TIP_DP__': (tidalData.tip_dp, 'tidal'),
        '__NBFR__': (tidalData.nbfr, 'tidal'),
        '__NETTYPE__': (tidalData.nettype, 'tidal'),
        '__NFLTYPE__': (tidalData.nfltype, 'tidal'),
        '__NETTYPE2__': (tidalData.nettype2, 'tidal'),
        '__NFLTYPE2__': (tidalData.nfltype2, 'tidal'),
        '__NNODE_ET__': (tidalData.nnode_et, 'tidal'),
        '__NNODE_FL__': (tidalData.nnode_fl, 'tidal')
    }
    
    # Verify all required parameters exist and have valid values
    missing_params = []
    for placeholder, (value, source) in constantsMapping.items():
        if value is None:
            missing_params.append(f"{placeholder} from {source}")
    
    if missing_params:
        raise ValueError(f"Missing required parameters: {', '.join(missing_params)}")
    
    # Read template file
    try:
        with open(templatePath, 'r') as f:
            template = f.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"Constants template file not found: {templatePath}")
    
    # Replace placeholders with values
    output = template
    for placeholder, (value, _) in constantsMapping.items():
        # Handle different types appropriately
        if isinstance(value, bool):
            value = '1' if value else '0'
        elif isinstance(value, float):
            value = f"{value:.15g}"  # Use g format for clean float representation
        else:
            value = str(value)
            
        output = output.replace(placeholder, value)
    
    # Write output file
    try:
        with open(outputPath, 'w') as f:
            f.write(output)
    except IOError:
        raise IOError(f"Failed to write constants file: {outputPath}")
    
    return True