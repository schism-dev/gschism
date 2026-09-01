from python.schism_consts import *

def processConstants(namelist, tidalData, ivcor: int, nvrt: int, nzlev: int, nsiglev: int,
                     nope: int, mnond: int, nland: int, mnlnd: int,
                     numPartNodes: list[int],
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
    nzlev : int
        Number of z-levels
    nsiglev : int
        Number of sigma levels
    templatePath : str
        Path to constants_template.h
    outputPath : str
        Path where constants.h will be generated
    """
   
    if ivcor != 1:
        useLocalizedSigma = 'DO_NOT_USE_LOCALIZED_SIGMA'
    else:
        useLocalizedSigma = 'USE_LOCALIZED_SIGMA'

    if namelist.opt.ics == 1:
        useFrames = 'DO_NOT_USE_FRAMES'
    elif namelist.opt.ics == 2:
        useFrames = 'USE_FRAMES'
    else:
        raise ValueError(f"Invalid ics = {namelist.opt.ics}, not in [1, 2]")

    # Turbulence submodel and stability model
    turbSubmodelMap = {'MY': 0, 'KL': 1, 'KE': 2, 'KW': 3, 'UB': 4}
    stabilityModelMap = {'GA': 0, 'KC': 1}
    turbSubmodelVal = -1      # TurbulenceSubmodel::NotSet
    stabilityModelVal = -1    # StabilityModel::NotSet
    if hasattr(namelist.opt, 'mid') and namelist.opt.mid is not None:
        mid = str(namelist.opt.mid).strip().upper()
        if mid not in turbSubmodelMap:
            raise ValueError(f"Invalid mid = {namelist.opt.mid!r}, expected one of {list(turbSubmodelMap)}")
        turbSubmodelVal = turbSubmodelMap[mid]
    if hasattr(namelist.opt, 'stab') and namelist.opt.stab is not None:
        stab = str(namelist.opt.stab).strip().upper()
        if stab not in stabilityModelMap:
            raise ValueError(f"Invalid stab = {namelist.opt.stab!r}, expected one of {list(stabilityModelMap)}")
        stabilityModelVal = stabilityModelMap[stab]

    nGPUs = len(numPartNodes)  # or deviceCount
    offsets = [0]
    for n in numPartNodes[:-1]:
        offsets.append(offsets[-1] + int(n))

    # Define mapping of template placeholders to values
    # Format: 'PLACEHOLDER': (value, source_description)
    constantsMapping = {
        # From namelist CORE section
        '__USE_LOCALIZED_SIGMA__': (useLocalizedSigma, 'namelist.core'),
        '__USE_FRAMES__': (useFrames, 'namelist.core'),
        '__BLOCK_SIZE__': (blockSize, 'cuda'),
        '__N_TRACERS__': (namelist.core.ntracer_gen, 'namelist.core'),
        '__IBC__': (namelist.core.ibc, 'namelist.core'),
        '__DT__': (namelist.core.dt, 'namelist.core'),
        
        # From namelist OPT section
        '__ICS__': (namelist.opt.ics, 'namelist.opt'),
        '__IMM__': (namelist.opt.imm, 'namelist.opt'),
        '__INDVEL__': (namelist.opt.indvel, 'namelist.opt'),
        '__IVEG__': (namelist.opt.iveg, 'namelist.opt'),
        '__IHHAT__': (namelist.opt.ihhat, 'namelist.opt'),
        '__IHYDRAULICS__': (namelist.opt.ihydraulics, 'namelist.opt'),
        '__NCHI__': (namelist.opt.nchi, 'namelist.opt'),
        '__HMIN_MAN__': (namelist.opt.hmin_man, 'namelist.opt'), 
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
        '__RMAXVEL__': (namelist.opt.rmaxvel, 'namelist.opt'),
        '__RHO0__': (namelist.opt.rho0, 'namelist.opt'),
        '__DFV0__': (namelist.opt.dfv0, 'namelist.opt'),
        '__DTB_MIN__': (namelist.opt.dtb_min, 'namelist.opt'),
        '__DTB_MAX__': (namelist.opt.dtb_max, 'namelist.opt'),
        '__NADV__': (namelist.opt.nadv, 'namelist.opt'),
        '__VELMIN_BTRACK__': (namelist.opt.velmin_btrack, 'namelist.opt'),
        '__DZB_MIN__': (namelist.opt.dzb_min, 'namelist.opt'),
        '__H0__': (namelist.opt.h0, 'namelist.opt'),
        '__ITUR__': (namelist.opt.itur, 'namelist.opt'),
        '__ISLIP__': (namelist.opt.islip, 'namelist.opt'),
        '__IUPWIND_MOM__': (namelist.opt.iupwind_mom, 'namelist.opt'),
        '__IHORCON__': (namelist.opt.ihorcon, 'namelist.opt'),
        '__HVIS_COEF0__': (namelist.opt.hvis_coef0, 'namelist.opt'),
        '__ISHAPIRO__': (namelist.opt.ishapiro, 'namelist.opt'),
        '__SHAPIRO0__': (namelist.opt.shapiro0, 'namelist.opt'),
        '__NITER_SHAP__': (namelist.opt.niter_shap, 'namelist.opt'),
        '__NCOR__': (namelist.opt.ncor, 'namelist.opt'),
        '__CORICOEF__': (namelist.opt.coricoef, 'namelist.opt'),
        '__SFEA0__': (namelist.opt.sfea0, 'namelist.opt'),
        '__NWS__': (namelist.opt.nws, 'namelist.opt'),
        '__WTIMINC__': (namelist.opt.wtiminc, 'namelist.opt'),
        '__IWIND_FORM__': (namelist.opt.iwind_form, 'namelist.opt'),
        '__IWINDOFF__': (namelist.opt.iwindoff, 'namelist.opt'),
        '__DRAMPWIND__': (namelist.opt.drampwind, 'namelist.opt'),
        '__IHCONSV__': (namelist.opt.ihconsv, 'namelist.opt'),
        '__ISCONSV__': (namelist.opt.isconsv, 'namelist.opt'),
        # From grid parameters
        '__IVCOR__': (ivcor, 'grid'),
        '__NVRT__': (nvrt, 'grid'),
        '__NZLEV__': (nzlev, 'grid'), 
        '__NSIGLEV__': (nsiglev, 'grid'),
        '__NOPEN_BOUNDARY_SEGMENTS__': (nope, 'grid'),
        '__MAX_NODES_PER_OPEN_BOUNDARY_SEGMENT__': (mnond, 'grid'),
        '__N_LAND_BOUNDARY_SEGMENTS__': (nland, 'grid'),
        '__MAX_NODES_PER_LAND_BOUNDARY_SEGMENT__': (mnlnd, 'grid'),
        '__N_GPUS__': (nGPUs, 'partition'),
        '__NUM_NODES_PER_PARTITION__': (', '.join(str(int(n)) for n in numPartNodes), 'partition'),
        '__NODE_OFFSET_PER_PARTITION__': (', '.join(str(o) for o in offsets), 'partition'),

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


    # Add var codes for getElementVar, getNodeVar...
    # Create enum string with integer values

    enum_str = ', '.join([f"{code}={i}" for i, code in enumerate(varInfoDict.keys())])
    output = output.replace("__VARCODES__", enum_str)
    
    # Write output file
    try:
        with open(outputPath, 'w') as f:
            f.write(output)
    except IOError:
        raise IOError(f"Failed to write constants file: {outputPath}")
    
    return True