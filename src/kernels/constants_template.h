#define __USE_LOCALIZED_SIGMA__
#define __USE_FRAMES__

constexpr double deg2rad = 3.14159265358979323846/180.0;
constexpr double rad2deg = 180.0/3.14159265358979323846;

constexpr int BLOCK_SIZE = __BLOCK_SIZE__;

constexpr int N_TRACERS = __N_TRACERS__;

constexpr int ivcor = __IVCOR__;              // Vertical coordinate system
constexpr int ics = __ICS__;                 // 1: Local coordinates, 2: 3D global cartesian coordinates
constexpr int imm = __IMM__;                 // 0: without bed deformation; 1: with bed deformation (e.g., tsunami); 2: 3D bed deformation model (needs user coding)
constexpr int indvel = __INDVEL__;           // Option to compute nodal velocities from element velocities
constexpr int iveg = __IVEG__;               // Vegetation option, 0: no vegetation, 1: vegetation
constexpr int ihhat = __IHHAT__;             // Option to limit \hat{H} to enhance stability for large friction in shallow area
constexpr int ihydraulics = __IHYDRAULICS__; // Option of including hydraulics
constexpr int nhtblocks = 0;                 // Number of hydraulic structures (zero for the moment)
constexpr int ibc = __IBC__;                 // Option to include baroclinic term (0 means yes)
constexpr int itur = __ITUR__;               // Turbulence closure model
constexpr int islip = __ISLIP__;             // Option to include slip horizontal boundary conditions
constexpr int ihorcon = __IHORCON__;         // Option to include horizontal viscosity
constexpr double hvis_coef0 = __HVIS_COEF0__;// Horizontal viscosity coefficient
constexpr int ishapiro = __ISHAPIRO__;       // Option to include Shapiro filter
constexpr double shapiro0 = __SHAPIRO0__;    // Shapiro filter parameter
constexpr int niter_shap = __NITER_SHAP__;   // Number of iterations for Shapiro filter
constexpr int ncor = __NCOR__;               // Option to include Coriolis force
constexpr double coricoef = __CORICOEF__;    // Coriolis force coefficient
constexpr double sfea0 = __SFEA0__*deg2rad; // Reference latitude for Coriolis force
constexpr int nws = __NWS__;               // Option to include wind forcing
constexpr double wtiminc = __WTIMINC__;    // Time increment for wind forcing
constexpr int iwind_form = __IWIND_FORM__; // Option to include wind forcing
constexpr int iwindoff = __IWINDOFF__;     // Option to include wind forcing
constexpr double drampwind = __DRAMPWIND__; // Option to include wind forcing

constexpr int ihconsv = __IHCONSV__;     // Option to include conservative wind forcing
constexpr int isconsv = __ISCONSV__;     // Option to include conservative wind forcing


constexpr int nchi = __NCHI__;             // Option to compute drag coefficients chi and chi2
constexpr double hmin_man = __HMIN_MAN__;  // Minimum water depth for Manning's friction

constexpr int nramp_elev = __NRAMP_ELEV__; // Option to use ramp for elevation boundary conditions

constexpr int    iloadtide = __ILOADTIDE__;       // Option to use polynomial for discharge boundary
constexpr double loadtide_coef = __LOADTIDE_COEF__; // Load tide coefficient

constexpr int    inv_atm_bnd = __INV_ATM_BND__;     // Option to use inverse barometric effects at boundary
constexpr double prmsl_ref = __PRMSL_REF__; // Reference pressure for inverse barometric effects

constexpr double veg_lai = __VEG_LAI__; // Vegetation drag coefficient
constexpr double veg_cw = __VEG_CW__; // Vegetation drag coefficient

constexpr int ntip = __NTIP__;            // Number of tidal potential components
constexpr double tip_dp = __TIP_DP__;   // cut-off depth for applying tidal potential

// Boundary forcing frequencies
constexpr int nbfr = __NBFR__;            // Number of boundary forcing frequencies

// Boundary condition types
constexpr int nettype = __NETTYPE__;      // Total number of type I elevation boundaries
constexpr int nfltype = __NFLTYPE__;      // Total number of type I flow boundaries
constexpr int nettype2 = __NETTYPE2__;    // Total number of type IV elevation boundaries
constexpr int nfltype2 = __NFLTYPE2__;    // Total number of type IV flow boundaries

// Boundary nodes
constexpr int nnode_et = __NNODE_ET__;    // Number of nodes requiring elev2D.th.nc
constexpr int nnode_fl = __NNODE_FL__;    // Number of nodes requiring flow2D.th.nc

constexpr int iunder_deep = __IUNDER_DEEP__;     // Option of under-deep methods
constexpr double h1_bcc   = __H1_BCC__;  // Parameters for under-deep methods
constexpr double h2_bcc   = __H2_BCC__; 
constexpr double hw_depth = __HW_DEPTH__; 
constexpr double hw_ratio = __HW_RATIO__;

constexpr int nvrt    = __NVRT__; // Number of vertical levels
constexpr int nzlev   = __NZLEV__; // Number of z levels (in contrast to sigma levels)
constexpr int nsiglev = __NSIGLEV__; // Number of sigma levels (in contrast to z levels)

//  To simplify comparisons, I define index versions of count variables.
constexpr int ivrt    = nvrt - 1;
constexpr int izlev   = nzlev - 1;
constexpr int isiglev = nsiglev - 1;

constexpr int nope  = __NOPEN_BOUNDARY_SEGMENTS__; // Number of open boundary segments
constexpr int mnond = __MAX_NODES_PER_OPEN_BOUNDARY_SEGMENT__; // Maximum number of nodes per open boundary segment
constexpr int nland = __N_LAND_BOUNDARY_SEGMENTS__; // Number of land boundary segments
constexpr int mnlnd = __MAX_NODES_PER_LAND_BOUNDARY_SEGMENT__; // Maximum number of nodes per land boundary segment

constexpr double rearth_eq   = __REARTH_EQ__; // [m] 
constexpr double rearth_pole = __REARTH_POLE__; // [m]

constexpr double thetai = __THETAI__;

constexpr double rmaxvel = __RMAXVEL__; // Maximum velocity for the x-component
constexpr double rmaxvelx = rmaxvel;       // Maximum velocity for the x-component
constexpr double rmaxvely = rmaxvel*1.013; // Maximum velocity for the y-component

constexpr double rho0 = __RHO0__; // [kg/m^3]
constexpr double dtb_min = __DTB_MIN__; //XXXX
constexpr double dtb_max = __DTB_MAX__; //XXXX
constexpr double velmin_btrack = __VELMIN_BTRACK__; // Minimum velocity for backtracking (m/s)
constexpr double btrack_nudge = __BTRACK_NUDGE__;   // Nudging factor for backtracking

constexpr int nadv = __NADV__; // Advection flag

constexpr double dfv0 = __DFV0__; 

constexpr double dzb_min = __DZB_MIN__; // Min. bottom boundary layer thickness [m]

constexpr double dt = __DT__;  // Time step  XXX 

constexpr double h0 = __H0__;  // Minimum water depth for wetting and drying [m]

enum class VarCode {__VARCODES__}; // Var codes for getElementVar, getNodeVar...
