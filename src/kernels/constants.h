#define USE_LOCALIZED_SIGMA
#define USE_FRAMES

constexpr double deg2rad = 3.14159265358979323846/180.0;
constexpr double rad2deg = 180.0/3.14159265358979323846;

constexpr int BLOCK_SIZE = 128;

constexpr int N_TRACERS = 2;

constexpr int nGPUs = 1; // nGPUs is always equal to the number of partitions.
constexpr int nPartitions = nGPUs;

// Turbulence model and submodel
enum class TurbulenceModel : int {
    VertVarying  = -1,
    HorzVarying  = -2,
    Constant     =  0,
    Pacanowski   =  2,
    GLS          =  3,
    GOTM         =  4,
    SedimentMix  =  5,
};

enum class TurbulenceSubmodel : int {
    TurbNotSet     = -1,
    MellorYamada   = 0,
    KL             = 1,
    KE             = 2,
    KW             = 3,
    UmlaufBurchard = 4,
};

enum class StabilityModel : int {
    StabNotSet    = -1,
    Galperin      = 0,   // stab = 'GA'
    KanthaClayson = 1,   // stab = 'KC'
};


constexpr TurbulenceModel turbulenceModel = static_cast<TurbulenceModel>(3);
constexpr TurbulenceSubmodel turbSubmodel = static_cast<TurbulenceSubmodel>(1);
constexpr StabilityModel stabilityModel   = static_cast<StabilityModel>(1);


// Nodes owned by each partition
constexpr int numNodesPerPartition[nPartitions] = {20641};
// Cumulative counts: offset[0]=numNodesPerPartition[0], offset[p]=offset[p-1]+numNodesPerPartition[p]
constexpr int nodeOffsetPerPartition[nPartitions] = {0};



constexpr int ivcor = 1;              // Vertical coordinate system
constexpr int ics = 2;                 // 1: Local coordinates, 2: 3D global cartesian coordinates
constexpr int imm = 0;                 // 0: without bed deformation; 1: with bed deformation (e.g., tsunami); 2: 3D bed deformation model (needs user coding)
constexpr int indvel = 0;           // Option to compute nodal velocities from element velocities
constexpr int iveg = 0;               // Vegetation option, 0: no vegetation, 1: vegetation
constexpr int ihhat = 1;             // Option to limit \hat{H} to enhance stability for large friction in shallow area
constexpr int ihydraulics = 0; // Option of including hydraulics
constexpr int nhtblocks = 0;                 // Number of hydraulic structures (zero for the moment)
constexpr int ibc = 1;                 // Option to include baroclinic term (0 means yes)
constexpr int itur = 3;               // Turbulence closure model
constexpr int icompute_cpsi3 = 0;  // 0: use tabulated cpsi3; !=0: compute (not implemented)
constexpr int islip = 0;             // Option to include slip horizontal boundary conditions
constexpr int iupwind_mom = 0; // Option on how to include upwind momentum advection
constexpr int ihorcon = 0;         // Option to include horizontal viscosity
constexpr double hvis_coef0 = 0.025;// Horizontal viscosity coefficient
constexpr int ishapiro = 1;       // Option to include Shapiro filter
constexpr double shapiro0 = 0.5;    // Shapiro filter parameter
constexpr int niter_shap = 1;   // Number of iterations for Shapiro filter
constexpr int ncor = 1;               // Option to include Coriolis force
constexpr double coricoef = 0;    // Coriolis force coefficient
constexpr double sfea0 = 46.25*deg2rad; // Reference latitude for Coriolis force
constexpr int nws = 2;               // Option to include wind forcing
constexpr double wtiminc = 90;    // Time increment for wind forcing
constexpr int iwind_form = 0; // Option to include wind forcing
constexpr int iwindoff = 0;     // Option to include wind forcing
constexpr double drampwind = 1; // Option to include wind forcing

constexpr int ihconsv = 0;     // Option to include conservative wind forcing
constexpr int isconsv = 0;     // Option to include conservative wind forcing


constexpr int nchi = 0;             // Option to compute drag coefficients chi and chi2
constexpr double hmin_man = 1;  // Minimum water depth for Manning's friction

constexpr int nramp_elev = 0; // Option to use ramp for elevation boundary conditions

constexpr int    iloadtide = 0;       // Option to use polynomial for discharge boundary
constexpr double loadtide_coef = 0.1; // Load tide coefficient

constexpr int    inv_atm_bnd = 0;     // Option to use inverse barometric effects at boundary
constexpr double prmsl_ref = 101325; // Reference pressure for inverse barometric effects

constexpr double veg_lai = 1; // Vegetation drag coefficient
constexpr double veg_cw = 1.5; // Vegetation drag coefficient

constexpr int ntip = 0;            // Number of tidal potential components
constexpr double tip_dp = 40;   // cut-off depth for applying tidal potential

// Boundary forcing frequencies
constexpr int nbfr = 8;            // Number of boundary forcing frequencies

// Boundary condition types
constexpr int nettype = 0;      // Total number of type I elevation boundaries
constexpr int nfltype = 0;      // Total number of type I flow boundaries
constexpr int nettype2 = 0;    // Total number of type IV elevation boundaries
constexpr int nfltype2 = 0;    // Total number of type IV flow boundaries

// Boundary nodes
constexpr int nnode_et = 0;    // Number of nodes requiring elev2D.th.nc
constexpr int nnode_fl = 0;    // Number of nodes requiring flow2D.th.nc

constexpr int iunder_deep = 0;     // Option of under-deep methods
constexpr double h1_bcc   = 50;  // Parameters for under-deep methods
constexpr double h2_bcc   = 100; 
constexpr double hw_depth = 1000000; 
constexpr double hw_ratio = 0.5;

constexpr int nvrt    = 46; // Number of vertical levels
constexpr int nzlev   = 0; // Number of z levels (in contrast to sigma levels)
constexpr int nsiglev = 46; // Number of sigma levels (in contrast to z levels)

//  To simplify comparisons, I define index versions of count variables.
constexpr int ivrt    = nvrt - 1;
constexpr int izlev   = nzlev - 1;
constexpr int isiglev = nsiglev - 1;

constexpr int nope  = 4; // Number of open boundary segments
constexpr int mnond = 88; // Maximum number of nodes per open boundary segment
constexpr int nland = 19; // Number of land boundary segments
constexpr int mnlnd = 875; // Maximum number of nodes per land boundary segment

constexpr double rearth_eq   = 6378206.4; // [m] 
constexpr double rearth_pole = 6378206.4; // [m]

constexpr double thetai = 0.502;

constexpr double rmaxvel = 10; // Maximum velocity for the x-component
constexpr double rmaxvelx = rmaxvel;       // Maximum velocity for the x-component
constexpr double rmaxvely = rmaxvel*1.013; // Maximum velocity for the y-component

constexpr double rho0 = 1000; // [kg/m^3]
constexpr double dtb_min = 15; 
constexpr double dtb_max = 30; 
constexpr double velmin_btrack = 0.0001; // Minimum velocity for backtracking (m/s)
constexpr double velmin_btrack_sq = velmin_btrack*velmin_btrack;

constexpr int nadv = 1; // Advection flag

constexpr double dfv0 = 0.01; // turbulent diffusion for momentum
constexpr double dfh0 = 0.0001; // turbulent diffusion for tracers
constexpr double xlsc0 = 0.1; // surface/bottom mixing length scale (itur=3,5)

constexpr double dzb_min = 0.5; // Min. bottom boundary layer thickness [m]

constexpr double dt = 90;  // Time step  XXX 

constexpr double h0 = 0.01;  // Minimum water depth for wetting and drying [m]

enum class VarCode {x=0, z=1, u=2, w=3, U=4, nx=5, ny=6, hhat=7, ghat=8, chi=9, ssign=10, btu=11, bcc=12, uth=13, ubc=14, dfh=15, dfv=16, island=17, ifltype=18, elbc=19, isdry=20, eta=21, area=22, qel=23, depth=24, radius=25, len=26, lonlat=27, kb=28, n1=29, n2=30, idx=31, htot=32, gradPr=33, gradEta1=34, gradEta2=35, gradETP=36, FaceIndices=37, u_air=38, outputCoords=39}; // Var codes for getElementVar, getNodeVar...
