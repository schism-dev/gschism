#define USE_LOCALIZED_SIGMA
#define USE_FRAMES

constexpr double deg2rad = 3.14159265358979323846/180.0;
constexpr double rad2deg = 180.0/3.14159265358979323846;

constexpr int BLOCK_SIZE = 128;

constexpr int N_TRACERS = 2;

constexpr int ivcor = 1;              // Vertical coordinate system
constexpr int ics = 2;                 // 1: Local coordinates, 2: 3D global cartesian coordinates
constexpr int imm = 0;                 // 0: without bed deformation; 1: with bed deformation (e.g., tsunami); 2: 3D bed deformation model (needs user coding)
constexpr int indvel = 0;           // Option to compute nodal velocities from element velocities
constexpr int iveg = 0;               // Vegetation option, 0: no vegetation, 1: vegetation
constexpr int ihhat = 1;             // Option to limit \hat{H} to enhance stability for large friction in shallow area
constexpr int ihydraulics = 0; // Option of including hydraulics
constexpr int nhtblocks = 0;                 // Number of hydraulic structures (zero for the moment)
constexpr int ibc = 1;                 // Option to include baroclinic term (0 means yes)
constexpr int itur = 0;               // Turbulence closure model
constexpr int islip = 0;             // Option to include slip horizontal boundary conditions
constexpr int ihorcon = 0;         // Option to include horizontal viscosity
constexpr double hvis_coef0 = 0.025;// Horizontal viscosity coefficient
constexpr int ishapiro = 1;       // Option to include Shapiro filter
constexpr double shapiro0 = 0.5;    // Shapiro filter parameter
constexpr int niter_shap = 1;   // Number of iterations for Shapiro filter
constexpr int ncor = 1;               // Option to include Coriolis force
constexpr double coricoef = 0;    // Coriolis force coefficient
constexpr double sfea0 = 24*deg2rad; // Reference latitude for Coriolis force
constexpr int nws = 2;               // Option to include wind forcing
constexpr double wtiminc = 3600;    // Time increment for wind forcing
constexpr int iwind_form = 1; // Option to include wind forcing
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

constexpr int ntip = 8;            // Number of tidal potential components
constexpr double tip_dp = 50;   // cut-off depth for applying tidal potential

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

constexpr int nvrt    = 50; // Number of vertical levels
constexpr int nzlev   = 0; // Number of z levels (in contrast to sigma levels)
constexpr int nsiglev = 50; // Number of sigma levels (in contrast to z levels)

//  To simplify comparisons, I define index versions of count variables.
constexpr int ivrt    = nvrt - 1;
constexpr int izlev   = nzlev - 1;
constexpr int isiglev = nsiglev - 1;

constexpr int nope  = 1; // Number of open boundary segments
constexpr int mnond = 1199; // Maximum number of nodes per open boundary segment
constexpr int nland = 43; // Number of land boundary segments
constexpr int mnlnd = 28091; // Maximum number of nodes per land boundary segment

constexpr double rearth_eq   = 6378206.4; // [m] 
constexpr double rearth_pole = 6378206.4; // [m]

constexpr double thetai = 1;

constexpr double rmaxvel = 4; // Maximum velocity for the x-component
constexpr double rmaxvelx = rmaxvel;       // Maximum velocity for the x-component
constexpr double rmaxvely = rmaxvel*1.013; // Maximum velocity for the y-component

constexpr double rho0 = 1000; // [kg/m^3]
constexpr double dtb_min = 10; //XXXX
constexpr double dtb_max = 30; //XXXX
constexpr double velmin_btrack = 0.0001; // Minimum velocity for backtracking (m/s)
constexpr double velmin_btrack_sq = velmin_btrack*velmin_btrack;
constexpr double btrack_nudge = 0.009013;   // Nudging factor for backtracking

constexpr int nadv = 1; // Advection flag

constexpr double dfv0 = 0.01; 

constexpr double dzb_min = 0.5; // Min. bottom boundary layer thickness [m]

constexpr double dt = 120;  // Time step  XXX 

constexpr double h0 = 0.0001;  // Minimum water depth for wetting and drying [m]

enum class VarCode {x=0, z=1, u=2, w=3, U=4, nx=5, ny=6, hhat=7, ghat=8, chi=9, ssign=10, btu=11, bcc=12, uth=13, ubc=14, island=15, ifltype=16, elbc=17, isdry=18, eta=19, area=20, qel=21, depth=22, radius=23, len=24, lonlat=25, kb=26, n1=27, n2=28, idx=29, htot=30, gradPr=31, gradEta1=32, gradEta2=33, gradETP=34, FaceIndices=35, u_air=36, outputCoords=37}; // Var codes for getElementVar, getNodeVar...
