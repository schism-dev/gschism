constexpr int BLOCK_SIZE = 256;

constexpr int N_TRACERS = 2;

constexpr int ivcor = 2;              // Vertical coordinate system
constexpr int ics = 1;                 // 1: global, 2: local
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
constexpr double sfea0 = 24;          // Reference latitude for Coriolis force


constexpr int nchi = 0;  // Option to compute drag coefficients chi and chi2

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

constexpr int nvrt    = 2; // Number of vertical levels
constexpr int nzlev   = 1; // Number of z levels (in contrast to sigma levels)
constexpr int nsiglev = 2; // Number of sigma levels (in contrast to z levels)

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
constexpr double velmin_btrack = 0.0001; //XXXX

constexpr double dzb_min = 0.5; // Min. bottom boundary layer thickness [m]

constexpr double dt = 120;  // Time step  XXX 

constexpr double h0 = 0.0001;  // Minimum water depth [m]

enum class VarCode {x=0, y=1, z=2, u=3, v=4, w=5, U=6, V=7, nx=8, ny=9, hhat=10, ghatx=11, ghaty=12, ssign=13, btu=14, btv=15, bccx=16, bccy=17, uth=18, vth=19, ubc=20, vbc=21, island=22, ifltype=23, elbc=24, eta=25, area=26, qel=27, depth=28, radius=29, xlon=30, ylat=31, kb=32, n1=33, n2=34, htot=35, gradEta1x=36, gradEta1y=37, gradEta2x=38, gradEta2y=39, gradETPx=40, gradETPy=41, gradPrx=42, gradPry=43, idx=44}; // Var codes for getElementVar, getNodeVar...
