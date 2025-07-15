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


constexpr int nchi = 0;  // Option to compute drag coefficients chi and chi2

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
constexpr int nbfr = 1;            // Number of boundary forcing frequencies

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

constexpr int nope  = 1; // Number of open boundary segments
constexpr int mnond = 13; // Maximum number of nodes per open boundary segment
constexpr int nland = 1; // Number of land boundary segments
constexpr int mnlnd = 31; // Maximum number of nodes per land boundary segment

constexpr double rearth_eq   = 6378206.4; // [m] 
constexpr double rearth_pole = 6378206.4; // [m]

constexpr double thetai = 0.6;

constexpr double rmaxvel = 10; // Maximum velocity for the x-component
constexpr double rmaxvelx = rmaxvel;       // Maximum velocity for the x-component
constexpr double rmaxvely = rmaxvel*1.013; // Maximum velocity for the y-component

constexpr double rho0 = 1000; // [kg/m^3]
constexpr double dtb_min = 10; //XXXX
constexpr double dtb_max = 30; //XXXX
constexpr double velmin_btrack = 0.0001; //XXXX

constexpr double dzb_min = 0.5; // Min. bottom boundary layer thickness [m]

constexpr double dt = 300;  // Time step  XXX 

constexpr double h0 = 0.01;  // Minimum water depth [m]

enum class VarCode {x=0, y=1, z=2, u=3, v=4, w=5, nx=6, ny=7, hhat=8, ssign=9, btu=10, btv=11, bccx=12, bccy=13, uth=14, vth=15, ubc=16, vbc=17, island=18, ifltype=19, elbc=20, eta=21, area=22, qel=23, depth=24, radius=25, xlon=26, ylat=27, kb=28, n1=29, n2=30}; // Var codes for getElementVar, getNodeVar...
