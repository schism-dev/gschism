// BEGIN INCLUDE: constants.h
#define USE_LOCALIZED_SIGMA
#define DO_NOT_USE_FRAMES

constexpr double deg2rad = 3.14159265358979323846/180.0;
constexpr double rad2deg = 180.0/3.14159265358979323846;

constexpr int BLOCK_SIZE = 128;

constexpr int N_TRACERS = 2;

constexpr int nGPUs = 1; // nGPUs is always equal to the number of partitions.
constexpr int nPartitions = nGPUs;

// Nodes owned by each partition
constexpr int numNodesPerPartition[nPartitions] = {1457889};
// Cumulative counts: offset[0]=numNodesPerPartition[0], offset[p]=offset[p-1]+numNodesPerPartition[p]
constexpr int nodeOffsetPerPartition[nPartitions] = {1457889};

constexpr int ivcor = 1;              // Vertical coordinate system
constexpr int ics = 1;                 // 1: Local coordinates, 2: 3D global cartesian coordinates
constexpr int imm = 0;                 // 0: without bed deformation; 1: with bed deformation (e.g., tsunami); 2: 3D bed deformation model (needs user coding)
constexpr int indvel = 0;           // Option to compute nodal velocities from element velocities
constexpr int iveg = 0;               // Vegetation option, 0: no vegetation, 1: vegetation
constexpr int ihhat = 1;             // Option to limit \hat{H} to enhance stability for large friction in shallow area
constexpr int ihydraulics = 0; // Option of including hydraulics
constexpr int nhtblocks = 0;                 // Number of hydraulic structures (zero for the moment)
constexpr int ibc = 1;                 // Option to include baroclinic term (0 means yes)
constexpr int itur = 0;               // Turbulence closure model
constexpr int islip = 0;             // Option to include slip horizontal boundary conditions
constexpr int iupwind_mom = 0; // Option on how to include upwind momentum advection
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
constexpr double dtb_min = 10; 
constexpr double dtb_max = 30; 
constexpr double velmin_btrack = 0.0001; // Minimum velocity for backtracking (m/s)
constexpr double velmin_btrack_sq = velmin_btrack*velmin_btrack;

constexpr int nadv = 1; // Advection flag

constexpr double dfv0 = 0.01; 

constexpr double dzb_min = 0.5; // Min. bottom boundary layer thickness [m]

constexpr double dt = 120;  // Time step  XXX 

constexpr double h0 = 0.0001;  // Minimum water depth for wetting and drying [m]

enum class VarCode {x=0, z=1, u=2, w=3, U=4, nx=5, ny=6, hhat=7, ghat=8, chi=9, ssign=10, btu=11, bcc=12, uth=13, ubc=14, island=15, ifltype=16, elbc=17, isdry=18, eta=19, area=20, qel=21, depth=22, radius=23, len=24, lonlat=25, kb=26, n1=27, n2=28, idx=29, htot=30, gradPr=31, gradEta1=32, gradEta2=33, gradETP=34, FaceIndices=35, u_air=36, outputCoords=37}; // Var codes for getElementVar, getNodeVar...
// END INCLUDE: constants.h
// BEGIN INCLUDE: meshdata.h
using int16_t = short;
using uintptr_t = unsigned long long;   // CUDA device pointers are 64-bit

constexpr double small1 = 1.0e-7;       // Small non-negative number
constexpr double small2 = small1*100.0;  // Slightly larger than small1
constexpr double small = 1.0e-7;            // Small positive number or 0

constexpr int internalSegment = -1;  // It is not really a segment, but we keep the name.
constexpr int landSegment     = -2;

constexpr double pi = 3.14159265358979323846;

constexpr double eps = 1.019e-2; // Nudging factor XXX WHY THIS VALUE?

constexpr double grav        = 9.81;      // [m/s^2]
constexpr double omega_e     = 7.292e-5;  // [1/s]

// eta coordinates of the nodes (idxNode = 0,1,2,3)
__device__ __forceinline__ double eta_n(int idxNode) {
   // bit1: 0 -> -1, 1 -> +1
   assert(idxNode >= 0 && idxNode < 4 && "etaSign: idxNode out of range");
   return (idxNode & 2) ? 1.0 : -1.0;     // [-1,-1,+1,+1]
}

// xi coordinates of the nodes (idxNode = 0,1,2,3)
__device__ __forceinline__ double xi_n(int idxNode) {
   // +1 if (bit0 XOR bit1), else -1
   assert(idxNode >= 0 && idxNode < 4 && "xiSign: idxNode out of range");
   return ((idxNode & 1) ^ (idxNode >> 1)) ? 1.0 : -1.0;   // [-1,+1,+1,-1]

}

#define INT_MAX 2147483647



// Time ramps ---------------------------------------------------------------------------

__device__ double rampbc;
__device__ double rampwind;
__device__ double rampwafo;
__device__ double ramp_ss;
__device__ double ramp;


// Tidal Harmonic constituents ---------------------------------------------------------------------------

__device__ double *amig; // Frequency of the harmonic constituent
__device__ double *face; // Nodal factor of the harmonic constituent
__device__ double *ff;   // Earth equilibrium argument of the harmonic constituent



struct Bool4 {
   unsigned int m; // bit i = lane i

   __device__ __forceinline__ bool get(int i) const { return (m >> i) & 1u; }
   __device__ __forceinline__ void set(int i, bool v) {m = (m & ~(1u << i)) | ((unsigned int)(v) << i); }

   __device__ __forceinline__ bool all(void) const { return m == 0x0F; }
   __device__ __forceinline__ bool any(void) const { return m != 0; }
   __device__ __forceinline__ bool none(void) const { return m == 0; }
   __device__ __forceinline__ int firstOne(void) const { return __ffs(m) - 1; }

};

__device__ __forceinline__ Bool4 operator!(Bool4 x) { x.m = (~x.m) & 0x0F; return x; }
__device__ __forceinline__ Bool4 operator&(Bool4 a, Bool4 b) { return Bool4{(unsigned int)(a.m & b.m)}; }
__device__ __forceinline__ Bool4 operator|(Bool4 a, Bool4 b) { return Bool4{(unsigned int)(a.m | b.m)}; }


// Union to store different types of variables in the same memory location. This is used to define temporary 
// variables that can be used to store different types of variables.
union TempVar {
   double d;
   int i;
   bool b;
};


extern "C" __global__ void initTidalHarmonicConstituents(double *_amig, double *_face, double *_ff) {
// Initialize the pointers to the tidal harmonic constituents

   amig = _amig;
   face = _face;
   ff   = _ff;
}

// matrix solver scalars ---------------------------------------------------------------------------
__device__ double *rMr_old;
__device__ double *rMr_new;
__device__ double *pAp;


// Simple math functions ---------------------------------------------------------------------------

__device__ __forceinline__ void solveQuadratic(double A, double B, double C, double &x1, double &x2) {
// Solves the quadratic equation A*x^2 + B*x + C = 0 and returns the two solutions in x1 and x2.
// Outputs:
//        x1: first solution
//        x2: second solution

       // Use float-specific absolute value for speed on GPU
    const double absA = fabs(A);
    const double absB = fabs(B);

    // 1. Near-linear case: a is effectively 0.
    // We use a small epsilon to prevent division by zero and handle degenerate quads.
    // This branch is almost never taken, so it won't hurt performance.
    if (absA < 1e-14 * (absB + 1e-30)) {
        const double x = -C / B;
        x1 = x;
        x2 = x;
        return;  // Degenerate case: A = 0.
    }

   const double delta = max(0.0, B*B - 4.0*A*C);

   assert(delta >= 0.0 && "solveQuadratic: negative discriminant");
   const double s = sqrt(delta);

   // Numerically stable quadratic roots:
   // q = -0.5*(B + sign(B)*sqrt(delta))
   const double q = -0.5 * (B + copysign(s, B));

   x1 = q/A;
   x2 = (fabs(q) > 1e-30) ? (C/q) : x1;   // Moving predicate.

}

template<int grade>
__device__ __forceinline__ double pow(double x);

template<>
__device__ __forceinline__ double pow<2>(double x) {
    return x*x;
}

template<>
__device__ __forceinline__ double pow<3>(double x) {
    return x*x*x;
}


__device__ __forceinline__ double clamp(double x, double minval, double maxval) {
    return fmin(fmax(x, minval), maxval);
}

__device__ __forceinline__ int clamp(int x, int minval, int maxval) {
    return max(minval, min(x, maxval));
}




// VECTOR/POINT operations ---------------------------------------------------------------------------

using Vector2 = double2; 
using Vector3 = double3;
using Point2  = double2;
using Point3  = double3;


template<typename T>
struct IsVector2 { static constexpr bool value = false; };
template<>
struct IsVector2<Vector2> { static constexpr bool value = true; };


__device__ __forceinline__ double2 operator+(const double2& a, const double2& b) {
    return make_double2(a.x + b.x, a.y + b.y);
}



__device__ __forceinline__ double2 operator+=(double2& a, const double2& b) {
    a.x += b.x;
    a.y += b.y;
    return a;
}

__device__ __forceinline__ double2 operator-(const double2& a, const double2& b) {
    return make_double2(a.x - b.x, a.y - b.y);
}


__device__ __forceinline__ double2 operator-(const double2& a) {
    return make_double2(-a.x, -a.y);
}


__device__ __forceinline__ double2 operator-=(double2& a, const double2& b) {
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

__device__ __forceinline__ double2 operator*(const double2& v, const double s) {
    return make_double2(v.x * s, v.y * s);
}

__device__ __forceinline__ double2 operator*(const double s, const double2& v) {
    return v * s;  // reuse the previous operator
}

__device__ __forceinline__ double2 operator/(const double2& v, double s) {
    return make_double2(v.x / s, v.y / s);
}

__device__ __forceinline__ double norm2(const double2& v) {
    return v.x*v.x + v.y*v.y;
}

__device__ __forceinline__ double length(const double2& v) {
    return sqrt(v.x*v.x + v.y*v.y);
}

__device__ __forceinline__ double2 normalize(const double2& v) {
    assert(length(v) > 0.0 && "normalize: length = 0");
    return v / length(v);
}

__device__ __forceinline__ double dot(const double2& v, const double2& w) {
    return v.x*w.x + v.y*w.y;
}

__device__ __forceinline__ double cross(const double2& v, const double2& w) {
    return v.x*w.y - v.y*w.x;
}

__device__ __forceinline__ double2 rot90(const double2& a) {
    return make_double2(a.y, -a.x);
}

__device__ __forceinline__ double3 operator+(const double3& a, const double3& b) {
    return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__device__ __forceinline__ double3 operator+=(double3& a, const double3& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

__device__ __forceinline__ double3 operator-(const double3& a, const double3& b) {
    return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__device__ __forceinline__ double3 operator-(const double3& a) {
    return make_double3(-a.x, -a.y, -a.z);
}

__device__ __forceinline__ double3 operator*(const double3& v, double s) {
    return make_double3(v.x * s, v.y * s, v.z * s);
}

__device__ __forceinline__ double3 operator*(double s, const double3& v) {
    return v * s;
}

__device__ __forceinline__ double3 operator/(const double3& v, double s) {
    return make_double3(v.x / s, v.y / s, v.z / s);
}

__device__ __forceinline__ double dot(const double3& v, const double3& w) {
    return v.x*w.x + v.y*w.y + v.z*w.z;
}

__device__ __forceinline__ double3 cross(const double3& a, const double3& b) {
    return make_double3(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    );
}

__device__ __forceinline__ double norm2(const double3& v) {
    return v.x*v.x + v.y*v.y + v.z*v.z;
}

__device__ __forceinline__ double length(const double3& v) {
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

__device__ __forceinline__ double3 normalize(const double3& v) {
    assert(length(v) > 0.0 && "normalize: length = 0");
    return v / length(v);
}

__device__ __forceinline__ Vector2 toVector2(const double3& v) {
    return make_double2(v.x, v.y);
}

__device__ __forceinline__ double3 toVector3(const double2& v, double w) {
    return make_double3(v.x, v.y, w);
}

__device__ __forceinline__ Point2 toPoint2(const double3& p) {
    return make_double2(p.x, p.y);
}

__device__ __forceinline__ Point3 toPoint3(const double2& p, double z) {
    return make_double3(p.x, p.y, z);
}

__device__ __forceinline__ Vector2& toVector2(double3& v) {
    return *reinterpret_cast<Vector2*>(&v);
}


__device__ __forceinline__ Point2& toPoint2(double3& p) {
    return *reinterpret_cast<Point2*>(&p);
}


// ------------------------------------------------------------------------------------------------

// Polynomial class
template<int grade>
struct Polynomial {
    // Static assertion to give a clear error message at compile time
    static_assert(grade == 4, "Only 4th order polynomials are implemented");
};


// Explicit specialization for grade 4
template<>
struct Polynomial<4> {
   double a0, a1, a2, a3, a4; // Coefficients a0, a1, a2, a3, a4 for a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0

   // Constructor from coefficients
   __device__ __forceinline__ Polynomial(double a0, double a1, double a2, double a3, double a4) : a0(a0), a1(a1), a2(a2), a3(a3), a4(a4) {}

   __device__ __forceinline__ double eval(double x) const {
      return a0 + x*(a1 + x*(a2 + x*(a3 + x*a4)));
   }
};

// ------------------------------------------------------------------------------------------------

__forceinline__ __device__ double signedArea(double x1, double x2, double x3, double y1, double y2, double y3) {
  // Compute signed area formed by pts 1,2,3 (positive counter-clockwise)
  return 0.5*((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3));
}

__forceinline__ __device__ double signedArea(Point2 p1, Point2 p2, Point2 p3) {
  // Compute signed area formed by pts 1,2,3 (positive counter-clockwise)
  return 0.5*((p1.x-p3.x)*(p2.y-p3.y) - (p2.x-p3.x)*(p1.y-p3.y));
}


__forceinline__ __device__ void crossProduct(double x1, double y1, double z1, double x2, double y2, double z2, double &x3, double &y3, double &z3) {
   x3 = y1*z2 - y2*z1;
   y3 = x2*z1 - x1*z2;
   z3 = x1*y2 - x2*y1;
// XXXXXX IS THIS REPEATED????
}


__forceinline__ __device__ void computeLatLon(Point3 const &p, double &rlat, double &rlon) {
// Given a 3D global point (may not be on the surface of the Earth), find its longitude and latitude in radians
//

   const double radius2 = dot(p, p);
   const double radius = sqrt(radius2);
 
   assert(radius != 0.0      && "computeLatLon: radius = 0");
   assert(p.z*p.z <= radius2 && "computeLatLon: p.z <= radius");

   rlon = atan2(p.y, p.x); // (-pi,pi]

   if (abs(rearth_pole - rearth_eq) < 1.0e-2) { // for backward compatibility
      rlat = asin(p.z/radius);
   } else {
      rlat = asin(p.z/rearth_pole);
   }
  
}


__forceinline__ __device__ Point3 computeXYZfromLatLon(double lat, double lon) {
// Compute the 3D Cartesian coordinates from the latitude and longitude
// Inputs:
//        lat: latitude in radians
//        lon: longitude in radians
// Output:
//        x, y, z: 3D Cartesian coordinates

   double sin_lat, cos_lat, sin_lon, cos_lon;

   assert(lat >= -pi/2.0 && lat <= pi/2.0 && "computeXYZfromLatLon: lat out of range");
   assert(lon >= -pi     && lon <= pi &&     "computeXYZfromLatLon: lon out of range");

   sincos(lat, &sin_lat, &cos_lat);
   sincos(lon, &sin_lon, &cos_lon);

   const double x = rearth_eq*cos_lat*cos_lon;
   const double y = rearth_eq*cos_lat*sin_lon;
   const double z = rearth_pole*sin_lat;

   return make_double3(x, y, z);
}

__inline__ __device__ void Error(const char* msg) {
   printf("Error: %s\n", msg);
   assert(false);
}

struct Tracers {
    double tracer[N_TRACERS];
};

// ------------------------------LocalSigma------------------------------

template<int IVCOR>
struct LocalSigma {};   

template<>
struct LocalSigma<2> {
   __device__ double operator[](int k) const { return 0.0; }
   __device__ double& operator[](int k) { assert(false && "LocalSigma<2>: operator[] not implemented"); return (double &)*((double *)nullptr); }
};

template<>
struct LocalSigma<1> {
   double localSigma[nvrt];
   __device__ double operator[](int k) const { return localSigma[k]; }
   __device__ double& operator[](int k) { return localSigma[k]; }
};


// ----------------------------- Frames -----------------------------

// Frames are used to transform vectors between different coordinate systems.
// They are templated in such a way that when ics==1, they have (most likely) no overhead.
template<int ICS>
struct Frame {};

template<>
struct Frame<1> {
   __forceinline__ __device__ void update(double, double) {}
   __forceinline__ __device__ Vector3 G2L(Vector3 vi) const { return vi; }
   __forceinline__ __device__ Vector3 L2G(Vector3 vi) const { return vi; }
   __forceinline__ __device__ Vector2 projHorVec(Vector2 v, const Frame<1>&) const { return v; }
   __forceinline__ __device__ Vector3 nUp(void)         const { return {0.0, 0.0, 0.0}; }
   __forceinline__ __device__ Vector3 nZonal(void)      const { return {0.0, 0.0, 0.0}; }
   __forceinline__ __device__ Vector3 nMeridional(void) const { return {0.0, 0.0, 0.0}; }
};

template<>
struct Frame<2> {
   // Column 0 (East):  (-sin_lon, cos_lon, 0)              - store 2 values, t20=0 implicit
   // Column 1 (North): Normalized oblate ellipsoid vector  - store 3 values
   // Column 2 (Up):    cross(col0, col1)                   - recomputed on demand
    double t00, t10;        // col0 (t20=0 implicit)
    double t01, t11, t21;   // col1

   __device__ void update(double lat, double lon) { // XXX use lat, lon instead of ar1, ar2?

      assert(lat >= -pi/2.0 && lat <= pi/2.0 && "GRID: lat out of range");
      assert(lon >= -pi     && lon <= pi     && "GRID: lon out of range");

      double sin_lon, cos_lon, sin_lat, cos_lat;
      sincos(lon, &sin_lon, &cos_lon); // These functions are faster than sin(ar1) and cos(ar1).
      sincos(lat, &sin_lat, &cos_lat);

      t00 = -sin_lon;
      t10 =  cos_lon;
      t01 = -cos_lon*sin_lat;
      t11 = -sin_lon*sin_lat;
      t21 = (rearth_pole/rearth_eq)*cos_lat;

      const double col1_norm = sqrt(t01*t01 + t11*t11 + t21*t21);
      assert (col1_norm != 0.0 && "GRID: col1_norm=0");

      const double col1_norm_inv = 1.0/col1_norm;
      t01 = t01*col1_norm_inv;
      t11 = t11*col1_norm_inv;
      t21 = t21*col1_norm_inv;

      // The 3rd column (col2) is the cross product of the first two and is computed when needed.


   }

   __forceinline__ __device__ void print(void) const {
      printf("frame = [%.16f, %.16f, %.16f, %.16f, %.16f]\n", t00, t10, t01, t11, t21);
   }

   __forceinline__ __device__ void getCol2(double &t02, double &t12, double &t22) const {
      // cross(col0, col1), with t20=0 implicitly:
      // |  i       j       k   |
      // | t00     t10      0   |
      // | t01     t11     t21  |
      t02 =  t10*t21;
      t12 = -t00*t21;
      t22 =  t00*t11 - t10*t01;
   }


   __forceinline__ __device__ Vector2 projHorVec(Vector2 const &vec0, Frame<2> const &frameout) const {
   // Transforms a horizontal vector from one frame to another and projects it to the new local horizontal plane.
   // Inputs:
   //        vec0: vel. in this frame;
   //        frameout: Output frame;
   // Output: 
   //        vel. in new frame  
   //
   // Projects a horizontal vector from frame 1 (this) to frame 2 (frameout),
   // discarding the z-component of the result to keep only the local horizontal
   // plane of frame 2.
   //
   // Since F1 and F2 are global-to-local rotation matrices, transforming from
   // frame 1 to frame 2 requires F2 * F1^{-1} = F2 * F1^T, as rotation matrices
   // are orthogonal. This product is computed as dot products between columns of
   // F1 and columns of F2, which is equivalent to multiplying rows of F1^T by
   // columns of F2.// 

      auto const fo = frameout; // Just an alias for the output frame.
      Vector2 res;

      // Since frame tensors are orthogonal, the transpose is the inverse. That's why the indices
      // of t and fo.t are the same.
      const double d00 = t00*fo.t00 + t10*fo.t10;
      const double d10 = t01*fo.t00 + t11*fo.t10;

      const double d01 = t00*fo.t01 + t10*fo.t11;
      const double d11 = t01*fo.t01 + t11*fo.t11 + t21*fo.t21;

      res.x = vec0.x * d00 + vec0.y * d10;
      res.y = vec0.x * d01 + vec0.y * d11;

      return res;
   }

   __forceinline__ __device__ Vector3 G2L(Vector3 const &vi) const {
      Vector3 vo;
      double t02, t12, t22;

      // Compute the 3rd column of the frame tensor
      getCol2(t02, t12, t22);

      vo.x = vi.x*t00 + vi.y*t10;
      vo.y = vi.x*t01 + vi.y*t11 + vi.z*t21;
      vo.z = vi.x*t02 + vi.y*t12 + vi.z*t22;

      return vo;
   
   }

   __forceinline__ __device__ Vector3 L2G(Vector3 const &vi) const {
      Vector3 vo;
      double t02, t12, t22;

   // Compute the 3rd column of the frame tensor
      getCol2(t02, t12, t22);

      // F * vi (non-transposed, so we use rows of F)
      vo.x = vi.x*t00 + vi.y*t01 + vi.z*t02;
      vo.y = vi.x*t10 + vi.y*t11 + vi.z*t12;
      vo.z =            vi.y*t21 + vi.z*t22;  // t20=0 term dropped

      return vo;
   }

   // Unit vectors that define the local frame
   __forceinline__ __device__ Vector3 nUp(void) const {
      // Column 2 = cross(column0, column1)
      return { t10*t21, -t00*t21, t00*t11 - t10*t01 };
   }

   __forceinline__ __device__ Vector3 nZonal(void) const {
      // Column 0
      return { t00, t10, 0.0 };
   }

   __forceinline__ __device__ Vector3 nMeridional(void) const {
      // Column 1
      return { t01, t11, t21 };
   }

};


// Forward declarations
struct Side;
struct Element;
struct Node;
struct Mesh;
struct NodeBoundary;

void showPoints(Point3 x1, Point3 x2, int jlev, const Side *origSide, unsigned int color = 0x0000FFu, char marker = 'o', int width = 1, int dashLen = 0, int label = -1);

__device__ __forceinline__ bool celestialWalk(Element *&curElem, Point3 x0, Point3 &xt, const Side *origSide, Side *& lastCrossSide, int jlev);


struct CubicSpline { 
   double *x   = nullptr;
   double *y   = nullptr;
   double *yp  = nullptr;
   double *ypp = nullptr;
};



// ------------------------------------------------------------
// Data structures for the system matrix and vector
// ------------------------------------------------------------

// Since every device owns a partition, we can divide a matrix A into row blocks A_p, where p is the device/partition index.
// Here we are abusing notation assuming that the global indices are consecutive per partition.
// In general, partition p will contain indices that are not consecutive in the global index space.

__device__ __forceinline__ int localIdx(Node* node);

// Wrapper of CuPyx CSR matrix. Used to store individual A_p rowblocks of the matrix.
struct CsrBlockMatrix {
   double* data;
   int*    indices;
   int*    indptr;
   int nrows;
   int ncols;

   __device__ CsrBlockMatrix(double* data, int* indices, int* indptr, int nrows, int ncols) : data(data), indices(indices), indptr(indptr), nrows(nrows), ncols(ncols) {}

   __device__ __forceinline__ void clearRow(Node* node) { 
      const int i = localIdx(node);
      for (int k = indptr[i]; k < indptr[i + 1]; ++k) data[k] = 0.0;
    }

   __device__ __forceinline__ void add(int i, int j, double val) {
      const int start = indptr[i];
      const int end   = indptr[i + 1];
      for (int k = start; k < end; ++k) {
         if (indices[k] == j) {
            data[k] += val;
            return;
         }
      }
      assert(false && "CSR entry (i,j) missing");
   };
};


// A block row matrix. All coefficients stored in a single device.
struct RowBlockMatrix {
   CsrBlockMatrix Ap[nPartitions];

   __device__ __forceinline__ RowBlockMatrix(CsrBlockMatrix const* src) {
      for (int q = 0; q < nPartitions; ++q)
         this->Ap[q] = src[q];
   }

   __device__ __forceinline__ void clearRow(Node* n) {
      for (int q = 0; q < nPartitions; ++q)
         Ap[q].clearRow(n);
   }

   __device__ __forceinline__ void add(Node* row, Node* col, double val)
   {
         assert(row->part == mesh->partNum);

         const int i = localIdx(row);
         const int j = localIdx(col);
         const int p = row->part; // Unused. Included for symmetry. Doesn't generate any code.
         const int q = col->part;

         Ap[q].add(i, j, val);

   }

};


// A full vector. Coefficients can be stored in different devices.
struct FullVector {
   double* data[nPartitions];   // Pointer to vectors in different devices

   __device__ FullVector(double* const src[nPartitions]) {
      for (int p = 0; p < nPartitions; ++p)
         data[p] = src[p];
   }

   __device__ __forceinline__ double& operator()(Node* node) {
      const int part = node->part;
      const int locIdx = localIdx(node); // Local index of the node in the partition
      return data[part][locIdx];
   }

   __device__ __forceinline__ double  operator()(Node* node) const {
      const int part = node->part;
      const int localIdx = node - mesh->nodes[part]; // Local index of the node in the partition
      return data[part][localIdx];
   }
};


__global__ void setupSolverMatricesAndVectors(
   double** data_ptrs,
   int**    indices_ptrs,
   int**    indptr_ptrs,
   int*     nrows,
   int*     ncols,
   double** x_ptrs,
   double** b_ptrs);


// ------------------------------------------------------------


struct Element {

   // Geometry information
   int idx;                   // Index of the element
   int16_t i34;               // Number of sides for the element 
   Side *sides[4];            // Sides of the element
   Node *nodes[4];            // Nodes of the element
   Element *nghbElems[4];     // Neighbors of the element

   int16_t kb;                // Bottom index of the element
   double area;               // Area of the element
   double ssign[4];           // Sign associated with each side of an element
   double radius;             // Element equivalent radii

   Frame<ics> frame;

   Vector2 gradPhi[4];        // Derivatives of shape functions
   Vector2 gradEta;           // Gradient of free surface elevation
   Vector2 gradPr;            // Gradient of pressure    
   Vector2 gradETP;           // Gradient of Earth Tidal Potential    XXXXX Really used?

   bool isDry;                // True if the element is dry

   Point2 center;            // Element center coordinates (global coordinates)

   double depth;              // Depth at element (min of all nodes)
   double z[nvrt];            // z-coord  (local frame - vertical up)

   // Physical information
   Tracers tracers[nvrt];     // Tracer values at each vertical level  

   static const bool isblock = false;              // Always false

   double w[nvrt];            // Vertical velocity at element centers & whole levels, calculated using F.V.M. 
                              // For hydrostatic model, this is the same as we(); for non-hydrostatic model, this is only used in upwind transport

   Vector2 ghat1;             // Element averaged ghat (\breve{G} in notes)
   

   double vsource;            // Volumetric source/sink term [m^3/s]

   double rho[nvrt];          // Density at each prism center
   double rho_mean[nvrt];     // Mean density at each prism center
   CubicSpline rho_spl;
   double rho_zz[nvrt];       // Second derivative of density at each prism center
   double rho_z[nvrt];        // First derivative of density at each prism center
   Vector2 gradRho[nvrt];     // Gradient of de-meaned density at each prism center
   

   __inline__ __device__ Point2 outputCoords(void) const;

   __inline__ __device__ double demeanedDensity(int k) const { return rho[k] - rho_mean[k]; }  

   __device__ __forceinline__ Point3 toLocalFrame(Point3 const& p) const {
   // If ics == 2, convert the global point to the local frame of the element
      
      if constexpr (ics == 1) {
         return p;
      } else {
         const Point3 c = {center.x, center.y, centerZ()};
         return frame.G2L(p - c);
      }
   }

   __device__ __forceinline__ Point3 toGlobal(Point3 const& p) const {
   // If ics == 2, convert the local point to a global point
      
      if constexpr (ics == 1) {
         return p;
      } else {
         const Point3 c = {center.x, center.y, centerZ()};
         return c + frame.L2G(p);
      }
   }

   __inline__ __device__ bool isTri (void) const { return i34 == 3; }
   __inline__ __device__ bool isQuad(void) const { return i34 == 4; }

   __inline__ __device__ bool is2D(void) const { return kb == nvrt - 2; }
   __inline__ __device__ bool is3D(void) const { return kb <  nvrt - 2; }

   __inline__ __device__ int nxq(int index) const { assert(index + 16*i34 >= 0); return (index + 16*i34) % i34; }
   __inline__ __device__ int nodeLocalIdx(Node *node) const {
      // Find the local index of the node in the element (0,1,2,3)
      int res = __ffs((nodes[0] == node) | 
                     ((nodes[1] == node) << 1) | 
                     ((nodes[2] == node) << 2) | 
                     ((nodes[3] == node) << 3)) - 1;
      
      assert(res >= 0 && res < i34 && nodes[res] == node && "nodeLocalIdx: node not found in element");
      return res;  
   }

   __inline__ __device__ int sideLocalIdx(const Side *side) const {
      // Find the local index of the side in the element (0,1,2,3)
      int res = __ffs((sides[0] == side) | 
                     ((sides[1] == side) << 1) | 
                     ((sides[2] == side) << 2) | 
                     ((sides[3] == side) << 3)) - 1;
      
      assert(res >= 0 && res < i34 && sides[res] == side && "sideLocalIdx: side not found in element");
      return res;
   }

   __inline__ __device__ int findLevel(double z1, int k) const {
      // Starting at index k, find sequentially the level that contains z

      // For testing purposes, forgets k and performs a full search.
      for (k = kb; k < nvrt - 1; k++) {
         if (z1 >= z[k] && z1 < z[k + 1]) return k;
      }
      return nvrt - 2;


      // if z is already betweeen z[k] and z[k+1], return k
      if (z1 >= z[k] && z1 < z[k + 1]) return k;
      else if (z1 <  z[k]) { 
         // Search for increasing levels
         for (int kk = k; kk < nvrt - 1; kk++) {
            if (z1 >= z[kk] && z1 < z[kk + 1]) return kk;
         }
         return nvrt - 1;
      }
      else {
         // Search for decreasing levels
         for (int kk = k; kk >= kb; kk--) {
            if (z1 >= z[kk] && z1 < z[kk + 1]) return kk;
         }
         return 0;
      }
   }

   __inline__ __device__ double z_prism(int k) const { assert(k < nvrt - 1); return 0.5*(z[k+1] + z[k]); }
   __inline__ __device__ double Dz(int k)      const { assert(k < nvrt - 1); return z[k+1] - z[k]; }


   __device__ __forceinline__ double centerZ(void) const;

   // "Local" horizontal coordinates at the nodes. They are computed from the 3D global node coordinates. For ics == 1 and 2
   // Local is in quotes because for ics == 1, the coordinates are unchanged. It just means that in these coordinates, the element is horizontal.
   __device__ __forceinline__ Point2 getLocalHorizontalNodePosition(int idxNode) const;
   __device__ __forceinline__ void getLocalHorizontalNodePositionsTri (Point2 &x0, Point2 &x1, Point2 &x2, 
                                                                       const Vector3& c, const Frame<ics>& dstFrame,int idx = 0) const;
   __device__ __forceinline__ void getLocalHorizontalNodePositionsTri (Point2 &x0, Point2 &x1, Point2 &x2, int idx = 0) const;
   __device__ __forceinline__ void getLocalHorizontalNodePositionsQuad(Point2 &x0, Point2 &x1, Point2 &x2, Point2 &x3,
                                                                       const Vector3& c, const Frame<ics>& dstFrame) const;
   __device__ __forceinline__ void getLocalHorizontalNodePositionsQuad(Point2 &x0, Point2 &x1, Point2 &x2, Point2 &x3) const;


   // (eta,xi) coordinates of the Gauss-Legendre quadrature points (idxNode = 0,1,2,3)
   __device__ __forceinline__ Point2 GaussPoint(int idxNode) const {
       constexpr double g = 0.5773502691896257; // 1/sqrt(3)
       assert(idxNode >= 0 && idxNode < 4 && "GaussPoint: idxNode out of range");
       const double eta = (idxNode & 1) ? +g : -g;
       const double xi  = (idxNode & 2) ? +g : -g;
       return {eta, xi};
   }  



   // Shape function Phi at a given point p = (eta,xi), idxNode = node index 0,1,2,3
   __device__ __forceinline__ double Phi_quad(int idxNode, Point2 p) const {
      return 0.25*(1.0 + p.x*eta_n(idxNode))*(1.0 + p.y*xi_n(idxNode));
   }

   // Derivative of shape function Phi_idx/dxi at a given point p = (eta,xi)
   __device__ __forceinline__ double dPhi_dxi_quad(int idxNode, Point2 p) const {
      return 0.25*xi_n(idxNode)*(1.0 + p.x*eta_n(idxNode));
   }   

   // Derivative of shape function Phi_idx/deta at a given point p = (eta,xi) 
   __device__ __forceinline__ double dPhi_deta_quad(int idxNode, Point2 p) const {
      return 0.25*(1.0 + p.y*xi_n(idxNode))*eta_n(idxNode);
   }      
   
   __inline__ __device__ double detJ_eta_quad(Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const;

   __inline__ __device__ double detJ_xi_quad(Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const;

   __inline__ __device__ double detJ_quad(Point2 p, Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const 
         { return (1.0/4.0)*area + detJ_eta_quad(x0, x1, x2, x3)*p.x+ detJ_xi_quad(x0, x1, x2, x3)*p.y; }

   // Check if the element is dry using the more up-to-date elevation.
   __inline__ __device__ bool isDry2(void) const;

   __inline__ __device__ bool isEdgeBoundary(Side *side) const;

   __forceinline__ __device__ bool inverseBilinearQuad(Point2 p, double &xi, double &eta, double shp[4]) const;

   __inline__ __device__ double computeSignedArea(void) const;


   enum EFixIf { fiNo, fiOnlyInside, fiAnywhere };

   __device__ bool computeShapeCoords(Point2 &xt, double shpCoords[4], EFixIf fixNegShapeCoords = fiNo, const Side *origSide = nullptr) const;
   
   __device__ __forceinline__ void shapeCoordsToXY(double shpCoords[4], Vector2 &xt) const;


   __device__ __forceinline__ double intPhiPhiTri(int i, int j) const;

   __device__ __forceinline__ double intGradPhiGradPhiTri(int i, int j) const;

   __device__ __forceinline__ double intPhiPhiQuad(int i, int j, Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const;

   __device__ __forceinline__ double intGradPhiGradPhiQuad(int i1, int i2, Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const;

   __device__ __forceinline__ Vector2 grad2(double *a) const
   {
      double gx = a[0]*gradPhi[0].x + a[1]*gradPhi[1].x + a[2]*gradPhi[2].x;
      double gy = a[0]*gradPhi[0].y + a[1]*gradPhi[1].y + a[2]*gradPhi[2].y;
      if (isQuad()) {
         gx += a[3]*gradPhi[3].x;
         gy += a[3]*gradPhi[3].y;
      }
      return {gx, gy};
   }

   __device__ __forceinline__ Vector2 grad2(double Node::* field) const
   {
      double gx = (nodes[0]->*field)*gradPhi[0].x + (nodes[1]->*field)*gradPhi[1].x + (nodes[2]->*field)*gradPhi[2].x;
      double gy = (nodes[0]->*field)*gradPhi[0].y + (nodes[1]->*field)*gradPhi[1].y + (nodes[2]->*field)*gradPhi[2].y;
      if (isQuad()) {
         gx += (nodes[3]->*field)*gradPhi[3].x;
         gy += (nodes[3]->*field)*gradPhi[3].y;
      }
      return {gx, gy};
   }

   __device__ __forceinline__ Vector3 n(int k, double &area3D) const;
   __device__ __forceinline__ Vector3 n(int k) const {
      double area3D;
      return n(k, area3D);
    }



   __device__ void setFrame(void) {
      if constexpr(ics == 1) return;

      double lat, lon;
      const Point3 p = toPoint3(center, centerZ());

      // Compute the latitude and longitude of the element center
      computeLatLon(p, lat, lon);

      frame.update(lat, lon);        
   }


   __device__ double interpolate(const double shapeCoords[4], const double zAtNodes[4]) const
   {
      // Interpolate the value at the point (xt,yt) using the shape coordinates and the values at the nodes
      // This function works for both triangular and quadrilateral elements, since shapeCoords[3] is always 0.0 for triangular elements
      return shapeCoords[0]*zAtNodes[0] + shapeCoords[1]*zAtNodes[1] + shapeCoords[2]*zAtNodes[2] + shapeCoords[3]*zAtNodes[3];
      
   }

   __device__ Vector3 interpolateVelocities(Point3 xt, const Side *origSide = nullptr) const;

   // ------------------------------------------------------------------------------------------------
   // These are functions that work on all nodes/sides of the element.
   // They use a standard, but not well known, way of accessing the fields of the nodes.
   //
   // An example of how to use these functions is:
   //   
   //    getNodalValues(L, &Node::compZCoord, zLAtNodes);
   // 
   // Where L is the level at which the field is to be interpolated. The field can be either
   // a double (depth), an array (u[nvrt]), or a function that depends on the level (Node::compZCoord(int k)). 
   // There are different versions of getNodalValues for each of these cases.

   template<typename T>
   __device__ __forceinline__ void getNodalValues(T Node::* field, T res[4]) const;

   template<typename T>
   __device__ __forceinline__ void getNodalValues(int level, T (Node::* field)[nvrt], T res[4]) const;
  
   template<typename T>
   __device__ __forceinline__ void getNodalValues(int level, T (Node::*func)(int) const, T res[4]) const;

   template<typename T>
   __device__ __forceinline__ void getSideValues(T Side::* field, T res[4]) const;

   template<typename T>
   __device__ __forceinline__ void getSideValues(int level, T (Side::* field)[nvrt], T res[4]) const;

   template<typename T>
   __device__ __forceinline__ void getSideValues(int level, T (Side::*func)(int) const, T res[4]) const;

   
   __forceinline__ __device__ void projectNodeVectorsToElementFrame(Vector2 *v) const;
   __forceinline__ __device__ void projectSideVectorsToElementFrame(Vector2 *v) const;


   template<typename T>
   __device__ __forceinline__ T interpolate(int level, const double shapeCoords[4], T (Node::*func)(int) const) const;

   template<typename T>
   __device__ __forceinline__ T interpolate(const double shapeCoords[4], T Node::* field) const;

   template<typename T>
   __device__ __forceinline__ T interpolate(int level, const double shapeCoords[4], T (Node::* field)[nvrt]) const;

   template<typename T>
   __device__ __forceinline__ T mean(T Node::* field) const;

   template<typename T>
   __device__ __forceinline__ T mean(T Side::* field) const;

   template<typename T>
   __device__ __forceinline__ T mean(int level,  T (Side::*func)(int) const) const;

   template<typename T>
   __device__ __forceinline__ T mean(int level, T (Node::* field)[nvrt]) const;

   template<typename T>
   __device__ __forceinline__ T mean(int level, T (Side::* field)[nvrt]) const;



   // __device__ double mean(double (Node::*func)(int) const) const {
   //    double values[4];
   //    getNodalValues(level, func, values);
   //    return shapeCoords[0]*values[0] + shapeCoords[1]*values[1] + shapeCoords[2]*values[2] + shapeCoords[3]*values[3];
   // }
   // ------------------------------------------------------------------------------------------------




   __device__ __forceinline__ Point2 getLocalNodePosition(int idx, const Side* origSide) const;

   __device__ __forceinline__ Node *getNodeInOrder(int idx) const;

   __device__ __forceinline__ void computePredicates(const Point2 x0, const Point2& xt, double &fac, const Side* origSide,
                                                    Bool4& isLeft, Bool4& isObtuse, Bool4& isBisector, int jlev) const;
    

};


struct BoundarySegment {

   bool isLand;

   double clen;  // Length of the boundary segment
   double carea; // Area of the boundary segment

   int ifltype;  // Type of flow boundary condition
   
   // Elevation boundary condition
   int iettype;  // Type of elevation boundary condition
   double eth;   // Elevation time history

   // Flow boundary condition
   double qthcon;  // Flow time history constant
   double vobc1;   // Velocity coefficient 1
   double vobc2;   // Velocity coefficient 2

   // Tracer boundary conditions
   // double trobc[natrm];    // Tracer boundary condition nudging factor
   // int  itrtype[natrm];    // Tracer boundary condition type

   BoundarySegment(bool _isLand) {
      isLand = _isLand;
      clen = 0.0;
      carea = 0.0;
      ifltype = 0;
      iettype = 0;
      eth = 0.0;
   }

   BoundarySegment(void) {
      isLand = false;
      clen = 0.0; 
      carea = 0.0;
      ifltype = 0;
      iettype = 0;
      eth = 0.0;
   }

   
};


struct SideBoundary {
   Side *side; // Pointer to the side that this boundary segment belongs to
   BoundarySegment *segment; // Pointer to the boundary segment

   Side *leftSide, *rightSide; // Pointers to the left and right sides of the boundary segment

   Vector2 ubc[nvrt]; // velocity boundary condition
   
   __inline__ __device__ double getBoundaryArrayValueWeight(void) const;

   __inline__ __device__ double etaMeanBnd(void) const;

   __inline__ __device__ Vector2 umo (int jfr) const;
   __inline__ __device__ Vector2 ufa (int jfr) const;  
   __inline__ __device__ Vector2 uth (int k) const;
   __inline__ __device__ double  nuth(int k) const;
   __inline__ __device__ double  eth (int k) const;
};


struct Side {

   // Geometry information
   int idx;               // Index of the side
   Node *n1, *n2;
   Element *elems[2];     // Elements that share this side XXX (why not e1,e2?)

   double len;            // Length of the side. Used to be called distj; but I think this name is less clear now.
   Vector2 center;        // Location of side center. 
   double zc[nvrt];       // Side center z-coord at each vertical level

   Frame<ics> frame;
   
   SideBoundary *sideBnd; // Pointer to additional boundary information. For non-boundary sides, 
                          // this is nullptr
   
   int16_t kb;            // Bottom index of the side
   bool isDry;            // True if the side is dry
 
   Vector2 n;             // Side unit normal vector (for ics=2: these are done inside local lon/lat frame)

   Tracers tracers[nvrt]; // Tracer concentrations at nodes

   // Physical information
   Vector2 u2[nvrt];            // Unit normal vector (for ics=2: these are done inside local lon/lat frame)
   Vector2 btu[nvrt];           // Back-tracked velocity
   Vector2 U;                   // Vertically averaged velocity vector  

   Vector2 nuGradU[nvrt];          // nu*grad u

   Vector2 bcc[nvrt];           // Baroclinic term
   Vector2 bcc_TDz[nvrt];       // Baroclinic term at previous timestep

   Vector2 gradEta1;             // Gradient of elevation at previous timestep
   Vector2 gradEta2;             // Gradient of elevation at current timestep
   Vector2 gradPr;               // Gradient of pressure
   Vector2 gradETP;              // Gradient of Earth Tidal Potential


   double hhat;                 // \breve{H} in notes

   double shapiro;

   double flux_adv_hface[nvrt]; // Flux of horizontal advection at each vertical level

   double coriolis;             // Coriolis parameter

   // Vegetation parameters
   double veg_c2;
   double veg_beta;

   // Bottom drag parameters.
   double chi;
   double chi2;

   __inline__ __device__ Point2 outputCoords(void) const;

   __inline__ __device__ double x(double t) const;
   __inline__ __device__ double y(double t) const;
   __inline__ __device__ double z(double t) const;
   __inline__ __device__ double lon(void) const;
   __inline__ __device__ double lat(void) const;
   __inline__ __device__ Point2 p2_t(double t) const;

   __device__ void setFrame(void) {
      if constexpr(ics == 1) return;

      double lat, lon;
      const Point3 p = toPoint3(center, centerZ());

      // Compute the latitude and longitude of the element center
      computeLatLon(p, lat, lon);

      frame.update(lat, lon);        
   }
   

   __device__ __forceinline__ Point3 toLocalFrame(Point3 const& p) const {
   // If ics == 2, convert the global point to the local frame of the side
      
      if constexpr (ics == 1) {
         return p;
      } else {
         const Point3 c = {center.x, center.y, centerZ()};
         const Point3 local = frame.G2L(p - c);
         return local;
      }
   }

   __device__ __forceinline__ Point3 toGlobal(Point3 const& p) const {
      if constexpr (ics == 1) {
         return p;
      } else {
         const Point3 c = {center.x, center.y, centerZ()};
         return c + frame.L2G(p);
      }
   }


   __device__ __forceinline__ Vector2 projectToLocalFrame(Vector2 const& v, Frame<ics> const& v_frame) const {
      if constexpr(ics == 1) return v;
      else return v_frame.projHorVec(v, frame);
   }

   __inline__ __device__ double centerZ(void)      const;
   __inline__ __device__ Point3 center3D(int k)    const { assert(false); assert(k < nvrt); return {center.x, center.y, zc[k]}; }  // USED??? CONSISTENT?? XXXXX
   __inline__ __device__ double Dz(int k)          const { assert(k < nvrt); return zc[k] - zc[k - 1]; }
   __inline__ __device__ double area(int k)        const { assert(k < nvrt - 1); return len*Dz(k + 1); }
   __inline__ __device__ double areaClamped(int k) const { assert(k < nvrt - 1); return len*(zc[max(k + 1, kb)] - zc[max(k, kb)]); }

   __inline__ __device__ Side* neigh(int index) const { 
      // Each internal side has 4 neighbors, they are the closest sides to each of the two elements that share this side
      assert(index >= 0 && index < 4); 
      assert(elems[0] != nullptr && elems[1] != nullptr);

      const int ielem = index & 1;           // 0 or 1
      const int iside = 2*(index >> 1) - 1;  // -1 or 1

      const Element* elem = elems[ielem];

      const int locSideIdx = elem->sideLocalIdx(this); // Local index of this side in the element

      return elem->sides[elem->nxq(locSideIdx + iside)];

   }

   __device__ Vector2 grad2(double Node::* field, bool (*condition)(const Element*)) const;

   __forceinline__ __device__ Vector2 t(void) const { return {-n.y, n.x}; } // Tangent vector

   __inline__ __device__ bool isBnd(void) const { return sideBnd != nullptr; } // True if the side is a boundary side
   __inline__ __device__ bool isOpenBnd(void) const { return isBnd() && sideBnd->segment->ifltype == 0; } // True if the side is an open boundary side
   __inline__ __device__ bool isInternal(void) const { return elems[1] != nullptr; } // Internal side
   __inline__ __device__ int iadv(void) const;

   __inline__ __device__ double grav3(void) const;
   
   __inline__ __device__ double eta2(void) const;
   __inline__ __device__ double eta1(void) const;
   __inline__ __device__ double w2(int k) const;
   __inline__ __device__ Vector2 u2Bottom(int dummy = 0) const { return u2[kb]; }  // XXXX dummy argument to satisfy the template in mean()
   __inline__ __device__ double depth(void) const;
   __inline__ __device__ double htot2(void) const { return eta2() + depth(); }
   __inline__ __device__ double htot1(void) const { return zc[nvrt - 1] - zc[kb]; }

   __inline__ __device__ double dfv(int k) const;

   __inline__ __device__ double veg_h(void)  const;
   __inline__ __device__ double veg_nv(void) const;
   __inline__ __device__ double veg_alpha_vert_mean(void) const;
   __inline__ __device__ double veg_alpha3D(int k) const;
   __inline__ __device__ double veg_c(void)  const { return veg_c2/(1.0 + veg_c2); }

   __inline__ __device__ double Cd(void) const;

   __inline__ __device__ Vector2 tau(void) const;

   __inline__ __device__ double zb(void) const { return zc[kb]; }
   

   __inline__ __device__ double zTopCanopy(void) const { return zc[kb] + veg_h(); } 
   __inline__ __device__ double zLevelVegRatio(int k) const { 
      // Ratio ocupied by the vegetation at vertical level k. It is 1.0 if the canopy is fully above the level k and 0.0 if it is fully below the level k
      return max(0.0, min(1.0, (zTopCanopy() - zc[k-1])/Dz(k))); 
   } 
   __inline__ __device__ bool hasVegetation(int k) const { return zTopCanopy() > zc[k]; }
   __inline__ __device__ bool isSubmerged(void)    const { return veg_h() < 0.99*htot2(); }
   __inline__ __device__ bool is2D(void) const { assert(kb <= nvrt - 2); return kb == nvrt - 2; }
   __inline__ __device__ bool is3D(void) const { assert(kb <= nvrt - 2); return kb <  nvrt - 2; }

   // Get valid starting element (for backtracking)
   __inline__ __device__ Element* getStartElement() {
      return elems[0] ? elems[0] : elems[1];
   }

   
   __forceinline__ __device__ Element* otherElem(const Element* elem) const {
      // Given an element, returns the other element that shares this side
      return (elems[0] == elem) ? elems[1] : elems[0];
   }

   __forceinline__ __device__ Node* otherNode(Node* node) {
      // Given a node, returns the other node that shares this side
      return (n1 == node) ? n2 : n1;
   }
 
   __inline__ __device__ void btrack(int jlev, double dt, double dtbk, Vector3 &u_int) const;




};


// This struct is used to store the boundary information for a node.
struct NodeBoundary {   
   int idx;
   BoundarySegment *segment;
   Node *node;
   
   SideBoundary *sideBnd1; // Pointer to the boundary side object (can be nullptr)
   SideBoundary *sideBnd2; // Pointer to the boundary side object (can be nullptr)

   bool isLand;

   double etaic;  // Elevation at nodes used for essential boundary conditions XXXX Where is this value set???
   double elbc;   // Elevation boundary condition
   double eth;    // Elevation time history

   double eta_mean; // Mean elevation at nodes

   bool isEssentialBC;

   Vector2 uth[nvrt];  // "Time history" velocity

   // Arrays for harmonic constituents (if needed)
   double  emo[nbfr];  // Elevation amplitude
   double  efa[nbfr];  // Elevation phase
   Vector2 umo[nbfr];  // velocity amplitude
   Vector2 ufa[nbfr];  // velocity phase


   NodeBoundary(void) {
      isLand = false;
   }
   NodeBoundary(bool _isLand) {
      isLand = _isLand;
   }

   
   __forceinline__ __device__ int numSidesForAverage() const {
   // Returns the number of sides that are connected to the node and have the same segment.
   // This is used to compute the average velocity at the node for open boundary nodes.
      int res = 0;
      if (sideBnd1 != nullptr && sideBnd1->segment == segment) res++;
      if (sideBnd2 != nullptr && sideBnd2->segment == segment) res++;
printf("numSidesForAverage: res = %d, segment = %p, sideBnd1 = %p, sideBnd2 = %p\n", res, segment, sideBnd1, sideBnd2);      
      return res;
   }


   __forceinline__ __device__ SideBoundary *sideBnd(int index) { 
      assert(index == 0 || index == 1); 
      return (index == 0) ? sideBnd1 : sideBnd2; 
   }

   __forceinline__ __device__ void linkSideBoundary(SideBoundary *sb) {
      if (isLand) return;  // land never uses these links

      SideBoundary *nullSideBoundary = nullptr;
      uintptr_t null = (uintptr_t)nullSideBoundary;

      if (atomicCAS((uintptr_t *)&sideBnd1, null, (uintptr_t)sb) == null) return;
      if (atomicCAS((uintptr_t *)&sideBnd2, null, (uintptr_t)sb) == null) return;

      // both slots taken. OK if this side was already linked
      assert(sideBnd1 == sb || sideBnd2 == sb);
   }


};

struct Node {

   // Geometry information
   int idx;            // Global index of the node
   int16_t part;       // Partition of the node

   Node **nodes;       // Array of pointers to the nodes that are neighbors of this node (there is no guarantee that they are connected)
   int16_t numNodes;   // Number of nodes neighbors of this node

   Element **elements; // Array of pointers to the elements that share this node
   int16_t numElements;// Number of elements that share this node

   int16_t numSides;   // Number of sides that are connected to this node


   Point2 p;           // Horizontal coordinates at the node. (if ics == 2, this are the x,y coordinates of the 3d global cartesian coordinates. Together with gz, they form the 3D global cartesian coordinates.)
   double gz;          // Global vertical coordinate at the node.
   Vector2 u2[nvrt];
   double  w2[nvrt];
   Vector2 utmp[nvrt];  // Temporary array to store velocities
   double z[nvrt];      // Local vertical coordinate at each vertical level. 
   double lon, lat;
   double depth;

   Frame<ics> frame;

   LocalSigma<ivcor> localSigma;  

   int16_t kb;     // Bottom index of the node

   bool isDry;     // True if the node is dry

   int16_t iadv;   // Advection flag

   NodeBoundary *nodeBnd; // Pointer to additional boundary information. For non-boundary nodes, 
                          // this is nullptr

   TempVar tmp; // Temporary variable to store different types of variables

   // Physical information
   Tracers tr0[nvrt];       // Initial tracer conc. at nodes
   Tracers tracers[nvrt];   // Tracer concentrations at nodes

   double rho[nvrt];     // Density
   double dfv[nvrt];     // Vertical viscosity

   double eta1; // Elevation at nodes at previous timestep
   double eta2; // Elevation at nodes at current timestep

   double ETP;  // Earth Tidal Potential

   double pr;    // Air pressure
   double grav2; // Gravity 
   
   Vector2 dav;  // Depth averaged velocity

   Vector2 tau; // Wind stress vector

   // Vegetation parameters
   double veg_h;
   double veg_nv;
   double veg_di;
   double veg_Cd;
   double veg_h_unbent;
   double veg_di_unbent;
   double veg_nv_unbent;
   double veg_alpha_vert_mean;
   double veg_alpha0;
   double veg_alpha3D[nvrt];

   // Bottom drag parameters
   double Cd;       // Bottom drag coefficient   

   // These three are mutually exclusive
   union {
      double rmanning; // Bottom Manning's coefficient
      double rough;    // Bottom roughness length
      double drag;     // Bottom drag coefficient
   };
   Vector2 tau_bot; // Bottom stress vector

   // XXXX REMOVE FOR IMM>0???
   double bdef1 = 0.0; // Bottom deformation at previous timestep
   double bdef2 = 0.0; // Bottom deformation at current timestep

   
   double qel; // RHS
   double *coeffDiag;    // Pointer to the diagonal coefficient of the matrix
   double *coeffRow;     // Pointer to the coefficients of the row of the CRS matrix for this node.
   int16_t *csrSlot;     // csrSlot[i], is the slot/index of nodes[i] entry the CRS matrix row for this node.
                         // The index of the diagonal element is coeffDiag - coeffRow.   

   double mX; // Solution of the FEM system

   // Meteorological forcing
   Vector2 u_air;
   double p_air;  // XXXX is this the same as pr?
   double t_air;
   double q_air;
   float cRot, sRot;  // Cosine and sine of the mesh rotation matrix




   __device__ void setWind(Vector2 &u) {
   // Sets the wind vector taking into account the mesh angle
   
      const double speed = length(u);

      // Rotate wind vector to the node frame
      u_air = {u.x*cRot - u.y*sRot, u.x*sRot + u.y*cRot};

      // Sets the surface stress vector due to wind.
      if constexpr(nws == 0) {
         tau = {0.0, 0.0};

      } else if constexpr(nws == 2 && ihconsv == 1 && iwind_form == 0) {
         assert(false && "setWinds: nws==2, ihconsv==1, iwind_form==0 (tauxz/tauyz) not supported");

      } else {

         double dragcoef;

         if constexpr(iwind_form == -1) {
            const double speed_clamp = clamp(speed, 6.0, 50.0);
            dragcoef = 1.0e-3 * (0.61 + 0.063*speed_clamp);

            const double cff = 0.001293*dragcoef*speed*rampwind;
            tau = cff*u_air;
 
         } else if constexpr(iwind_form == 1) {
            if (speed <= 35.0) dragcoef = 1.00e-4 * (-0.016*speed*speed + 0.967*speed + 8.058);
            else               dragcoef = 2.23e-3 * 35.0/speed;

            dragcoef = max(dragcoef, 4.0e-4);
            const double cff = 0.001293*dragcoef*speed*rampwind;
            tau = cff*u_air;

         } else {
            tau = {0.0, 0.0};
         }
      }
   }


   __device__ void setFrame(double lat, double lon) {
      if constexpr(ics == 1) return;

      frame.update(lat, lon);        
    }

   //__device__ Point3 p3D(int k) const { return {p.x, p.y, z[k]}; }   // XXXX can this be confused with 3D as in ics=2???

   __device__ double sigmaTransitionDepth(void) const;

   __device__ double computeETP(int idxTidalPotential, double time, double* tfreq, double* tamp, double* tnf, double* tear, int* jspc) const;

   __inline__ __device__ Point2 outputCoords(void) const {
      if constexpr(ics == 1)  return p;
      else                    return {rad2deg*lon, rad2deg*lat};
   }

   __device__ __forceinline__ Point3 toLocalFrame(Point3 const& p) const {
   // If ics == 2, convert the global point to the local frame
      
      if constexpr (ics == 1) {
         return p;
      } else {
         const Point3 c = {p.x, p.y, gz};
         const Point3 local = frame.G2L(p - c);
         return local;
      }
   }

   __device__ __forceinline__ Point3 toGlobal(Point3 const& p) const {
   // If ics == 2, convert the global point to the local frame
      
      if constexpr (ics == 1) {
         return p;
      } else {
         const Point3 c = {p.x, p.y, gz};
         return c + frame.L2G(p);
      }
   }

   __inline__ __device__ double Dz(int k) const { return z[k] - z[k-1]; }
   __inline__ __device__ double htot2(void) const { return eta2 + depth; }
   __inline__ __device__ double zTopCanopy(void) const { return z[kb] + veg_h; } 
   
   __inline__ __device__ Vector3 uvw2(int k) const { return {u2[k].x, u2[k].y, w2[k]}; }

   __inline__ __device__ bool isDry2(void) const { return htot2() <= h0; } // isDry, but updated to use eta2

   __inline__ __device__ bool isBnd(void) const { return nodeBnd != nullptr; } // True if the node is a boundary node
   // __inline__ __device__ int numElements(void) const { return isNbd() ? numElements - 1 : numElements; } // Number of elements that are connected to this node

   __device__ __forceinline__ void clearMatrixRow(void) {
      // Clear the matrix row for this node

      for (int i = 0; i <= numNodes; i++) coeffRow[i] = 0.0;

   }

   __device__ void addToSparseMatrix(Node const *node, double value) {
   // Add the value to the sparse matrix for this node that corresponds to the file of node *this and the column of node *node
      
      if (node == this) {
         *coeffDiag += value;
         return;
      }

      int i;
      for (i = 0; i < numNodes; i++) {
         if (nodes[i] == node) break;
      }

      assert(i < numNodes && "addToSparseMatrix: node not found in array of nodes");
      coeffRow[csrSlot[i]] += value;
   }


   __device__ __inline__ void updateZCoordData(void);

   __device__ __forceinline__ double compZCoord(int k) const;

   __device__ void updateKb(void);

   __device__ __forceinline__ void computeBottomDragCd(void);



};




struct Mesh {
   public:  

   int partNum; // Partition number for this mesh.

   // Pointers to each GPU's local memory (size is the number of partitions).
   Node**    nodes;  
   Element** elems;  
   Side**    sides;  
   
   // Local memory
   Node*    localNodeData; 
   Element* localElemData;
   Side*    localSideData; 

   NodeBoundary *localNodeBoundaryMem; // Array of objects with the boundary information of each node
   SideBoundary *localSideBoundaryMem; // Array of objects with the boundary information of each side
   double       *coeffMatrixMem;       // Matrix coefficients for neighbors. [0..numNodes-1], doesn't include the diagonal
   BoundarySegment *localSegmentMem;   // Array of objects with the boundary segments
   
   const int npgb;  // Global number of horizontal nodes
   const int negb;  // Global number of elements
   const int nsgb;  // Global number of sides

   const int np;   // Local number of horizontal nodes 
   const int ne;   // Local number of elements
   const int ns;   // Local number of sides

   int nbp;  // Local number of boundary nodes
   int nbs;  // Local number of boundary sides

   int nEntries;  // Total number of matrix elements

   const int max_np; // Maximum number of horizontal nodes per partition (shared by all partitions)
   const int max_ne; // Maximum number of elements per partition (shared by all partitions)
   const int max_ns; // Maximum number of sides per partition (shared by all partitions)

  
   // Additional topology information
   Node    **nodeNeighNodes;     // Array of pointers to the nodes neighbors of each node
   Element **nodeNeighElems;     // Array of pointers to the elements neighbors of each node
   int16_t *nodePermutationSort; // Array of indices for the array nodes such that neighbors nodes are sorted in increasing order of their global index.

   BoundarySegment *landBoundarySegment;
   NodeBoundary *landBoundaryNode;

   // Vertical coordinate system
   // Sigma levels
   double h_s; 
   double h_c; 

   double theta_b;
   double theta_f, sinh_theta_f;
   double *sigma;  // There is a sigma_local on each node for the case ivcor=1
   double *ztot;
   double *cs;
   double *dcs;



   __device__ Mesh(int partNum,
                   int npgb, int negb, int nsgb, 
                   int np, int ne, int ns,
                   int max_np, int max_ne, int max_ns, int nseg,
                   NodeBoundary *plandNodeBndMem, BoundarySegment *plandSegMem,
                   Node **nodesPtr, Element **elemsPtr, Side **sidesPtr,
                   Node *nodesMem, Element *elemsMem, Side *sidesMem,
                   Node    **localNodeNeighNodesMem,
                   Element **localNodeNeighElemsMem,
                   int16_t *localNodePermutationSortMem,
                   double *partCoeffMatrixMem,
                   NodeBoundary *localNodeBoundaryMem,
                   SideBoundary *localSideBoundaryMem,
                   BoundarySegment *localSegmentMem,
                   int *kb, double h_c, double h_s, double theta_b, double theta_f, double *sigma, double *ztot,
                   double *cs, double *dcs
                   ) : partNum(partNum), npgb(npgb), negb(negb), nsgb(nsgb),
                       np(np), ne(ne), ns(ns),
                       max_np(max_np), max_ne(max_ne), max_ns(max_ns),
                       nodes(nodesPtr), elems(elemsPtr), sides(sidesPtr),
                       localNodeData(nodesMem), localElemData(elemsMem), localSideData(sidesMem),
                       nodeNeighNodes(localNodeNeighNodesMem), nodeNeighElems(localNodeNeighElemsMem),
                       nodePermutationSort(localNodePermutationSortMem),
                       localNodeBoundaryMem(localNodeBoundaryMem), localSideBoundaryMem(localSideBoundaryMem),
                       localSegmentMem(localSegmentMem),
                       h_c(h_c), h_s(h_s), theta_b(theta_b), theta_f(theta_f), sigma(sigma), ztot(ztot),
                       cs(cs), dcs(dcs),
                       coeffMatrixMem(partCoeffMatrixMem),  // XXX REMOVE?
                       landBoundarySegment(plandSegMem),
                       landBoundaryNode(plandNodeBndMem)
   {
      printf("Mesh constructor for partition %d\n", partNum);

      printf("npgb = %d, negb = %d, nsgb = %d, nseg = %d\n", npgb, negb, nsgb, nseg);
      printf("np = %d, ne = %d, ns = %d\n", np, ne, ns);
      printf("max_np = %d, max_ne = %d, max_ns = %d\n", max_np, max_ne, max_ns);
 
      printf("elemsPtr = %p, elems = %p\n", *elems, elemsMem);
      
      sinh_theta_f = sinh(theta_f);

      *landBoundaryNode    = NodeBoundary(true);
      *landBoundarySegment = BoundarySegment(true);   

   }

   __device__ void populateNodes(int *i34, int *isbnd, int *sharedIdxNodeBoundary,
                  double *node_coords, double *node_coords_ll, bool use_ll,
                  double *wind_rotate_angle, 
                  int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                  int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                  int16_t *nodePermutationSort,
                  int *localNodeNeighNodesStart, int *localNodeNeighElemsStart,
                  int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                  int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                  int *localNodePart,     int *localElemPart,     int *localSidePart,
                  double *drag, double *roughness, double *rmanning, int *kb, double *localSigma);
   __device__ void populateNodes2(void);
   __device__ void populateSides(int *i34, int *sharedIdxSideBoundary,
                  double *node_coords, double *node_coords_ll, bool use_ll,
                  double *wind_rotate_angle,   
                  int *isdel, int *isidenode, int *isbnd,
                  int *nond, int *nlnd, int *iond, int *ilnd,                                  
                  int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                  int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                  int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                  int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                  int *localNodePart,     int *localElemPart,     int *localSidePart,
                  double *drag, double *roughness);          

   __device__ void populateElements(int *i34, int *elnode, int *elside,
                     double *node_coords, double *node_coords_ll, bool use_ll,
                     double *wind_rotate_angle,   
                     int *isdel, int *isidenode, int *isbnd,
                     int *nond, int *nlnd, int *iond, int *ilnd,                                  
                     int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                     int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                     int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                     int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                     int *localNodePart,     int *localElemPart,     int *localSidePart,
                     double *drag, double *roughness);   
                     

   __device__ void populateBoundaryNodes(int *nond, int *iond,
                                          double *emo, double *efa, 
                                          double *umo, double *ufa, 
                                          double *vmo, double *vfa,
                                          int *ifltype, int *iettype); 

   __device__ void populateBoundarySideLinks(void);

   __device__ void populateMesh(int *i34, int *elnode, int *elside,
                   double *node_coords, double *node_coords_ll, bool use_ll,
                   double *wind_rotate_angle, int *sharedIdxNodeBoundary, int *sharedIdxSideBoundary,
                   int *isdel, int *isidenode, int *isbnd,
                   int *nond, int *nlnd, int *iond, int *ilnd,
                   int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                   int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                   int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                   int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                   int *localNodePart,     int *localElemPart,     int *localSidePart,
                   double *drag, double *roughness,
                   double *emo, double *efa, 
                   double *umo, double *ufa, 
                   double *vmo, double *vfa);


   __device__ ~Mesh(void)
   {
   }

 

}; // End of Mesh class

__device__ __inline__ double dotprod3(double* vec1, double* vec2) {
   return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}

__device__ __inline__ double dotprod4(double* vec1, double* vec2) {
   return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2] + vec1[3]*vec2[3];
}   

__device__ __inline__ double min3(double a, double b, double c) {
   return min(a, min(b, c));
}

__device__ __inline__ double max3(double a, double b, double c) {
   return max(a, max(b, c));
}

__device__ __inline__ double min4(double a, double b, double c, double d) {
   return min(a, min(b, min(c, d)));
}


// Helper functions for computing the local node index at a given offset from another.
__device__ __inline__ int nxq(int offset, int index, int size) {
   // In c++, the % operator can return negative numbers, for this reason we add 16*size 
   // to the index to make it positive
   assert(index + offset + 16*size >= 0);
   return (index + offset + 16*size) % size;  
}

__device__ __inline__ int nxq(int index, int size) {
   // In c++, the % operator can return negative numbers, for this reason we add 16*size 
   // to the index to make it positive
   assert(index + 16*size >= 0);
   return (index + 16*size) % size; 
}


extern "C" {
   __global__ void getStructSizes(int* sizes);

   __global__ void printMatrix(bool full = true);

   __global__ void returnMatrix(double *res);

   __global__ void returnX(double *res);

   __global__ void returnB(double *res);

   __global__ void computeTimeRamps(double time, int nws, bool if_source, double drampbc, double drampwind, double drampwafo, double dramp_ss, double dramp) {
      // Compute the different time ramps used in the code 
      // This function is performed by just one thread

      rampbc = 1.0;
      if constexpr(ibc == 0) {
         if (drampbc > 0.0) {
            rampbc = tanh(2.0*time/86400.0/drampbc);
         } else {
            rampbc = 1.0;
         }
      }

      if (nws != 0) {
         if (drampwind > 0.0) {
            rampwind = tanh(2.0*time/86400.0/drampwind);
         } else {
            rampwind = 1.0;
         }
      }

      if (drampwafo > 0.0) {
         rampwafo = tanh(2.0*time/86400.0/drampwafo);
      } else {
         rampwafo = 1.0;
      }

      // For source/sinks
      if (if_source && dramp_ss > 0.0) {
         ramp_ss = tanh(2.0*time/86400.0/dramp_ss);
      } else {
         ramp_ss = 1.0;
      }

      if (dramp > 0.0) {
         ramp = tanh(2.0*time/86400.0/dramp);
      } else {
         ramp = 1.0;
      }


   }
      

   __global__ void computeBedDeformation(void) {
      printf("NOT IMPLEMENTED: computeBedDeformation\n");

   }


   __global__ void setupCGSolver(double* data_d, double* b_d, int* indices, int* indptr);
   __global__ void setupMultiGPUCGSolver(double* data_d, double* b_d, int* indices, int* indptr);


}


// END INCLUDE: meshdata.h
#include <cuda_runtime.h>


__device__ double *debugPlot = nullptr;
__device__ int *debugPlotIdx = nullptr;
__device__ int debugPlotMutex = 0;



#define IDX(ptr) (ptr == nullptr ? -1 : (ptr->idx))

// Each device has a global mesh object
__device__ Mesh* mesh = nullptr; 



// Convenience macros for accessing the objects in kernels.

#define SIDE(side) const int idx_side = threadIdx.x + blockIdx.x * blockDim.x; \
                   Side *side = &(mesh->localSideData[idx_side]); \
                   if (idx_side >= mesh->ns) return;

#define ELEM(elem) const int idx_elem = threadIdx.x + blockIdx.x * blockDim.x; \
                   Element *elem = &(mesh->localElemData[idx_elem]); \
                   if (idx_elem >= mesh->ne) return;

#define NODE(node) const int idx_node = threadIdx.x + blockIdx.x * blockDim.x; \
                   Node *node = &(mesh->localNodeData[idx_node]); \
                   if (idx_node >= mesh->np) return;

#define NODE_ALL(node, valid) const int idx_node = threadIdx.x + blockIdx.x * blockDim.x;       \
                              Node *node = &(mesh->localNodeData[idx_node]); \
                              const bool valid = (idx_node < mesh->np);                   


#define NODE_BND(node_bnd) const int idx_node_bnd = threadIdx.x + blockIdx.x * blockDim.x; \
                           NodeBoundary *node_bnd = &(mesh->localNodeBoundaryMem[idx_node_bnd]); \
                           if (idx_node_bnd >= mesh->nbp) return;

#define SIDE_BND(side_bnd) const int idx_side_bnd = threadIdx.x + blockIdx.x * blockDim.x; \
                           SideBoundary *side_bnd = &(mesh->localSideBoundaryMem[idx_side_bnd]); \
                           if (idx_side_bnd >= mesh->nbs) return;  


__device__ inline void atomicAddInt16Unsafe(int16_t* address, int16_t val)
{
    uintptr_t p = (uintptr_t)address;
    assert((p & 1) == 0);

    unsigned int* base = (unsigned int*)(p & ~3);
    unsigned int addend = (p & 2) ? ((unsigned int)val << 16)
                                  :  (unsigned int)val;

    atomicAdd(base, addend);
}

Polynomial<4> *polyDischarge; // XXX read from file


//===============================================================================


/**
 * Find the element that contains the target point.
 * 
 * @param x0,y0,z0     Starting position coordinates
 * @param xt,yt,zt     Target position coordinates
 * @param dtb          Time step for tracking
 * @param u_int        Initial velocity components
 * @param nextElem     Starting/resulting element index
 * 
 * @return Status code of type trackExitCodes
 * 
 * Uses the Leaded Celestial Walk algorithm to find the element that contains the target point.
 * This is a variant of the Celestial Walk algorithm that goes from element A containing the starting 
 * point x0 (in one side) to element B containing the target point xt. To enforce that the next element is
 * as close the the segment AB as possible, the next point at each iteration is computed
 * along segment AB, but just after the current element.
 */


__device__ __forceinline__ bool celestialWalk(Element *&curElem, Point3 x0, Point3 &xt, const Side *origSide, Side *&lastCrossSide, int jlev) {

   constexpr int MAX_ITERATIONS = 40;

   int iter = 0;
   double fac = 0.00001;
   while (iter < MAX_ITERATIONS) {
      iter++;

      Bool4 isLeft, isObtuse, isBisector;
      curElem->computePredicates(toPoint2(x0), toPoint2(xt), fac, origSide, isLeft, isObtuse, isBisector, jlev);
// showPoints(x0, xt, jlev, origSide, 0xFF0000u, 'x', 1, 0, 0);
      Bool4 isValid = !isLeft & (!isObtuse | isBisector); // bitwise operation are overloaded


      if (isValid.none()) 
      // Inside the element if none of the edges is a valid continuation
      {
         // Just returns true and the element (curElem) is unchanged).
         // This means that the algorithm has found the element that contains the target point.
         return true;
      }
      else { // Outside the element

         // Outside this element. Pick the first valid exit side and cross.
         lastCrossSide = curElem->sides[isValid.firstOne()];
         Element *prevElem = curElem;
         curElem = lastCrossSide->otherElem(curElem);

         // If it crosses into a dry element or a boundary element, return false. This means that the track will now move 
         // along the boundary.
         if (curElem == nullptr || curElem->isDry) {
         
            curElem = prevElem;
            return false;
         }
         
      }

   }

  assert(false && "Leaded Celestial Walk: exceeded max iterations");
}



//===============================================================================

__device__ void lock(int *mutex) {
    while (atomicCAS(mutex, 0, 1) != 0); // spin until lock acquired
}

__device__ void unlock(int *mutex) {
    atomicExch(mutex, 0);
}

__device__ void addDebugPlot(Vector3 x1, Vector3 x2, unsigned int color = 0xFF0000, char marker = 'o', int width = 1, int dashLen = 0, int label = -1) {
// return;   
   // Lock the mutex
   lock(&debugPlotMutex);
 
   bool ok = true;
   if (debugPlotIdx == nullptr) ok = false;
   if (*debugPlotIdx >= 1000000) ok = false;
   if (ok) {
      debugPlot[(*debugPlotIdx) + 0*1000000] = x1.x;
      debugPlot[(*debugPlotIdx) + 1*1000000] = x1.y;
      debugPlot[(*debugPlotIdx) + 2*1000000] = x2.x;
      debugPlot[(*debugPlotIdx) + 3*1000000] = x2.y;
      debugPlot[(*debugPlotIdx) + 4*1000000] = (double)color;
      debugPlot[(*debugPlotIdx) + 5*1000000] = (double)(((int)(unsigned char)marker & 0xFF) | ((dashLen & 0xFF) << 8) | ((width & 0xFF) << 16));
      debugPlot[(*debugPlotIdx) + 6*1000000] = (double)label;
      (*debugPlotIdx)++;
   }

   __threadfence();               // ensure all above writes are visible
   // Unlock the mutex
   unlock(&debugPlotMutex);

}

__device__ void addDebugPlot(Vector2 x1, Vector2 x2, unsigned int color = 0xA0A0A0, char marker = 'o', int width = 1, int dashLen = 0, int label = -1) {
   addDebugPlot(toPoint3(x1, 0.0), toPoint3(x2, 0.0), color, marker, width, dashLen, label);
}



void showPoints(Point3 x1, Point3 x2, int jlev, const Side *origSide, unsigned int color, char marker, int width, int dashLen, int label) {
   // if ((origSide->idx == 1276362 || origSide->idx == 1276363 || origSide->idx == 1276361 || origSide->idx == 0 || origSide->idx == 9) && jlev == 48) {

   if (true) {
      // Convert points back to original frame
      x1.z = 0.0;
      x2.z = 0.0;
      Point3 X1 = origSide->frame.L2G(x1) + toPoint3(origSide->center, origSide->centerZ());// + origSide->center3D(jlev);
      Point3 X2 = origSide->frame.L2G(x2) + toPoint3(origSide->center, origSide->centerZ());// + origSide->center3D(jlev);
 //printf("wwwwww origSide->idx = %d, marker = %c, width = %d, dashLen = %d, x1 = (%.17f, %.17f, %.17f), x2 = (%.17f, %.17f, %.17f)\n", origSide->idx, marker, width, dashLen, x1.x, x1.y, x1.z, x2.x, x2.y, x2.z);

      addDebugPlot(X1, X2, color, marker, width, dashLen);
   }
}

//===============================================================================
//     
//===============================================================================


__device__ __forceinline__ Vector2 clampVelocity(Vector2 v) {
   return {clamp(v.x, -rmaxvelx, rmaxvelx), clamp(v.y, -rmaxvely, rmaxvely)};
}


__device__ __forceinline__ void update_bdef(double time, double x0, double y0, double &dep, Vector3 &vel) {
   // For imm=2, user needs to update bottom depth and velocity

   dep = min(1.0, 7.0 - (x0 + time));
   vel = { -1.0, 0.0, 0.0 };

}


__global__ void setupMatrices(
   double** data_ptrs,
   int**    indices_ptrs,
   int**    indptr_ptrs,
   int*     nrows,
   int*     ncols,
   double** x_ptrs,
   double** b_ptrs) {
      
   CsrBlockMatrix blocks[nPartitions];
   for (int q = 0; q < nPartitions; ++q)
      blocks[q] = CsrBlockMatrix(data_ptrs[q], indices_ptrs[q], indptr_ptrs[q], nrows[q], ncols[q]);
   mesh->A = new RowBlockMatrix(blocks);
   mesh->x = new FullVector(x_ptrs);
   mesh->b = new FullVector(b_ptrs);
}


__inline__ __device__ Point2 Element::outputCoords(void) const {
   if constexpr(ics == 1)  {
      return center;

   } else {
      const Point2 p0 = nodes[0]->outputCoords();
      const Point2 p1 = nodes[1]->outputCoords();
      const Point2 p2 = nodes[2]->outputCoords();

      if (isTri()) {
         return (1.0/3.0)*(p0 + p1 + p2);

      } else {
         const Point2 p3 = nodes[3]->outputCoords();
         return (1.0/4.0)*(p0 + p1 + p2 + p3);
      }

   }
}


__device__ __forceinline__ double Element::centerZ(void) const { 
// Z coordinate of the element center in the global cartesian coordinates.
// This works because we also compute the average of the x,y coordinates of the nodes to compute the 3D center.

   if constexpr(ics==1) { 
      return 0.0; 

   } else { 

      if (isTri())  return (1.0/3.0)*(nodes[0]->gz + nodes[1]->gz + nodes[2]->gz);
      else          return (1.0/4.0)*(nodes[0]->gz + nodes[1]->gz + nodes[2]->gz + nodes[3]->gz);
   } 
}


   template<typename T>
   __device__ __forceinline__ void Element::getNodalValues(T Node::* field, T res[4]) const {
      // Get the values of the field at the nodes of the element
      res[0] = nodes[0]->*field;
      res[1] = nodes[1]->*field;
      res[2] = nodes[2]->*field;
      if (isQuad()) res[3] = nodes[3]->*field;
      else          res[3] = T();
   }

   template<typename T>
   __device__ __forceinline__ void Element::getNodalValues(int level, T (Node::* field)[nvrt], T res[4]) const {

      assert(level < nvrt && "Level is out of range");
      // Get the values of the field at the nodes of the element
      res[0] = (nodes[0]->*field)[max(level, nodes[0]->kb)];
      res[1] = (nodes[1]->*field)[max(level, nodes[1]->kb)];
      res[2] = (nodes[2]->*field)[max(level, nodes[2]->kb)];
      if (isQuad()) res[3] = (nodes[3]->*field)[max(level, nodes[3]->kb)];
      else          res[3] = T();
   }    

   template<typename T>
   __device__ __forceinline__ void Element::getNodalValues(int level, T (Node::*func)(int) const, T res[4]) const {

      assert(level < nvrt && "Level is out of range");

      // Get the values of the func at the nodes of the element
      res[0] = (nodes[0]->*func)(max(level, nodes[0]->kb));
      res[1] = (nodes[1]->*func)(max(level, nodes[1]->kb));
      res[2] = (nodes[2]->*func)(max(level, nodes[2]->kb));
      if (isQuad()) res[3] = (nodes[3]->*func)(max(level, nodes[3]->kb));
      else          res[3] = T();  // This is 0.0 or a zero vector

   }      

   template<typename T>
   __device__ __forceinline__ void Element::getSideValues(T Side::* field, T res[4]) const {
      // Get the values of the field at the nodes of the element
      res[0] = sides[0]->*field;
      res[1] = sides[1]->*field;
      res[2] = sides[2]->*field;
      if (isQuad()) res[3] = sides[3]->*field;
      else          res[3] = T();
   }

   template<typename T>
   __device__ __forceinline__ void Element::getSideValues(int level, T (Side::* field)[nvrt], T res[4]) const {

      assert(level < nvrt && "Level is out of range");

      // Get the values of the field at the nodes of the element
      res[0] = (sides[0]->*field)[max(level, sides[0]->kb)];
      res[1] = (sides[1]->*field)[max(level, sides[1]->kb)];
      res[2] = (sides[2]->*field)[max(level, sides[2]->kb)];
      if (isQuad()) res[3] = (sides[3]->*field)[max(level, sides[3]->kb)];
      else          res[3] = T();
   }    

   template<typename T>
   __device__ __forceinline__ void Element::getSideValues(int level, T (Side::*func)(int) const, T res[4]) const {

      assert(level < nvrt && "Level is out of range");

      // Get the values of the func at the nodes of the element
      res[0] = (sides[0]->*func)(max(level, sides[0]->kb));
      res[1] = (sides[1]->*func)(max(level, sides[1]->kb));
      res[2] = (sides[2]->*func)(max(level, sides[2]->kb));
      if (isQuad()) res[3] = (sides[3]->*func)(max(level, sides[3]->kb));
      else          res[3] = T();  // This is 0.0 or a zero vector

   }      

   
   __forceinline__ __device__ void Element::projectNodeVectorsToElementFrame(Vector2 *v) const {
   // Project the vectors from the node frame to the element frame

 
      v[0] = nodes[0]->frame.projHorVec(v[0], frame);
      v[1] = nodes[1]->frame.projHorVec(v[1], frame);
      v[2] = nodes[2]->frame.projHorVec(v[2], frame);
      if (isQuad()) v[3] = nodes[3]->frame.projHorVec(v[3], frame);
      else          v[3] = Vector2();
   }

   __forceinline__ __device__ void Element::projectSideVectorsToElementFrame(Vector2 *v) const {
   // Project the vectors from the side frame to the element frame

      v[0] = sides[0]->frame.projHorVec(v[0], frame);
      v[1] = sides[1]->frame.projHorVec(v[1], frame);
      v[2] = sides[2]->frame.projHorVec(v[2], frame);
      if (isQuad()) v[3] = sides[3]->frame.projHorVec(v[3], frame);
      else          v[3] = Vector2();
   }


   template<typename T>
   __device__ __forceinline__ T Element::interpolate(int level, const double shapeCoords[4], T (Node::*func)(int) const) const {
      T values[4];
      getNodalValues(level, func, values);

      // If the field is a horizontal vector and the ics is 2, project the vector from the node frame to the element frame
      if constexpr (IsVector2<T>::value && ics == 2) {
// printf("interpolate 1: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
         projectNodeVectorsToElementFrame(values);
printf("interpolate 2: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
      }

      return shapeCoords[0]*values[0] + shapeCoords[1]*values[1] + shapeCoords[2]*values[2] + shapeCoords[3]*values[3];
   }

   template<typename T>
   __device__ __forceinline__ T Element::interpolate(const double shapeCoords[4], T Node::* field) const {
      T values[4];
      getNodalValues(field, values);

      // If the field is a horizontal vector and the ics is 2, project the vector from the node frame to the element frame
      if constexpr (IsVector2<T>::value && ics == 2) {
// printf("interpolate 1: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
         projectNodeVectorsToElementFrame(values);
printf("interpolate 3: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
      }
      
      return shapeCoords[0]*values[0] + shapeCoords[1]*values[1] + shapeCoords[2]*values[2] + shapeCoords[3]*values[3];
   }


   template<typename T>
   __device__ __forceinline__ T Element::interpolate(int level, const double shapeCoords[4], T (Node::* field)[nvrt]) const {
      T values[4];
      getNodalValues(level, field, values);
      
if (shapeCoords[0] > 1.1 || shapeCoords[1] > 1.1 || shapeCoords[2] > 1.1 || 
    shapeCoords[0] < -0.1 || shapeCoords[1] < -0.1 || shapeCoords[2] < -0.1) printf("XXX shapeCoords = (%.17f, %.17f, %.17f, %.17f)\n", shapeCoords[0], shapeCoords[1], shapeCoords[2], shapeCoords[3]);
      // If the field is a horizontal vector and the ics is 2, project the vector from the node frame to the element frame
      if constexpr (IsVector2<T>::value && ics == 2) {
// printf("interpolate 1: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
         projectNodeVectorsToElementFrame(values);
// printf("interpolate 4: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
      }

      return shapeCoords[0]*values[0] + shapeCoords[1]*values[1] + shapeCoords[2]*values[2] + shapeCoords[3]*values[3];
   }

   template<typename T>
   __device__ __forceinline__ T Element::mean(T Node::* field) const {
      T values[4];
      getNodalValues(field, values);

      // If the field is a horizontal vector and the ics is 2, project the vector from the node frame to the element frame
      if constexpr (IsVector2<T>::value && ics == 2) {
// printf("mean 1: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
         projectNodeVectorsToElementFrame(values);
// printf("mean 2: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
      }

      return (values[0] + values[1] + values[2] + values[3])/i34;
   }

   template<typename T>
   __device__ __forceinline__ T Element::mean(T Side::* field) const {
      T values[4];
      getSideValues(field, values);
      
      // If the field is a horizontal vector and the ics is 2, project the vector from the node frame to the element frame
      if constexpr (IsVector2<T>::value && ics == 2) {
// printf("mean 1: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
         projectSideVectorsToElementFrame(values);
// printf("mean 2: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
      }

      return (values[0] + values[1] + values[2] + values[3])/i34;
   }   

   template<typename T>
   __device__ __forceinline__ T Element::mean(int level,  T (Side::*func)(int) const) const {

      T values[4];


      values[0] = (sides[0]->*func)(max(level, sides[0]->kb));
      values[1] = (sides[1]->*func)(max(level, sides[1]->kb));
      values[2] = (sides[2]->*func)(max(level, sides[2]->kb));
      if (isQuad()) values[3] = (sides[3]->*func)(max(level, sides[3]->kb));
      else          values[3] = T();

      if constexpr (IsVector2<T>::value && ics == 2) projectSideVectorsToElementFrame(values);

      return (values[0] + values[1] + values[2] + values[3])/i34;
   }


   template<typename T>
   __device__ __forceinline__ T Element::mean(int level, T (Node::* field)[nvrt]) const {
      T values[4];
      getNodalValues(level, field, values);

      // If the field is a horizontal vector and the ics is 2, project the vector from the node frame to the element frame
      if constexpr (IsVector2<T>::value && ics == 2) {
// printf("mean 1: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
         projectNodeVectorsToElementFrame(values);
// printf("mean 2: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
      }

      return (values[0] + values[1] + values[2] + values[3])/i34;
   }

   template<typename T>
   __device__ __forceinline__ T Element::mean(int level, T (Side::* field)[nvrt]) const {
      T values[4];
      getSideValues(level, field, values);

      // If the field is a horizontal vector and the ics is 2, project the vector from the node frame to the element frame
      if constexpr (IsVector2<T>::value && ics == 2) {
// XXXXX printf("mean 1: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
         projectSideVectorsToElementFrame(values);
// printf("mean 2: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
      }

      return (values[0] + values[1] + values[2] + values[3])/i34;
   }   



__device__ __forceinline__ Point2 Element::getLocalHorizontalNodePosition(int idxNode) const {
// Local horizontal node position in the element frame. idxNode = 0,1,2,3
// "Local" here means that the coordinates are in the element frame where the element is horizontal.
// Only when ics == 2, the origin is at the element center.

   assert(idxNode < i34 && idxNode >= 0 && "getLocalHorizontalNodePosition: idxNode out of range");

   const Node *node = nodes[idxNode];

   if constexpr(ics==1) {
      return node->p;

   } else {
      const Vector3 c = {center.x, center.y, centerZ()}; // Element center

      const Point3 p = toPoint3(node->p, node->gz); // XXXXX use p3D or similar?

      return toVector2(frame.G2L(p - c));  // Vectors are implicitly cast into points (in this case at the local frame centered at c)
   }

}


__device__ __forceinline__ void Element::getLocalHorizontalNodePositionsTri(Point2 &x0, Point2 &x1, Point2 &x2, 
                                                                            const Vector3& c, const Frame<ics>& dstFrame, int idx) const {
// Performs "getLocalHorizontalNodePosition" for the three nodes of the triangle.

   const int id1 = idx;
   const int id2 = nxq(1 + id1);
   const int id3 = nxq(2 + id1);

   if constexpr(ics==1) {
      x0 = nodes[id1]->p;
      x1 = nodes[id2]->p;
      x2 = nodes[id3]->p;

   } else {

      const Point3 p0 = toPoint3(nodes[id1]->p, nodes[id1]->gz); // XXXXX use p3D or similar?
      const Point3 p1 = toPoint3(nodes[id2]->p, nodes[id2]->gz);
      const Point3 p2 = toPoint3(nodes[id3]->p, nodes[id3]->gz);

      x0 = toVector2(dstFrame.G2L(p0 - c));  // Vectors are implicitly cast into points at the local frame centered at c
      x1 = toVector2(dstFrame.G2L(p1 - c));
      x2 = toVector2(dstFrame.G2L(p2 - c));

   }

}

__device__ __forceinline__ void Element::getLocalHorizontalNodePositionsTri(Point2 &x0, Point2 &x1, Point2 &x2, int idx) const {
   getLocalHorizontalNodePositionsTri(x0, x1, x2, {center.x, center.y, centerZ()}, frame, idx);
}

__device__ __forceinline__ void Element::getLocalHorizontalNodePositionsQuad(
            Point2 &x0, Point2 &x1, Point2 &x2, Point2 &x3,
            const Vector3& c, const Frame<ics>& dstFrame) const {
// Performs "getLocalHorizontalNodePosition" for the four nodes of the quad.

   assert(isQuad() && "getLocalHorizontalNodePositionsQuad only works for quads");



   if constexpr(ics==1) {
      x0 = nodes[0]->p;
      x1 = nodes[1]->p;
      x2 = nodes[2]->p;  
      x3 = nodes[3]->p;

   } else {

      const Point3 p0 = toPoint3(nodes[0]->p, nodes[0]->gz); // XXXXX use p3D or similar?
      const Point3 p1 = toPoint3(nodes[1]->p, nodes[1]->gz);
      const Point3 p2 = toPoint3(nodes[2]->p, nodes[2]->gz);
      const Point3 p3 = toPoint3(nodes[3]->p, nodes[3]->gz);

      x0 = toVector2(dstFrame.G2L(p0 - c));  // Vectors are implicitly cast into points at the local frame centered at c
      x1 = toVector2(dstFrame.G2L(p1 - c));
      x2 = toVector2(dstFrame.G2L(p2 - c));
      x3 = toVector2(dstFrame.G2L(p3 - c));
   }

}

__device__ __forceinline__ void Element::getLocalHorizontalNodePositionsQuad(
   Point2& x0, Point2& x1, Point2& x2, Point2& x3) const
{
   getLocalHorizontalNodePositionsQuad(x0, x1, x2, x3, {center.x, center.y, centerZ()}, frame);
}


__device__ __forceinline__ Vector3 Element::n(int k, double &area3D) const {
// Computes the *local* (element frame) 3D normal vector to the element at level k
// Returns the normal vector and the area of the element at this level
   
   double nx, ny, nz;

   if (isTri()) {
      // For triangles...

      // Indices of the nodes at level k, taking care of not going below the bed
      const int k0 = max(k, nodes[0]->kb);
      const int k1 = max(k, nodes[1]->kb);
      const int k2 = max(k, nodes[2]->kb);

      // Gets the element frame coordinates of the nodes at level k
      Point3 x0, x1, x2;
      if constexpr(ics==1) { // XXX/ I think that we can use always ics==2 because it is more general.
         x0 = toPoint3(nodes[0]->p, nodes[0]->z[k0]);
         x1 = toPoint3(nodes[1]->p, nodes[1]->z[k1]);
         x2 = toPoint3(nodes[2]->p, nodes[2]->z[k2]);

      } else {  // ics == 2
         Point2 h0, h1, h2;

         getLocalHorizontalNodePositionsTri(h0, h1, h2);

         x0 = toPoint3(h0, nodes[0]->z[k0]);
         x1 = toPoint3(h1, nodes[1]->z[k1]);
         x2 = toPoint3(h2, nodes[2]->z[k2]);
      }

      const Vector3 n20 = x2 - x0;
      const Vector3 n10 = x1 - x0; 

      // Computes the cross product of vectors n20 and n10. 
      nx = (n10.y)*(n20.z) - (n20.y)*(n10.z);
      ny = (n20.x)*(n10.z) - (n10.x)*(n20.z);
      nz = 2.0*area; // = (n10.x)*(n20.y) - (n20.x)*(n10.y);

  
   }
   else {
      // For quads...
      Point2 h0, h1, h2, h3;

      getLocalHorizontalNodePositionsQuad(h0, h1, h2, h3);

      // Computes the cross product of the two diagonals (h2 - h0) and (h3 - h1). 
      // WARNING: For the vertical coordinate we are assuming that z = 0 in the element frame coincides with the 0 level of the ellipsoid.
      const Vector3 n20 = toVector3(h2 - h0, nodes[2]->z[k] - nodes[0]->z[k]);
      const Vector3 n31 = toVector3(h3 - h1, nodes[3]->z[k] - nodes[1]->z[k]);

      // Computes the cross product of the two diagonals (h2 - h0) and (h3 - h1). 
      nx = n20.y*n31.z - n20.z*n31.y;
      ny = n20.z*n31.x - n20.x*n31.z;
      nz = n20.x*n31.y - n20.y*n31.x;

   }

   const double len = sqrt(nx*nx + ny*ny + nz*nz);
   area3D = 0.5*len;

   const double inv_len = 1.0/len;

   return {nx*inv_len, ny*inv_len, nz*inv_len};

}


// Helper functions to get the 2D local node position for the Celestial walk
__device__ __forceinline__ Point2 Element::getLocalNodePosition(int idx, const Side* origSide) const {
   const int idx2 = (isTri() && idx == 3) ? 2 : idx;
   const Node *node = getNodeInOrder(idx2);

   double z = 0.0;
   if constexpr(ics == 2) z = node->gz;
   
   Point3 p = toPoint3(node->p, z);

   return toPoint2(origSide->toLocalFrame(p));
}

__device__ __forceinline__ Node *Element::getNodeInOrder(int idx) const {
   return (ssign[idx] > 0) ? sides[idx]->n1 : sides[idx]->n2;
}

__device__ __forceinline__ void Element::computePredicates(const Point2 x0, const Point2& xt, double &fac, const Side* origSide,
                                Bool4& isLeft, Bool4& isObtuse, Bool4& isBisector, int jlev) const
// Computes the three per-side predicates needed for the celestial walk
// (Kuijper, Ermolaev, Devillers, JCGT 7(3), 2018, arXiv:1710.01620).
//
// xt is in local frame coordinates for origSide if ics == 2.
//
// All work is done unconditionally on 4 sides to keep CUDA threads uniform.
// For triangles, side 3 is a sentinel: getNode(3) returns getNode(2), making
// corners[3] == corners[2] and s3 == (0, 0). Predicates that would be
// meaningless on triangles are masked to a safe value (see comments below).
{
   // Project all 4 corner positions once into the local 2D frame.
   // For triangles, C[3] == C[2] (no garbage, but geometrically degenerate).
   // getLocalNodePosition only transforms/projects if ics == 2.
   const Point2 C[4] = {getLocalNodePosition(0, origSide), 
                        getLocalNodePosition(1, origSide), 
                        getLocalNodePosition(2, origSide), 
                        getLocalNodePosition(3, origSide)};

   // Computes the maximum distance from the starting point to the farthest element node projected on the segment x0xt.
   // This has nothing to do with the original celestial walk. It is used in the Leaded Celestial 
   //Walk algorithm to find a point on the segment x0xt that is as close as possible to the starting point,
   // But still farther towards the target point from the current point.
   const Vector2 vt0 = xt - x0;
   double maxProjDist = max(dot(vt0, C[0] - x0), dot(vt0, C[1] - x0));
   maxProjDist = max(maxProjDist, dot(vt0, C[2] - x0)); 
   maxProjDist = isQuad() ? max(maxProjDist, dot(vt0, C[3] - x0)) : maxProjDist;

   fac = max(1.00001*fac, min(1.0, (maxProjDist/dot(vt0, vt0))*1.00001));
   const Vector2 ixt = x0 + (xt - x0)*fac;


// if (origSide->idx == 2860848 && jlev == 46) {
// showPoints(toPoint3(x0,0.0), toPoint3(xt,0.0), jlev, origSide, 0xFF0000u, 'x', 1, 0, 0);
// showPoints(toPoint3(ixt,0.0), toPoint3(ixt,0.0), jlev, origSide, 0x00FF00u, 'o', 1, 0, 0);
// printf("this->idx = %d, origSide->idx = %d, ixt = (%.17f, %.17f), x0 = (%.17f, %.17f), xt = (%.17f, %.17f), fac = %.17f\n", this->idx, origSide->idx, ixt.x, ixt.y, x0.x, x0.y, xt.x, xt.y, fac);
// }
   // Side I goes from corners[I] to corners[nxq(I+1)] in CCW perimeter order.
   // (Note: with corners built from getNode in CCW order, no ssign factor is needed
   //  here the orientation is baked into the corners array.)
   const int i3 = nxq(3);
   const Vector2 s0 = C[ 1] - C[ 0];
   const Vector2 s1 = C[ 2] - C[ 1];
   const Vector2 s2 = C[i3] - C[ 2];
   const Vector2 s3 = C[ 0] - C[i3];
   
   // Side I goes from C[I] to C[nxq(I+1)] in CCW perimeter order.
   // The corners array is built from getNode in CCW order, so the orientation
   // is already correct and no ssign factor is needed in the formulas below.
   //
   // i3 = nxq(3): the index of the corner after side 2.
   //   Triangle: i3 == 0 (the closing corner of the triangle).
   //   Quad:     i3 == 3.
   //
   // The four side direction vectors. For triangles, s3 == (0, 0) since
   // C[i3] == C[0] (the perimeter has already closed by side 2).
   isLeft.set(0,  cross(s0, ixt - C[ 0]) > -small);
   isLeft.set(1,  cross(s1, ixt - C[ 1]) > -small);
   isLeft.set(2,  cross(s2, ixt - C[ 2]) > -small);
   isLeft.set(3, (cross(s3, ixt - C[i3]) > -small) || isTri());  

// printf("idx = %d, cross(s0, xt - C[ 0] ) = %.17f, cross(s1, xt - C[ 1] ) = %.17f, cross(s2, xt - C[ 2] ) = %.17f, cross(s3, xt - C[i3] ) = %.17f\n", idx, cross(s0, xt - C[ 0] ), cross(s1, xt - C[ 1] ), cross(s2, xt - C[ 2] ), cross(s3, xt - C[i3] ));

   // For triangles, the corner after side 2 is the closing corner where
   // side 2 meets side 0 (not side 3, which doesn't exist). The "next side
   // tangent" used by isObtuse and isBisector at that corner is therefore
   // s0 for triangles, s3 for quads.
   const Vector2 s3_tri = isTri() ? s0 : s3;

   // isObtuse.b[I]: true iff the internal face angle at the corner after
   // side I is obtuse. Equivalent to dot of consecutive side directions > 0
   // (acute "turn angle" along the perimeter = obtuse internal angle).
   //
   // isObtuse.b[3] for triangles: dot(s3, s0) = dot(0, s0) = 0, so the
   // result is false naturally. No sentinel needed.
   isObtuse.set(0, dot(s0, s1) > 0.0);
   isObtuse.set(1, dot(s1, s2) > 0.0);
   isObtuse.set(2, dot(s2, s3_tri) > 0.0);  // triangle: closes back to s0
   isObtuse.set(3, dot(s3, s0) > 0.0);  // Always false for triangles (s3=(0,0))

   // isBisector.b[I]: at the corner after side I, this tests whether xt is
   // on the "exit through side I" half of the corner cone, separated by the
   // approximate angular bisector of the celestial walk. The bisector is
   // perpendicular to the chord across the corner; the test is the sign of
   // (xt - corner_vertex) projected onto that chord.
   //
   // isBisector.b[3] is always masked out by isLeft.b[3] = true for triangles
   // (the validity rule sets valid_3 = false), so its computed value for
   // triangles doesn't matter. No sentinel needed.
   isBisector.set(0, dot(ixt - C[ 1], s0 + s1) < 0.0);
   isBisector.set(1, dot(ixt - C[ 2], s1 + s2) < 0.0);
   isBisector.set(2, dot(ixt - C[i3], s2 + s3_tri) < 0.0);
   isBisector.set(3, dot(ixt - C[ 0], s3 + s0) < 0.0);   

}


__inline__ __device__ double Element::computeSignedArea(void) const {
   // Compute the area of the element

   Point2 x0, x1, x2, x3;
   getLocalHorizontalNodePositionsTri(x0, x1, x2);


   double res = signedArea(x0, x1, x2);

   // In the quad case we need to add the area of the other triangle.
   if (isQuad()) {
      x3 = getLocalHorizontalNodePosition(3);
      double ar1 = signedArea(x0, x2, x3);
      double ar2 = signedArea(x0, x1, x3);
      double ar3 = signedArea(x1, x2, x3);
      assert (ar1 >= 0.0 && ar2 >= 0.0 && ar3 >= 0.0 && "Negative area found in quad element");

      res += ar1;

   }

   return res;
}


__inline__ __device__ bool Element::isDry2(void) const { return (nodes[0]->isDry2()) ||
                                                                (nodes[1]->isDry2()) ||
                                                                (nodes[2]->isDry2()) ||
                                                                (isQuad() && nodes[3]->isDry2());
}

__inline__ __device__ double Element::detJ_xi_quad(Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const {
   return (1.0/8.0)*((x1.x - x0.x)*(x2.y - x3.y) - (x2.x - x3.x)*(x1.y - x0.y));
}

__inline__ __device__ double Element::detJ_eta_quad(Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const {
   return (1.0/8.0)*((x2.x - x1.x)*(x3.y - x0.y) - (x3.x - x0.x)*(x2.y - x1.y));
}

__inline__ __device__ bool Element::isEdgeBoundary(Side *side) const { 
   // In boundary sides, the element 1 is nullptr
   return side->elems[1] == nullptr; 
}

__forceinline__ __device__ bool Element::inverseBilinearQuad(Point2 p, double &xi, double &eta, double shp[4]) const {
   // Given a point p in the quad element returns the bilinear coordinates (xi,eta) and the shape coordinates shp[4]
   // For ics == 2, this function uses local frame coordinates.

   assert(isQuad());

   // In bilinear coordinates, a point of coordinates (xi,eta) is given by: 
   // p = c + 0.25*(xi*vxi + eta*veta + xi*eta*vb)
   // where vxi, veta, vb are vectors and c is the center of the element.

   Point2 p0, p1, p2, p3;

   getLocalHorizontalNodePositionsQuad(p0, p1, p2, p3);

   // A point of bilinear coordinates (xi,eta) has cartesian coordinates given by:
   // p = c + 0.25 * (xi*vxi + eta*veta + xi*eta*vb)
   // vectors vxi and veta connect the centers of opposite sides
   // Vector vb is the vector from actual p3 to the p3 that would complete a parallelogram
   // vb is zero for parallelogram elements.
   const Point2 c = 0.25*(p0 + p1 + p2 + p3); // Center of the element (in local frame coordinates if ics == 2)

   //  The bilinear mapping is:
   //      P(s,t) = C + s*U + t*V + s*t*W
   //
   //  where:
   //    C = centroid (average of four vertices)
   //    U = half-vector connecting midpoints of opposite sides (s-direction)
   //    V = half-vector connecting midpoints of opposite sides (t-direction)
   //    W = bilinear correction; measures deviation from parallelogram
   //        +W at P0,P2 (where s*t=+1), -W at P1,P3 (where s*t=-1)
   //        W = 0 for parallelogram elements.
   //
   //  Key idea: project H = P - C onto a direction nb parallel to W.
   //
   //  cross(nb, W) = 0 because nb || W, so the cross product kills the
   //  nonlinear s*t*W term completely, giving a PURE LINEAR equation:
   //
   //      cross(nb, H) = s * cross(nb, U) + t * cross(nb, V)      (L)
   //
   //  Geometrically, (L) projects everything onto the direction
   //  perpendicular to W, where the mapping is purely affine.
   //
   //  The dot product with nb retains the nonlinear term:
   //
   //      dot(nb, H) = s*dot(nb,U) + t*dot(nb,V) + s*t*dot(nb,W)  (Q)
   //
   //  Geometrically, (Q) projects onto the W direction, where
   //  the nonlinearity lives.
   //
   //  We solve (L) for whichever variable has the larger coefficient
   //  (better conditioned back-substitution), substitute into (Q),
   //  and get a scalar quadratic in the other variable.
   // ================================================================

   // ---- Basis vectors ----
   // U: s-direction. 2*U connects midpoint(P0,P3) to midpoint(P1,P2).
   const Vector2 vxi  = (p1 + p2) - (p0 + p3);   // = 4*U
   // V: t-direction. 2*V connects midpoint(P0,P1) to midpoint(P2,P3).
   const Vector2 veta = (p2 + p3) - (p0 + p1);   // = 4*V
   // W: nonlinear correction. Difference between quad and parallelogram.
   // |vb|/max(|vxi|,|veta|) measures the "non-parallelogramness".
   const Vector2 vb   = (p0 + p2) - (p1 + p3);   // = 4*W

   // ---- H: displacement from centroid (scaled by 4 for convenience) ----
   // Working with 4*H avoids the 0.25 multiplications on U, V, W.
   // The equation becomes:  4*H = xi*(4*U) + eta*(4*V) + xi*eta*(4*W)
   // which is just:  4*dp = xi*vxi + eta*veta + xi*eta*vb
   const Vector2 dp4 = 4.0*(p - c);   // = 4*H

   // ---- nb: direction parallel to W, L1-normalized ----
   // We don't need a unit vector -- just a well-scaled direction parallel
   // to W.  L1 norm avoids sqrt entirely (exact to 1 ulp, no rounding
   // from sqrt).  When W ~ 0 (parallelogram), the quadratic coefficient
   // Qw ~ 0 and the solver's linear branch handles it correctly.
   //
   // Branchless: nb = vb / (|vb.x| + |vb.y| + tiny).
   // When vb ~ 0, nb ~ 0, all L/Q coefficients vanish proportionally,
   // and the quadratic solver takes its linear path (a ~ 0).
   const double wlen1 = fabs(vb.x) + fabs(vb.y) + 1e-30;
   const double inv_wlen = 1.0 / wlen1;
   const Vector2 nb = vb * inv_wlen;

   // ---- Linear equation (L): Lp = Ls*xi + Lt*eta ----
   // cross(nb, .) projects onto the direction PERPENDICULAR to W.
   // Since nb || W, cross(nb, vb) = 0 exactly, killing the xi*eta term.
   // What remains is a pure affine relationship between xi and eta.
   const double Lp = cross(nb, dp4);     // cross(nb, 4*H)
   const double Ls = cross(nb, vxi);     // cross(nb, 4*U)
   const double Lt = cross(nb, veta);    // cross(nb, 4*V)

   // ---- Quadratic equation (Q): Qp = Qs*xi + Qt*eta + Qw*xi*eta ----
   // dot(nb, .) projects onto the W direction, where the nonlinearity lives.
   // dot(nb, vb) = |vb|^2/(|vb.x|+|vb.y|) != 0 for non-parallelograms.
   const double Qp = dot(nb, dp4);       // dot(nb, 4*H)
   const double Qs = dot(nb, vxi);       // dot(nb, 4*U)
   const double Qt = dot(nb, veta);      // dot(nb, 4*V)
   const double Qw = dot(nb, vb);        // dot(nb, 4*W)

   // ---- Solve: pick the better-conditioned variable from (L) ----
   // Back-substitute from (L) dividing by the LARGER coefficient,
   // minimizing round-off amplification.
   //
   // Commented below is the original code. It is more clear, but is branched. Below is a faster versions
   // that uses predicated moves.
   //
   // if (fabs(Ls) > fabs(Lt)) {
   //    // Solve quadratic in eta, back-substitute for xi from (L).
   //    // From (L):  xi = (Lp - Lt*eta) / Ls
   //    // Substitute into (Q), multiply by Ls, group by eta:
   //    const double A = -Qw * Lt;
   //    const double B = -Qs * Lt + Qt * Ls + Qw * Lp;
   //    const double C =  Qs * Lp - Qp * Ls;

   //    double eta1, eta2;
   //    solveQuadratic(A, B, C, eta1, eta2);

   //    // Pick root with smaller |eta| (domain centered at 0).
   //    // Ternary-> predicated move, no branch divergence.
   //    eta = (fabs(eta1) < fabs(eta2)) ? eta1 : eta2;

   //    // Back-substitute: divides by |Ls| >= |Lt|, best conditioned.
   //    xi = (Lp - Lt * eta) / Ls;

   // } else {
   //    // Solve quadratic in xi, back-substitute for eta from (L).
   //    // From (L):  eta = (Lp - Ls*xi) / Lt
   //    // Substitute into (Q), multiply by Lt, group by xi:
   //    const double A = -Qw * Ls;
   //    const double B =  Qs * Lt - Qt * Ls + Qw * Lp;
   //    const double C =  Qt * Lp - Qp * Lt;

   //    double xi1, xi2;
   //    solveQuadratic(A, B, C, xi1, xi2);

   //    xi = (fabs(xi1) < fabs(xi2)) ? xi1 : xi2;

   //    eta = (Lp - Ls * xi) / Lt;
   // }
   const bool solve_s = fabs(Lt) >= fabs(Ls);   // true-> quadratic gives xi
   const double La = solve_s ? Lt : Ls;          // larger  |L|, for back-sub
   const double Lb = solve_s ? Ls : Lt;          // smaller |L|, in quadratic
   const double Qa = solve_s ? Qt : Qs;          // Q-coeff matching La
   const double Qb = solve_s ? Qs : Qt;          // Q-coeff matching Lb

   // Quadratic:  A*r^2 + B*r + C = 0
   // (same formula for both cases, just with swapped coefficients)
   const double A = -Qw * Lb;
   const double B = -Qa * Lb + Qb * La + Qw * Lp;
   const double C =  Qa * Lp - Qp * La;

   double r1, r2;
   solveQuadratic(A, B, C, r1, r2);

   // Pick root closer to 0 (domain is [-1,1] centered at origin)
   const double r = (fabs(r1) < fabs(r2)) ? r1 : r2;

   // Back-substitute into (L): other = (Lp - Lb*r) / La
   const double other = (Lp - Lb * r) / La;

   // Assign to (xi, eta) depending on which variable the quadratic solved
   xi  = solve_s ? r     : other;
   eta = solve_s ? other : r;


   // ---- Clamp to [-1, 1] (safety net for points slightly outside) ----
   xi  = clamp(xi,  -1.0, 1.0);
   eta = clamp(eta, -1.0, 1.0);

   // ---- Shape functions ----
   // N0 = (1-xi)(1-eta)/4  at P0 (xi=-1, eta=-1)
   // N1 = (1+xi)(1-eta)/4  at P1 (xi=+1, eta=-1)
   // N2 = (1+xi)(1+eta)/4  at P2 (xi=+1, eta=+1)
   // N3 = (1-xi)(1+eta)/4  at P3 (xi=-1, eta=+1)
   shp[0] = (1.0 - xi) * (1.0 - eta) * 0.25;
   shp[1] = (1.0 + xi) * (1.0 - eta) * 0.25;
   shp[2] = (1.0 + xi) * (1.0 + eta) * 0.25;
   shp[3] = (1.0 - xi) * (1.0 + eta) * 0.25;

   return true;
}

/**
   Compute the shape coordinates of the point (xt,yt) in the element.
   @param xt (input): x-coordinate of the point
   @param shpCoords (output): shape coordinates of the point
   @param fixNegshpCoord (input): if true, fix negative shape coordinates
   @return true if the point is inside the element, false otherwise
**/
__device__ bool Element::computeShapeCoords(Point2 &xt, double shpCoords[4], EFixIf fixNegshpCoord, const Side *origSide) const {

   // Compute the shape coordinates of the point xt in the element.
   //   If fixNegshpCoord = fiOnlyInside, only fix the shape coordinates if the point is *almost*inside the element.
   //   If fixNegshpCoord = fiAnywhere, we fix the shape coordinates if they are negative.
   // It's not reliable to use the area coord for ics=2 when the pt is far away
   // from the local frame, so use additional check.


   Point3 origSideCenter = toPoint3(origSide->center, origSide->centerZ());

   if (isTri()) { // Triangular element

      // Horizontal local node coordinates. 
      Point2 x0, x1, x2;
      getLocalHorizontalNodePositionsTri(x0, x1, x2, origSideCenter, origSide->frame);

      shpCoords[0] = signedArea(xt,x1,x2)/area;
      shpCoords[1] = signedArea(x0,xt,x2)/area; 
      shpCoords[2] = 1.0 - shpCoords[0] - shpCoords[1];
      shpCoords[3] = 0.0;

      double tmpmin = min3(shpCoords[0], shpCoords[1], shpCoords[2]);
      const bool res = (tmpmin >= -small);  
      
      // When to fix the negative shape coordinates?
      const bool fix = (fixNegshpCoord == EFixIf::fiOnlyInside && res) || (fixNegshpCoord == EFixIf::fiAnywhere);

      // If at least one of the shape coordinates is negative, and we have the flag fixNegshpCoord set,
      // then we fix the negative shape coordinates.
      if ((tmpmin < 0) && fix) {

         // NOTE: I have changed this algorithm. Now it scales and offsets all the shape coordinates by the same factor,
         // keeping the sum of the shape coordinates equal to 1.
         // I have tested and compared the new algorithm with the old one, and they give
         // very similar results in terms of how close the "fixed" point is to the original point in cartesian coordinates.
         // The new algorithm should be faster since it is branchless.
         const double a = - tmpmin + 1.0e-8;  // XXX Is 1.0e-8 the best value?
         const double b = 1.0/(shpCoords[0] + shpCoords[1] + shpCoords[2] + 3.0*a);

         shpCoords[0] = (shpCoords[0] + a)*b;
         shpCoords[1] = (shpCoords[1] + a)*b;
         shpCoords[2] = (shpCoords[2] + a)*b;          

         assert(!(shpCoords[0] < 0.0 || shpCoords[1] < 0.0 || shpCoords[2] < 0.0 ||
                  shpCoords[0] > 1.0 || shpCoords[1] > 1.0 || shpCoords[2] > 1.0) &&
                "Negative shape coordinates in computeShapeCoords, when fixNegshpCoord is true");

         // Update the point coordinates for the new shape coordinates
         double xcoords[3] = {x0.x, x1.x, x2.x};
         double ycoords[3] = {x0.y, x1.y, x2.y};

         xt = {dotprod3(xcoords, shpCoords), dotprod3(ycoords, shpCoords)};

      }

      // Return true if the point is inside the element
      return res;      
   }
   else {
      assert(isQuad() && "computeShapeCoords");

      // Horizontal local node coordinates with respect to the original side.
      Point2 x0, x1, x2, x3;
      getLocalHorizontalNodePositionsQuad(x0, x1, x2, x3, origSideCenter, origSide->frame);

      const double a0 = signedArea(x0, x1, xt);
      const double a1 = signedArea(x1, x2, xt);
      const double a2 = signedArea(x2, x3, xt);
      const double a3 = signedArea(x3, x0, xt);

      const double minArea = min4(a0, a1, a2, a3);
      
      const bool res = (minArea > -small*area);

      // When to fix the negative shape coordinates?
      const bool fix = (fixNegshpCoord == EFixIf::fiOnlyInside && res) || (fixNegshpCoord == EFixIf::fiAnywhere);

      double xi, eta;    

      // Converts the point xt to the element frame coordinates and then to bilinear+shape coordinates.
      // Does nothing if ics == 1.
      assert(origSide != nullptr);
      Point2 xt_elem = toPoint2(toLocalFrame(origSide->toGlobal(toPoint3(xt, 0.0))));

      inverseBilinearQuad(xt_elem, xi, eta, shpCoords);

      if (fix) {
         double xn[4] = {x0.x, x1.x, x2.x, x3.x};
         double yn[4] = {x0.y, x1.y, x2.y, x3.y};
         xt.x = interpolate(shpCoords, xn);
         xt.y = interpolate(shpCoords, yn);
      }

      return res;
   }


   assert(false && "computeShapeCoords: should never reach this point");
   return false;  // Never reached

}


__device__ __forceinline__ void Element::shapeCoordsToXY(double shpCoords[4], Vector2 &xt) const {
// Convert shape coordinates to cartesian coordinates

   if (isTri()) {
      // Triangular element
assert(false && "shapeCoordsToXY: should it use local frame coordinates?");
      double xcoords[3] = {nodes[0]->p.x, nodes[1]->p.x, nodes[2]->p.x};
      double ycoords[3] = {nodes[0]->p.y, nodes[1]->p.y, nodes[2]->p.y};
   }
   else {
      assert(false && "shapeCoordsToXY: quadrilateral element not implemented");
   }

}



__device__ __forceinline__ double Element::intPhiPhiTri(int i, int j) const {
   // return \int \phi_i*\phi_j dA for tris, where phi_i and phi_j are the shape functions of the nodes i and j.
   
   return ((i == j) ? (2.0/12.0) : (1.0/12.0))*area;
}

__device__ __forceinline__ double Element::intGradPhiGradPhiTri(int i, int j) const {
   // return \int (grad phi_i).(grad phi_j) dA for tris, where phi_i and phi_j are the shape functions of the nodes i and j.
   // gradPhi[i] is constant over a linear triangle, so the integral is area * dot(gradPhi[i], gradPhi[j]).
   assert(isTri() && i >= 0 && j >= 0 && i < 3 && j < 3);
   return area * dot(gradPhi[i], gradPhi[j]);
}

__device__ __forceinline__ double Element::intPhiPhiQuad(int i, int j, Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const {
   // return \int \phi_i*\phi_j dA for quads, where phi_i and phi_j are the shape functions of the nodes i and j.
   // The integral can be computed analytically, using the shape functions in xi and eta coordinates and the Jacobian.
   // The formula for the integral is:
   // (area/16)(1+\xi_i\xi_j/3)(1+\eta_i\eta_j/3) + (detJ_xi/12)(1+\eta_i\eta_j/3)(\xi_i+\xi_j) + (detJ_eta/12)(1+\xi_i\xi_j/3)(\eta_i+\eta_j)\]
   
   const double xi_term  = 1.0 + (xi_n(i) * xi_n(j))/3.0;
   const double eta_term = 1.0 + (eta_n(i)*eta_n(j))/3.0;
   const double eta_sum  = eta_n(i) + eta_n(j);
   const double xi_sum   = xi_n(i)  + xi_n(j);

   const double detJ_xi  = detJ_xi_quad (x0, x1, x2, x3);  // XXXXX put the 8 inside the function??? xi derivative of the Jacobian determinant
   const double detJ_eta = detJ_eta_quad(x0, x1, x2, x3);  // eta derivative of the Jacobian determinant


   return (1.0/16.0)*area*xi_term*eta_term + (1.0/12.0)*(detJ_eta*xi_term*eta_sum + detJ_xi*eta_term*xi_sum);
}
   

__device__ __forceinline__ double Element::intGradPhiGradPhiQuad(int i1, int i2, Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const {
   // return \int (grad phi_i1).(grad phi_i2) dA, where phi_i1 and phi_i2 are the shape functions of the nodes i1 and i2.
   // The integral is computed by using the 2x2 Gauss-Legendre quadrature points in xi and eta coordinates.
   // Again, it uses the determinant of the Jacobian to compute the integral in x,y coordinates.
   // The formula is also based on the identity:
   //   dPhi/dx = (1/det J) * [ (dPhi/dxi)*(dy/deta) - (dPhi/deta)*(dy/dxi)]
   //   dPhi/dy = (1/det J) * [-(dPhi/dxi)*(dx/deta) + (dPhi/deta)*(dx/dxi)]

   assert(isQuad() && i1 >= 0 && i2 >= 0 && i1 < 4 && i2 < 4 && "intGradPhiGradPhi: bad ranges");



   double res = 0.0;
   
   // Quadrature points in xi and eta coordinates
   #pragma unroll
   for (int i = 0; i < 4; i++) {
      // For each quadrature point...

      const Point2 p = GaussPoint(i);  // Quadrature point in (eta,xi) coordinates
      const double eta = p.x;
      const double xi  = p.y;

      // Jacobian entries at the Gauss point p
      const double dx_dxi  = 0.25 * ( (x1.x + x2.x - x0.x - x3.x) + eta*(x0.x + x2.x - x1.x - x3.x) );
      const double dx_deta = 0.25 * ( (x2.x + x3.x - x0.x - x1.x) + xi *(x0.x + x2.x - x1.x - x3.x) );
      const double dy_dxi  = 0.25 * ( (x1.y + x2.y - x0.y - x3.y) + eta*(x0.y + x2.y - x1.y - x3.y) );
      const double dy_deta = 0.25 * ( (x2.y + x3.y - x0.y - x1.y) + xi *(x0.y + x2.y - x1.y - x3.y) );      

      const double detJ = detJ_quad(p, x0, x1, x2, x3);
      assert(detJ > 0.0 && "Negative or zero determinant of the Jacobian in intGradPhiGradPhi");

      // Derivatives of the shape functions i1, i2 at the Gauss point p
      const double dPhi_dxi_1  = dPhi_dxi_quad (i1, p);
      const double dPhi_deta_1 = dPhi_deta_quad(i1, p);
      const double dPhi_dxi_2  = dPhi_dxi_quad (i2, p);
      const double dPhi_deta_2 = dPhi_deta_quad(i2, p);

      // Transform the derivatives of the shape functions to x,y coordinates (this formula leaves a multiplicative factor detJ)
      const double dPhi_dx_1_detJ =  dy_deta * dPhi_dxi_1 - dy_dxi * dPhi_deta_1;
      const double dPhi_dy_1_detJ = -dx_deta * dPhi_dxi_1 + dx_dxi * dPhi_deta_1;
      const double dPhi_dx_2_detJ =  dy_deta * dPhi_dxi_2 - dy_dxi * dPhi_deta_2;
      const double dPhi_dy_2_detJ = -dx_deta * dPhi_dxi_2 + dx_dxi * dPhi_deta_2;

      // Contribution to the integral from the Gauss point for grad(phi_i1).grad(phi_i2)
      res += (dPhi_dx_1_detJ*dPhi_dx_2_detJ + dPhi_dy_1_detJ*dPhi_dy_2_detJ)/detJ;
   }

   return res;
}


/**  Interpolate velocities at a given point (xt,yt,zt)

   @param xt, yt, zt: coordinates of the point to interpolate
   @param uuint: velocity interpolated at the final position (xt,yt,zt)

   This is done in four steps:

   1. Find the shape coordinates for the final position in the element.
   2. Interpolate the z level of the top and bottom levels at xt, yt. 
   3. Perform a binary search to find the vertical level containing zt. It will have to interpolate the 
      z levels at xt, yt at each step.   
   4. Interpolate the velocities at the final position (xt,yt,zt). This is done by interpolating the velocities
      horizontally and vertically.
**/
__device__ Vector3 Element::interpolateVelocities(Point3 xt, const Side *origSide) const {
   // Interpolate the velocities at the point xt. 
   // This function works the same for triangular and quadrilateral elements.
   double shpCoord[4];
   Vector3 u_int;


   // First step is to compute the shape coordinates of (xt,yt)
// printf("interpolateVelocities 1: xt = (%f, %f, %f)\n", xt.x, xt.y, xt.z);
   computeShapeCoords(toVector2(xt), shpCoord, fiAnywhere, origSide);

   assert((min4(shpCoord[0], shpCoord[1], shpCoord[2], shpCoord[3]) >= -small) && "Invalid shape coordinates");

   // Second step is to find the level that contains zt
   // Uses a binary search to find the level
   // Note that z coordinates are interpolated from their values at the nodes of the element.
   // I use L(eft) and R(ight) as usual, but here L is the bottom level and R is the top level
   int L = kb - 1; // In some cases, the bottom level of the element is 1 more than the bottom level of the node XXXX (Check)
   int R = nvrt - 1;

   // Get the z coordinates at the nodes of the element at the levels L and R
   double zLAtNodes[4], zRAtNodes[4], zMAtNodes[4];
   getNodalValues(L, &Node::compZCoord, zLAtNodes);
   getNodalValues(R, &Node::compZCoord, zRAtNodes);


   double zL = interpolate(shpCoord, zLAtNodes);
   double zR = interpolate(shpCoord, zRAtNodes);

   // Interpolates vertically the velocities
   double wght1, wght2;
   if (xt.z <= zL) {
      wght1 = 1.0;
      wght2 = 0.0;
      xt.z = zL;
   }
   else if (xt.z >= zR) {
      wght1 = 0.0;
      wght2 = 1.0;
      xt.z = zR;
   }
   else {
      // Perform a binary search to find the levels between L and R that contain xt.z
      while (R - L > 1) {
         int M = (L + R) / 2.0;

         getNodalValues(M, &Node::compZCoord, zMAtNodes);
         double zM = interpolate(shpCoord, zMAtNodes);
         if (zM < xt.z) 
         {
            L = M;
            zL = zM;
         }
         else 
         {
            R = M;
            zR = zM;
         }
      }
      assert(L == R - 1 && "Binary search failed to find the levels containing zt");
      // At this point, L and R are the two levels that contain zt, and zL and zR are the z coordinates at the nodes of the element at these levels

      wght2 = (xt.z - zL)/(zR - zL);
      wght1 = 1.0 - wght2;

   }
// printf("XXYYX element %d, L = %d, R = %d, zL = %.17f, zR = %.17f, xt.z = %.17f, wght1 = %.17f, wght2 = %.17f\n", idx, L, R, zL, zR, xt.z, wght1, wght2);
// printf("HEREEEEEE 3, wght1 = %.17f, wght2 = %.17f, shpCoord = (%.17f, %.17f, %.17f, %.17f)\n", wght1, wght2, shpCoord[0], shpCoord[1], shpCoord[2], shpCoord[3]);
   // Interpolate horizontally and vertically the velocities at the two levels L, R
   // Remember that interpolate takes frames into account for Vector2 fields and ics==2
   Vector2 uv = wght1*interpolate(L, shpCoord, &Node::u2) + wght2*interpolate(R, shpCoord, &Node::u2);
   double w   = wght1*interpolate(L, shpCoord, &Node::w2) + wght2*interpolate(R, shpCoord, &Node::w2);

   // Project the velocity from the local frame to the origSide frame
   if constexpr(ics == 2) uv = frame.projHorVec(uv, origSide->frame);
   u_int = {uv.x, uv.y, w};

   // Neither of the interpolated velocities should be NaN
   assert(!isnan(u_int.x) && !isnan(u_int.y) && !isnan(u_int.z) && "NaN velocities in interpolateVelocities");

   return u_int;

}




__inline__ __device__ double SideBoundary::getBoundaryArrayValueWeight(void) const {

   // For a side boundary, computes the weight w for which the value at node n1 is weighted by to
   // compute the value at the side center. The one for n2 is computed from 1 - w
   const Node *n1 = side->n1;
   const Node *n2 = side->n2;

   bool n1_match = (n1->nodeBnd->segment == segment);
   bool n2_match = (n2->nodeBnd->segment == segment);

   assert(n1_match || n2_match && "getBoundaryArrayValue: side with no segment boundary");
   if (n1_match && n2_match) {
      return 0.5;
   } else if (n1_match) {
      return 1.0;
   } else if (n2_match) {
      return 0.0;
   }

   return 0.0;
}

__inline__ __device__ double  SideBoundary::etaMeanBnd(void) const { return 0.5*(side->n1->nodeBnd->eta_mean + side->n2->nodeBnd->eta_mean); }         // XXXX Should I use something like getBoundaryArrayValueWeight here?


__inline__ __device__ Vector2 SideBoundary::umo(int jfr) const {

   const double w1 = getBoundaryArrayValueWeight();
   const double w2 = 1.0 - w1;
   const NodeBoundary *n1Bnd = side->n1->nodeBnd;
   const NodeBoundary *n2Bnd = side->n2->nodeBnd;

   return w1*n1Bnd->umo[jfr] + w2*n2Bnd->umo[jfr];
}


__inline__ __device__ Vector2 SideBoundary::ufa(int jfr) const {
   const double w1 = getBoundaryArrayValueWeight();
   const double w2 = 1.0 - w1;
   const NodeBoundary *n1Bnd = side->n1->nodeBnd;
   const NodeBoundary *n2Bnd = side->n2->nodeBnd;

   return w1*n1Bnd->ufa[jfr] + w2*n2Bnd->ufa[jfr];

}


__inline__ __device__ Vector2 SideBoundary::uth(int k) const {
   // Returns the value of the velocity at the side center
   const double w1 = getBoundaryArrayValueWeight();
   const double w2 = 1.0 - w1;
   const NodeBoundary *n1Bnd = side->n1->nodeBnd;
   const NodeBoundary *n2Bnd = side->n2->nodeBnd;

   return w1*n1Bnd->uth[k] + w2*n2Bnd->uth[k];
}


__inline__ __device__ double  SideBoundary::nuth(int k) const {
   // Returns the value of the normal component of the velocity at the side center
   const double w1 = getBoundaryArrayValueWeight();
   const double w2 = 1.0 - w1;
   const NodeBoundary *n1Bnd = side->n1->nodeBnd;
   const NodeBoundary *n2Bnd = side->n2->nodeBnd;

   return w1*n1Bnd->uth[k].x + w2*n2Bnd->uth[k].x;
}


__inline__ __device__ Point2 Side::outputCoords(void) const {
   if constexpr(ics == 1) {
      return center;

   } else {
      const Point2 p0 = n1->outputCoords();
      const Point2 p1 = n2->outputCoords();
      return (1.0/2.0)*(p0 + p1);

   }
}

__inline__ __device__ double  Side::x(double t)   const { assert(t >= 0.0 && t <= 1.0); return (1-t)*n1->p.x + t*n2->p.x; }
__inline__ __device__ double  Side::y(double t)   const { assert(t >= 0.0 && t <= 1.0); return (1-t)*n1->p.y + t*n2->p.y; }
__inline__ __device__ double  Side::z(double t)   const { assert(t >= 0.0 && t <= 1.0); return (1-t)*n1->gz  + t*n2->gz; }
__inline__ __device__ double  Side::lon(void)     const { return 0.5*(n1->lon + n2->lon); }
__inline__ __device__ double  Side::lat(void)     const { return 0.5*(n1->lat + n2->lat); }
__inline__ __device__ Point2  Side::p2_t(double t) const { return {x(t), y(t)}; }
__inline__ __device__ double  Side::centerZ(void) const { return 0.5*(n1->gz + n2->gz); }
__inline__ __device__ Vector2 Side::grad2(double Node::* field, bool (*condition)(const Element*)) const {
// Computes the gradient of the field at the side center.
   Vector2 grad = {0.0, 0.0};
   int count = 0;

   // Goes over the two elements that share this side, computes the gradient of the field at the side center,
   // as the average of the gradients at each element.
   // First element
   if (elems[0] != nullptr && condition(elems[0])) {
      Vector2 elemGrad = elems[0]->grad2(field);

      // For ics=2, projects the element gradient from its frame to the side frame
      if constexpr(ics == 2) elemGrad = elems[0]->frame.projHorVec(elemGrad, frame);

      grad += elemGrad;

      count++;
   }

   // Second element
   if (elems[1] != nullptr && condition(elems[1])) {

      Vector2 elemGrad = elems[1]->grad2(field);

      // For ics=2, projects the element gradient from its frame to the side frame
      if constexpr(ics == 2) elemGrad = elems[1]->frame.projHorVec(elemGrad, frame);

      grad += elemGrad;

      count++;
   }

   // If there are valid elements, returns the average of the gradients at the side's frame.
   if (count > 0) {
      grad = grad*(1.0/count);
   }
   
   return grad;
}

__inline__ __device__ double  Side::grav3(void)   const { return 0.5*(n1->grav2  + n2->grav2); }
__inline__ __device__ int     Side::iadv(void)    const { return  min(n1->iadv,    n2->iadv); }   
__inline__ __device__ double  Side::eta2(void)    const { return 0.5*(n1->eta2   + n2->eta2); }
__inline__ __device__ double  Side::eta1(void)    const { return 0.5*(n1->eta1   + n2->eta1); }
__inline__ __device__ double  Side::w2(int k)     const { return 0.5*(n1->w2[k]  + n2->w2[k]); }
__inline__ __device__ double  Side::depth(void)   const { return 0.5*(n1->depth  + n2->depth); } 
__inline__ __device__ double  Side::dfv(int k)    const { return 0.5*(n1->dfv[k] + n2->dfv[k]); }

__inline__ __device__ double  Side::veg_h(void)   const { return 0.5*(n1->veg_h  + n2->veg_h); }
__inline__ __device__ double  Side::veg_nv(void)  const { return 0.5*(n1->veg_nv + n2->veg_nv); }
__inline__ __device__ double  Side::veg_alpha_vert_mean(void) const { return 0.5*(n1->veg_alpha_vert_mean + n2->veg_alpha_vert_mean); }
__inline__ __device__ double  Side::veg_alpha3D(int k) const { return 0.5*(n1->veg_alpha3D[k] + n2->veg_alpha3D[k]); }


__inline__ __device__ double  Side::Cd(void)     const { if (!isDry) return 0.5*(n1->Cd + n2->Cd); else return 0.0; }

__inline__ __device__ Vector2 Side::tau(void)    const { 
   // Returns the value of the tau at the side center, taking care of the projection from the node frames to the side frame (ics==1 or ics==2)
   if constexpr(ics == 1) { 
      return 0.5*(n1->tau + n2->tau);

   } else {
      // Project the tau from the node frames to the side frame
      const Vector2 stau1 = n1->frame.projHorVec(n1->tau, frame);
      const Vector2 stau2 = n2->frame.projHorVec(n2->tau, frame);

      return 0.5*(stau1 + stau2);
   }
}

/**
 * This method performs backtracking from the side center until a time -dt back in time 
 * from the current time. The backtracking is performed using the velocities (u,v,w) at each
 * intermediate point.
 *
 * We approximate the backtracking curve as a piecewise linear curve. Each segment has
 * a length dtb*vel_step, where vel_step is the module of the velocity at that step and dtb is 
 * the timestep (of the backtracking). dtb is computed as the minimum of dtbk and the 
 * remaining time.
 * 
 * NOTE: Runge-Kutta is not implemented yet. XXX
 * 
 * @param uuint, vvint, wwint: Velocities at the final point
 * @param jlev: Level at which the backtracking is performed
 * @param dt, dtbk: Tracking timesteps
  */
__inline__ __device__ void Side::btrack(
      int jlev,             // Level at which the backtracking is performed XXX rename to k???
      double dt,            // Tracking time
      double dtbk,          // Tracking timestep
      Vector3 &u_int        // Velocities at the initial point
) const {

   // Start from side center
   Point3 x0;

   if constexpr(ics == 1) x0 = {center.x, center.y, zc[jlev]};
   else                   x0 = {0.0, 0.0, zc[jlev]};  // for ics == 2, the origin is the side centern which is (0,0) at the local frame

 
   double time_rm = dt; // Remaining time
   double time_rm_1 = 1e10;
   double time_rm_2 = 2e10;
   double time_rm_3 = 3e10;


   // Initial velocities
   u_int = {u2[jlev].x, u2[jlev].y, w2(jlev)};

   // Find a non dry element for starting the backtracking. At this point there should be at least one.
   Element* curElem = (!elems[0]->isDry) ? elems[0] : elems[1];  // elems[0] is always non-null 

   assert(curElem != nullptr && !curElem->isDry && "Starting element is dry or boundary side");

// showPoints(x0, x0, jlev, this, 0xFF0000u, 'x', 1, 0, 0);


   constexpr int MAX_ITERATIONS = 100; 
   int iter = 0;
   Point3 xt = {-1.0e25, -1.0e25, -1.0e25};
   while(true) {

      // Calculate timestep
      const double dtb = min(dtbk, time_rm);

      // Next point in the backtracking curve
      xt = x0 - dtb*u_int;

      double t, s; // parameters for the line segments p0pt and p1p2
      Side *lastCrossSide = nullptr;
      const bool found = celestialWalk(curElem, x0, xt, this, lastCrossSide, jlev);


      if (time_rm_3 <= time_rm) {
      // If backtracking gets stuck, we stop the backtracking and use the velocity at the last x0.
         xt = x0;
         break;
      }


      if (!found) {
      // This means that the walk reached a boundary or a dry element.

         assert(lastCrossSide!=nullptr && "lastCrossSide is nullptr");

         // Finds the point at the boundary where the line segment x0xt intersects.
         // Here we will assume that it does intersect the boundary segment, even if in some cases,
         // due to floating point errors, it does not exactly intersect.
         // We call p1 and p2 the 2D points of the boundary segment. p0, pt are the 2D versions of x0 and xt.
         // v1 is the vector from p0 to pt, v2 is the vector from p1 to p2.
         // To compute the intersection point we start with the parametric equations.
         // p1 + t*v2 = p0 + s*v1,    t, s in [0,1]
         // We cross-prouct by v1 to cancel the s term:
         // cross(v2, v1)*t = cross(p0 - p1, v1)
         // t = cross(p0 - p1, v1) / cross(v2, v1)

         // Side coordinates in the origSide frame.
         const Point2 p1 = toPoint2(toLocalFrame(toPoint3(lastCrossSide->n1->p, lastCrossSide->n1->gz)));
         const Point2 p2 = toPoint2(toLocalFrame(toPoint3(lastCrossSide->n2->p, lastCrossSide->n2->gz)));

         Point2  p0 = {x0.x, x0.y};
         Point2  pt = {xt.x, xt.y};

         const Vector2 v1 = pt - p0;
         const Vector2 v2 = p2 - p1;



         // We are assuming intersection, so we clamp t to the range [0,1].
         const double denom = cross(v2, v1);
         



         if (abs(denom) < small) {
            // small denominators are fairly common. They occur when the line segment x0xt 
            // is parallel to the boundary segment. This occurs after a "collision" with the
            // boundary, where the backtracking curve startt sliding along the boundary.
            // In this case, we simply assume that the two lines are very close to each other
            // and simply compute the closest point to xt (s = 1.0) on the boundary segment (clamped to p1p2)
            t = clamp(dot(pt - p1, v2) / dot(v2, v2), 0.0, 1.0);
            if (t<0.0000000001) t = 0.0;
            if (t>0.9999999999) t = 1.0;

            // s = clamp(dot(pt - p1, v2) / dot(v2, v2), 0.0, 1.0);
            s = clamp(dot(p1 + t*v2 - p0, v1)/dot(v1, v1), 0.0, 1.0);

         }
         else {
            
            t = clamp(cross(p0 - p1, v1) / denom, 0.0, 1.0);
            if (t<0.0000000001) t = 0.0; // XXXXX use constants
            if (t>0.9999999999) t = 1.0;
            // From the originaL equation we can solve for s:
            // p1 + t*v2 = p0 + s*v1    =>   s*v1 = p1 + t*v2 - p0 
            // Computing the dot product by v1, we get:
            // s = dot((p1 + t*v2 - p0), v1) / dot(v1, v1)
            s = clamp(dot(p1 + t*v2 - p0, v1)/dot(v1, v1), 0.0, 1.0);
            if (s<0.0000000001) s = 0.0; // XXXXX use constants

         }

         assert(curElem!=nullptr && !curElem->isDry && "curElem is nullptr or dry");


         // We update the new target point to be right at the intersection point,
         // and the remaining time must include the time we haven't "used".
         time_rm -= s*dtb;
         xt = toPoint3(p1 + t*v2, x0.z + s*(xt.z - x0.z)); // I use p1,p2 to compute the intersection point to ensure that it is on the boundary, but the z-coordinate is the same as x0.z


         // If it is an open boundary, just exits the loop. No sliding along the boundary is performed.
         if (lastCrossSide->sideBnd!=nullptr && !lastCrossSide->sideBnd->segment->isLand) break;

         // If it is not an open boundary, the velocity for the next step is computed by interpolating the velocities at the intersection point.
         // but extracting only the component along the tangent vector v2.
         Vector2 tangent = normalize(v2);
         const Vector3 u = curElem->interpolateVelocities(xt, this);

         const Vector2 u2D = make_double2(u.x, u.y);
         u_int = toVector3(dot(u2D, tangent)*tangent, u.z);


      } else {

         assert(curElem!=nullptr && !curElem->isDry && "Final element is dry or boundary side");

         time_rm -= dtb; 

         u_int = curElem->interpolateVelocities(xt, this);

      }

      if (time_rm <= 0.0) break; 

      // If the velocity is too small, we stop the backtracking
      const double vmag2 = norm2(toVector2(u_int));
      if (vmag2 <= velmin_btrack_sq)  {
         time_rm = 0.0; 
         return; 
      }

      // The current point becomes the next segment's starting point.
      x0 = xt;
 
      iter++;

// XXXXXXXXXXXXXXXXX
if (iter == 38) {printf("*****iter = %d, time_rm = %.17f\n", iter, time_rm); return;}
      assert(iter < MAX_ITERATIONS && "Max iterations reached in btrack");

      assert(time_rm_3 > time_rm && "time_rm is not decreasing");
      time_rm_3 = time_rm_2;
      time_rm_2 = time_rm_1;
      time_rm_1 = time_rm;
   }

   assert(curElem != nullptr && "Final element is nullptr");
   assert(!curElem->isDry && "Final element is dry");
   


   // Interpolate the velocities at the final point with respect to the current side frame
   u_int = curElem->interpolateVelocities(xt, this);


}




__device__ double Node::sigmaTransitionDepth(void) const { return min(depth, mesh->h_s); }  

__device__ double Node::computeETP(int idxComp, double time, double* tfreq, double* tamp, double* tnf, double* tear, int* jspc) const {

   assert(lat >= -pi/2.0 && lat <= pi/2.0 && "computeETP: lat out of range");
   assert(lon >= -pi && lon <= pi         && "computeETP: lon out of range");

   const int ncyc = int(tfreq[idxComp]*time/(2.0*pi));
   const double arg = tfreq[idxComp]*time - ncyc*(2.0*pi) + jspc[idxComp]*lon + tear[idxComp];

   double fun_lat = 0.0;
   switch (jspc[idxComp]) {
      case 0:
         fun_lat = 3.0*pow<2>(sin(lat)) - 1.0;
         break;
      case 1:
         fun_lat = sin(2.0*lat);
         break;
      case 2:
         fun_lat = pow<2>(cos(lat));
         break;
      default:
         assert(false && "computeETP: invalid jspc");
   }
   
   return 0.69*ramp*tamp[idxComp]*tnf[idxComp]*fun_lat*cos(arg);

}

__device__ void Node::updateKb(void) {
// Updates the node's kb index. Must be called before compZCoord(), preferibly right after Depth and/or eta2 are updated

   if constexpr(ivcor == 1) return; // When ivcor=1, the kb index is fixed and read from the vgrid file  XXXX What if the bed changes???


   // assert(!isDry && "ZCOOR: node is dry"); XXXX

   if constexpr(ivcor == 2) {  // SZ coordinates
      
      // Compute bottom level index
      if (depth <= mesh->h_s) {
         kb = izlev;   
         
      } else {
         kb = -1;
         int L = 0;
         int R = nzlev; 
         
         while (L <= R) {
            int M = (R + L) / 2;


            assert((-depth >= mesh->ztot[L] && -depth < mesh->ztot[R]) && "ZCOOR: invalid depth\n");
            
            if (-depth >= mesh->ztot[L] && -depth < mesh->ztot[M]) {
               R = M;
            } else {
               L = M;
            }

            if (R - L == 1) {
               kb = L;
               break;
            }
         }

         // Checks that the binary search was successful
         assert(-depth >= mesh->ztot[kb] && -depth < mesh->ztot[kb + 1]);
       
      }
   } 
}


__device__ __inline__ void Node::updateZCoordData(void) {
// XXX: To be called before compZCoord(), preferibly right after Depth and/or eta2 are updated

   // assert(!isDry2() && "updateZCoordData: node is dry");
   // if (isDry2()) {
   //    kb = -1;  // XXXXXXXX KKKKK
   //    return;
   // }


   if constexpr(ivcor == 2) {  // SZ coordinates
      
      // Compute bottom level index
      if (!isDry) updateKb();  // Why is this call repeated? XXXXXXXX It is also called in ul0UpdateZCoords
   } 
   else {  // if (ivcor == 1) Localized sigma
      // Do nothing
   }


}


__device__ __forceinline__ double Node::compZCoord(int k) const {
   // Computes the z coordinate at level k.
   // Depending on the vertical coordinate system, the z coordinate is computed differently.
   // When ivcor=1, the z coordinate is computed using the localized sigma levels.
   // When ivcor=2, the z coordinate is computed using the zeta/sigma levels. Z levels are
   // used 

   assert(k >= 0 && k < nvrt && "compZCoord: invalid level\n");
   
   // NOTE: Assumes that updateZCoordData() has been called before this function.

   if (isDry2()) return 0.0;  

   assert(htot2() > h0 && "compZCoord: invalid depth\n");

   if constexpr(ivcor == 2) {

      const double hmod2 = sigmaTransitionDepth(); // Depth transition between z and sigma levels
      const double h_c = mesh->h_c;
      const double theta_f = mesh->theta_f;
      const double sinh_theta_f = mesh->sinh_theta_f;


      // TODO: Here we are assuming that the bottom level, kbpl, was found in a previous step.
      // The FORTAN code computes it every time this function is called, but it should be done only 
      // once per time step/node. XXXX
      // CHECK: updateZCoordData

      //  To simplify comparisons, I define index versions of count variables.
      if      (k == ivrt)     return eta2; 
      else if (k == izlev)    return -hmod2;
      else if (k == kb)       return -depth;
      else if (k < izlev)     return mesh->ztot[k];
      else if (k > izlev) {
         const int ksig = k - izlev;

         if (hmod2 <= h_c) {
            return mesh->sigma[ksig] * (hmod2 + eta2) + eta2;

         } 
         else if (eta2 <= - h_c - (depth - h_c)*theta_f/sinh_theta_f) {
            printf("Please choose larger h_c: %f %f\n", eta2, h_c);
            assert(false && "compZCoord: invalid h_c\n");
         } 
         else {
            return eta2*(1.0 + mesh->sigma[ksig]) + 
                   h_c * mesh->sigma[ksig] + 
                   (hmod2 - h_c) * mesh->cs[ksig];
         }
      } else { 
         assert(false && "compZCoord: invalid level or coordinate system\n");
      }

   }
   else {  // if constexpr(ivcor == 1)
      if (k <= kb)        return -depth; 
      else if (k == ivrt) return eta2;
      else if (k > kb)    return (eta2 + depth)*localSigma[k] + eta2;

   }
   
   assert(false && "compZCoord: invalid level or coordinate system\n");
   return 0.0;
}


__device__ __forceinline__ void Node::computeBottomDragCd(void) {
   // Computes the bottom drag Cd based on the selected nchi

   if constexpr(nchi == 0) {
      // This line is repeated over and over again every time step, but it is very short, so it is not worth optimizing.
      Cd = drag;  // Read from drag file

   } else if constexpr(nchi == 1) {
      // 

      bool ltmp = false; // for WBL iteration

      assert(rough >= 0.0 && "computeBottomDragVegetation: rough < 0");
      if (isDry || (rough == 0.0)) {
         Cd = 0.0;

      } else { //roughness >0 
         
         const double bthick = max(dzb_min, z[kb + 1] - z[kb]); // thickness of bottom bnd layer 

         assert(bthick > rough && "computeBottomDragVegetation: bthick <= rough");

         Cd = 1.0/pow<2>((2.5*log(bthick/rough)));       

      } // if rough

   } else if constexpr(nchi == -1) {
   // Manning 
 
      const double htot = max(hmin_man, htot2()); 
      Cd = isDry ? 0.0 : grav2*rmanning*rmanning/pow(htot, 0.333);

   } else {
      assert(false && "computeBottomDragVegetation:selected nchi not implemented.");
      
   } // if nchi

}


__device__ int findSideBoundaryType(Side &side, int *iond, int *nond, int *ilnd, int *nlnd, bool &isLand) {
   // Checks if a side in on the boundary
   // For this it checks if it belongs to any of the boundary segments. They can be either land or open boundary segments.

   // Returns -1 if the side is not on the boundary, otherwise returns the index of the boundary segment.
   // isLand is set to true if the side is on a land boundary, otherwise it is set to false.

   // Indices of the nodes that form the side
   //
   // NOTE: This function is slow. It's to be called from "populate" functions, not from the main loop.
   const int n1 = side.n1->idx;
   const int n2 = side.n2->idx;
   isLand = false;
   

   // Check if the side is on a land boundary
   for (int i = 0; i < nland; i++) {
      for (int j = 0; j < nlnd[i] - 1; j++) {
         // Consecutive nodes in the land boundary segment
         const int b1 = ilnd[i + j*nland];
         const int b2 = ilnd[i + (j + 1)*nland];

         assert(b1 != -1 && b2 != -1);

         if ((b1 == n1 && b2 == n2) || (b1 == n2 && b2 == n1)) {
            
            isLand = true;
            return i;
         }
      }
   }

   // Check if the side is on an open boundary
   for (int i = 0; i < nope; i++) {
      for (int j = 0; j < nond[i] - 1; j++) {
         // Consecutive nodes in the open boundary segment
         const int b1 = iond[i + j*nope];
         const int b2 = iond[i + (j + 1)*nope];

         assert(b1 != -1 && b2 != -1);

         if ((b1 == n1 && b2 == n2) || (b1 == n2 && b2 == n1)) {
            isLand = false;
            return i;
         }
      }
   }

   // Check that the nodes are not both on the boundary
   assert( (side.n1->nodeBnd == nullptr || side.n2->nodeBnd == nullptr)
        && "side not in ilnd/iond but both nodes are boundary" );

   return -1;
}


__device__ void Mesh::populateNodes(int *i34, int *isbnd, int *sharedIdxNodeBoundary,
                                    double *node_coords, double *node_coords_ll, bool use_ll,
                                    double *wind_rotate_angle, 
                                    int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                                    int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                                    int16_t *nodePermutationSort, 
                                    int *localNodeNeighNodesStart, int *localNodeNeighElemsStart, 
                                    int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                                    int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                                    int *localNodePart,     int *localElemPart,     int *localSidePart,
                                    double *crs_data, 
                                    double *drag, double *roughness, double *rmanning, int *kb, double *localSigma) {

   NODE(pnode);  
   const int i = idx_node;
   const int iglobal = globalNodeIndices[i];
   Node &node(*pnode);

   node.idx = iglobal;
   node.part = localNodePart[iglobal];

   // It is important to differentiate between the global and local "start" indices because they
   // play different, but related, roles in this code. The global "start" indices are used
   // for the CSR-like indexing of the node and element neighbor arrays in the "global"
   // (unpartitioned) mesh. The local "start" indices are used for the CSR-like indexing of the 
   // node and element neighbor arrays in the "local" (partitioned) mesh.
   const int idxNodeStart = nodeNeighNodesStart[iglobal];
   const int idxNodeEnd   = nodeNeighNodesStart[iglobal + 1];
   const int idxNodeLocalStart = localNodeNeighNodesStart[i];
   node.numNodes = idxNodeEnd - idxNodeStart;
   node.nodes = nodeNeighNodes + idxNodeLocalStart;

   const int rowStart = idxNodeLocalStart + i; // +i is to account for the diagonal element.
   node.coeffRow = crs_data + rowStart; 



   const int idxElemStart = nodeNeighElemsStart[iglobal];
   const int idxElemEnd   = nodeNeighElemsStart[iglobal + 1]; 
   const int idxElemLocalStart = localNodeNeighElemsStart[i];
   node.numElements = idxElemEnd - idxElemStart;
   node.elements = nodeNeighElems + idxElemLocalStart;


   for (int idx = idxNodeStart, idxLocal = idxNodeLocalStart; idx < idxNodeEnd; idx++, idxLocal++) {
      const int idxNodeGlobal = nodeNeighNodesIndices[idx];
      const int partNode = localNodePart[idxNodeGlobal]; 
      const int idxNode  = localNodeIndices[idxNodeGlobal];

      Node *node_ptr = nodes[partNode] + idxNode;
      
      nodeNeighNodes[idxLocal] = node_ptr;
   }


   for (int idx = idxElemStart, idxLocal = idxElemLocalStart; idx < idxElemEnd; idx++, idxLocal++) {
      const int idxElemGlobal = nodeNeighElemsIndices[idx];
      const int partElem = localElemPart[idxElemGlobal];
      const int idxElem  = localElemIndices[idxElemGlobal];

      Element *elem_ptr = elems[partElem] + idxElem;

      nodeNeighElems[idxLocal] = elem_ptr;
   }      

   // Compute the cosine and sine of the mesh rotation angle. We use float here because memory efficiency is more important than precision in this case.
   assert(wind_rotate_angle[i] >= -pi && wind_rotate_angle[i] <= pi && "wind_rotate_angle out of range");

   sincosf(-(float)wind_rotate_angle[i], &node.sRot, &node.cRot);

   // Fills the node coordinates. If ics == 1, the coordinates are in the local frame. If ics == 2, the coordinates are in the global frame.
   if constexpr(ics == 1) {
      node.p.x   = node_coords[i * 3];
      node.p.y   = node_coords[i * 3 + 1];
      node.gz = 0.0;

      if (use_ll) {                                    // <-- add
         node.lon = node_coords_ll[i * 2    ]*deg2rad; // <-- add
         node.lat = node_coords_ll[i * 2 + 1]*deg2rad; // <-- add
      }  

   } else {
      const double lon = deg2rad*node_coords[i * 3];
      const double lat = deg2rad*node_coords[i * 3 + 1];

      Point3 p = computeXYZfromLatLon(lat, lon);
      node.p = toVector2(p);
      node.gz = p.z;
      node.lon = lon;
      node.lat = lat;

      #ifdef USE_FRAMES
      node.setFrame(lat, lon);
      #endif


   }

   node.depth = node_coords[i * 3 + 2];

   node.numSides = 0;

   // drag, rough, and manning are actually stored in the same memory location. This is just sugar to make the code more readable without using unnecessary memory.
   if      constexpr(nchi ==  0) node.drag     = drag[iglobal];
   else if constexpr(nchi ==  1) node.rough    = roughness[iglobal];
   else if constexpr(nchi == -1) node.rmanning = rmanning[iglobal];
   else assert(false && "populateNodes: selected nchi not implemented.");

   if constexpr(ivcor == 1) 
   {  
      node.kb = kb[iglobal];

      for (int k = 0; k < nvrt; k++) node.localSigma[k] = localSigma[k + iglobal*nvrt];
   }

   for (int k = 0; k < nvrt; k++) node.rho[k] = rho0;
   for (int k = 0; k < nvrt; k++) node.z[k] = 0.0;
   for (int k = 0; k < nvrt; k++) node.dfv[k] = dfv0;   // XXXX 

   
   // Links the node to the boundary information
   if (isbnd[iglobal] >= 0) {
     
      // Get the index of the node boundary
      const int idxNodeBoundary = atomicAdd(sharedIdxNodeBoundary, 1) + 1;

      node.nodeBnd = &(localNodeBoundaryMem[idxNodeBoundary]);
      node.nodeBnd->node = &node;

      node.nodeBnd->sideBnd1 = nullptr;
      node.nodeBnd->sideBnd2 = nullptr;

      if (isbnd[iglobal]>=0) {

         node.nodeBnd->segment = &(localSegmentMem[isbnd[iglobal]]);
         
         node.nodeBnd->segment->ifltype = 0; //(isbnd[i]==0)?3:0; // XXXXX KKKKK
         node.nodeBnd->segment->iettype = 0; // XXXXXXXXX

         // Initialize the length and area of the segment to nan to indicate that they are not yet computed
         node.nodeBnd->segment->clen  = nan("1");
         node.nodeBnd->segment->carea = nan("1");

      }
   } else {

      node.nodeBnd = nullptr; // For internal nodes

      if (isbnd[iglobal] == landSegment) {
         node.nodeBnd = landBoundaryNode;
         node.nodeBnd->segment = landBoundarySegment;
      }
      
   }
   
   // XXXX GAUSSIAN HILL VELOCITIES
   #ifdef GAUSSIAN_HILL_VELOCITIES
         if (i == 0) printf("GAUSSIAN HILL VELOCITIES!!!\n");
         constexpr double rot_per = 3000.0;
         constexpr double rot_f = 0.0*2.0*pi/rot_per;
         for (int k = 0; k < nvrt; k++) { node.u2[k] = {-node.p.y*rot_f, node.p.x*rot_f}; node.w2[k] = 0.0; }
      // XXXXX
   #endif

}

__device__ void Mesh::populateNodes2(int *localNodeNeighNodesStart, int *crs_indices, int *crs_indptr) {
   NODE(node);


   // Computes csrSlot, the array of indices/slots of the neighbor nodes in the CRS matrix.
   int16_t *csrSlot = node->csrSlot;  

   // Simply count the number of nodes with a smaller global index than the current node.
   // It accounts for the diagonal element.
   for (int i = 0; i < node->numNodes; ++i) {
      const int col = node->nodes[i]->idx;
      int16_t idx = 0;

      //Takes care of the diagonal element.
      if (node->idx < col) ++idx;          

      for (int j = 0; j < node->numNodes; ++j) if (node->nodes[j]->idx < col) ++idx;
      csrSlot[i] = idx;
   }

   // Find the diagonal element in the CRS matrix.
   int16_t idxDiag = 0;
   for (int j = 0; j < node->numNodes; ++j)
      if (node->nodes[j]->idx < node->idx) ++idxDiag;

   node->coeffDiag = &(node->coeffRow[idxDiag]);

   // Fills the CRS indices and indptr arrays for the cupyz's crs matrix in the solver.
   const int idxNodeLocalStart = localNodeNeighNodesStart[i];

   const int rowStart = idxNodeLocalStart + i;  // +i is to account for the diagonal element.
   const int rowNNZ = node.numNodes + 1;        // +1 is to account for the diagonal element.
   crs_indptr[i] = rowStart;

   const int d = (int)(node->coeffDiag - node->coeffRow);
   crs_indices[rowStart + d] = node->idx;  
   for (int j = 0; j < node->numNodes; ++j) crs_indices[rowStart + node->csrSlot[j]] = node->nodes[j]->idx;


}




__device__ void Mesh::populateSides(int *i34, int *sharedIdxSideBoundary,
                                    double *node_coords, double *node_coords_ll, bool use_ll,
                                    double *wind_rotate_angle,   
                                    int *isdel, int *isidenode, int *isbnd,
                                    int *nond, int *nlnd, int *iond, int *ilnd,                                  
                                    int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                                    int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                                    int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                                    int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                                    int *localNodePart,     int *localElemPart,     int *localSidePart,
                                    double *drag, double *roughness) {


   SIDE(pside);
   Side &side(*pside);
   const int i = idx_side;
   const int iglobal = globalSideIndices[i];
   side.idx = iglobal;  
   
   // Parameters used to compute a linear approximation of the Coriolis coefficient around the reference latitude sfea0
   const double fc   = 2.0*omega_e*sin(sfea0);
   const double beta = 2.0*omega_e*cos(sfea0);

   assert(sfea0 >= -pi/2.0 && sfea0 <= pi/2.0 && "GRID: sfea0 out of range");

   // finding the pointer to the nodes that form the side
   const int idx_n1 = isidenode[iglobal*2];
   const int idx_n2 = isidenode[iglobal*2 + 1];

if (i % 100 == 0) printf("populateSides: iglobal = %d, idx_n1 = %d, idx_n2 = %d\n", iglobal, idx_n1, idx_n2);

   const int locidx_n1 = localNodeIndices[idx_n1];
   const int locidx_n2 = localNodeIndices[idx_n2];
   const int part_n1 =   localNodePart[idx_n1];
   const int part_n2 =   localNodePart[idx_n2];

if (i % 100 == 0) printf("populateSides: i = %d, part_n1 = %d, part_n2 = %d, locidx_n1 = %d, locidx_n2 = %d\n", i, part_n1, part_n2, locidx_n1, locidx_n2);
if (i % 100 == 0) printf("populateSides: i = %d, nodes[part_n1] = %p, nodes[part_n2] = %p\n", i, mesh->nodes[part_n1], mesh->nodes[part_n2]);

   side.n1 = (Node *)&(mesh->nodes[part_n1][locidx_n1]);
   side.n2 = (Node *)&(mesh->nodes[part_n2][locidx_n2]);

if (i % 100 == 0) printf("populateSides: 1, i = %d, side.n1 = %p, side.n2 = %p\n", i, side.n1, side.n2);

   // Counts the number of sides connected to each node
   atomicAddInt16Unsafe(&side.n1->numSides, 1);  // Unsafe here means that it doesn't check for overflow. This function is defined in this file.
   atomicAddInt16Unsafe(&side.n2->numSides, 1);

   // Side boundary type
   NodeBoundary *n1bnd = side.n1->nodeBnd;
   NodeBoundary *n2bnd = side.n2->nodeBnd;
if (i == 0) printf("populateSides: 2\n");

   if (n1bnd != nullptr && n2bnd != nullptr) 
   // If both nodes are on the boundary, then the side may be on the boundary (Notice that there are cases where the side is not on the boundary, but both nodes are on the boundary)
   {
      bool isLand = false;
      int isbs = findSideBoundaryType(side, iond, nond, ilnd, nlnd, isLand);
      
      if (isbs > -1) {
         // assert(isLand || isbs>=0);
         const int idxSideBoundary = atomicAdd(sharedIdxSideBoundary, 1) + 1;

         side.sideBnd = &(localSideBoundaryMem[idxSideBoundary]);

         if   (isLand) side.sideBnd->segment = landBoundarySegment;
         else          side.sideBnd->segment = &(localSegmentMem[isbs]);

         assert (isLand == side.sideBnd->segment->isLand);
         
         side.sideBnd->side = &side;

         // Links each node boundary to the side boundary using the first available side boundary pointer
         n1bnd->linkSideBoundary(side.sideBnd);
         n2bnd->linkSideBoundary(side.sideBnd);

      }
      else {
         side.sideBnd = nullptr;
      }
   }
   else 
   {
      side.sideBnd = nullptr;
   }

if (i % 100 == 0) printf("populateSides: 3, i = %d\n", i);

   // finding the pointer to the elements that share the side
   const int idx_e1 = isdel[iglobal*2];
   const int idx_e2 = isdel[iglobal*2 + 1];

   if (idx_e1 > -1) {
      const int locidx_e1 = localElemIndices[idx_e1];
      const int part_e1 =  localElemPart[idx_e1];
      side.elems[0] = elems[part_e1] + locidx_e1;

   }
   else {
      side.elems[0] = nullptr;
   }

if (i % 100 == 0) printf("populateSides: 4, i = %d\n", i);

   if (idx_e2 > -1) {
      const int locidx_e2 = localElemIndices[idx_e2];
      const int part_e2 =  localElemPart[idx_e2];
      side.elems[1] = elems[part_e2] + locidx_e2;
   }
   else {
      side.elems[1] = nullptr;
   }

if (i % 100 == 0) printf("populateSides: 5, i = %d\n", i);

   if (side.elems[1] == nullptr && side.sideBnd == nullptr) {
      
      NodeBoundary *n1bnd = side.n1->nodeBnd;
      NodeBoundary *n2bnd = side.n2->nodeBnd;
      printf("Side %d has no elements and no boundary. Types = %p, %p\n", i, n1bnd, n2bnd);
         
      // printf("idx_e1: %d, idx_e2: %d\n", idx_e1, idx_e2);
      // printf("elems[0]: %p, elems[1]: %p\n", side.elems[0], side.elems[1]);
      // printf("sideBnd: %p\n", side.sideBnd);
      // assert(false && "Side has no elements and no boundary");
   }
   // assert(side.elems[1] != nullptr || side.sideBnd != nullptr);

   // Compute the center of the side
   side.center = 0.5*(side.n1->p + side.n2->p);

   #ifdef USE_FRAMES
   side.setFrame();
   #endif

if (i % 100 == 0) printf("populateSides: 6, i = %d\n", i);

   // Normal, unitary vector to the side
   if constexpr(ics == 1) {
      side.n = { -(side.n1->p.y - side.n2->p.y), (side.n1->p.x - side.n2->p.x)};
      side.n = normalize(side.n);
   }
   else {
      // Unit side vector in 3D cartesian coordinates
      Vector3 nSide = { side.n1->p.x - side.n2->p.x, side.n1->p.y - side.n2->p.y, side.n1->gz - side.n2->gz };

      nSide = normalize(nSide);

      // Its cross product with the unit local z-axis vector (nUp) is the
      // normal vector we are looking for.
      Vector3 n3D = normalize(cross(nSide, side.frame.nUp()));  

      // The 2D components can be obtained by projecting the 3D normal vector to the plane defined by the unit vectors nZonal and nMeridional
      const Vector2 n2D = {dot(n3D, side.frame.nZonal()), dot(n3D, side.frame.nMeridional())};
      side.n = -normalize(n2D); 

   }

if (i % 100 == 0) printf("populateSides: 7, i = %d\n", i);

   assert(ishapiro == 1 && "Only shapiro = 1 is implemented");
   side.shapiro = shapiro0;

   #ifdef GAUSSIAN_HILL_VELOCITIES
   // XXXX
   if (i == 0) printf("GAUSSIAN HILL VELOCITIES!!!\n");
   constexpr double rot_per=3000.0;
   constexpr double rot_f=0.0*2.0*pi/rot_per;/// XXXX
   for (int k = 0; k < nvrt; k++) side.u2[k] = {-side.center.y*rot_f , side.center.x*rot_f};
   // XXXXX      
   #endif
   // Compute the z level at the center of the side
   for (int k = 0; k < nvrt; k++) side.zc[k] = 0.5*(side.n1->z[k] + side.n2->z[k]);  // XXXX substitute by computation?
   
   // Compute the length of the side
   side.len = length(toPoint3(side.n2->p, side.n2->gz) - toPoint3(side.n1->p, side.n1->gz));

if (i % 100 == 0) printf("populateSides: 8, i = %d\n", i);

   // Compute the Coriolis coefficient

   if constexpr(ncor == -1 || ncor == 0) {
      side.coriolis = coricoef;
   }
   else {
      const double sphi = side.lat();  // Latitude at the center of the side

      assert(use_ll && "use_ll must be true when ncor > 0");

      if constexpr(ics == 1) {
         side.coriolis = fc + beta*(sphi - sfea0);  // First order Taylor expansion of the Coriolis coefficient around the reference latitude sfea0
      }
      else { // ics == 2
         assert(sphi >= -pi/2.0 && sphi <= pi/2.0 && "GRID: sphi out of range");
         side.coriolis = 2.0*omega_e*sin(sphi);  
      }

   }

if (i % 100 == 0) printf("populateSides: 9, i = %d\n", i);

}


__device__ void Mesh::populateElements(int *i34, int *elnode, int *elside,
   double *node_coords, double *node_coords_ll, bool use_ll,
   double *wind_rotate_angle,   
   int *isdel, int *isidenode, int *isbnd,
   int *nond, int *nlnd, int *iond, int *ilnd,                                  
   int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
   int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
   int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
   int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
   int *localNodePart,     int *localElemPart,     int *localSidePart,
   double *drag, double *roughness) {

   ELEM(pelem);
   Element &elem(*pelem);
   const int i = idx_elem;

   elem.idx = i;

   for (int j = 0; j < i34[i]; j++) {
      const int idxNode = elnode[i*4 + j];  
      const int idxSide = elside[i*4 + j];
      
      elem.nodes[j] = (Node *)&(nodes[localNodePart[idxNode]][localNodeIndices[idxNode]]);
      elem.sides[j] = (Side *)&(sides[localSidePart[idxSide]][localSideIndices[idxSide]]);

      // Find the neighboring element  
      if (elem.sides[j]->elems[0] == &elem) {
         elem.nghbElems[j] = elem.sides[j]->elems[1];
      }
      else {
         elem.nghbElems[j] = elem.sides[j]->elems[0];
      }
      
   }

   elem.i34 = i34[i];

   // Compute the barycenter of the element and its depth
   elem.center = {0.0, 0.0};
   elem.depth = 1e10;
   for (int j = 0; j < i34[i]; j++) {

      elem.center += elem.nodes[j]->p/i34[i];
      elem.depth  = min(elem.depth, elem.nodes[j]->depth);
   }

   #ifdef USE_FRAMES
   elem.setFrame();
   #endif

   elem.area = elem.computeSignedArea(); 
   elem.radius = sqrt(elem.area)/pi; 
   
   // Determines which way the normal vector to each side is oriented.
   // The orientation means that for each side, elems[0] points to elems[1]. In the boundary this guarantees that the normal vector is pointing outwards.
   for (int j = 0; j < i34[i]; j++)  
   {
      if (elem.sides[j]->elems[0] == &elem) {
         elem.ssign[j] =  1.0;
      }
      else {
         elem.ssign[j] = -1.0;
      }

   }

  

   // Derivatives of shape functions
   //  For ics=2, this is done in the element frame
   //  For quads, the derivative is evaluated at centroid
   for (int j = 0; j < elem.i34; j++) {
      const double area = elem.area;
      const int id2 = elem.nxq(1 + j);
      const int id3 = elem.nxq(2 + j);
      const int id4 = elem.nxq(3 + j);
      if (elem.isTri()) {
         // This expression can be derived from the derivatives of the shape functions in area coordinates.
         // phi_0 = (x2 - x0)^(x - x0)/((x2 - x0)^(x1 - x0)), where ^ is the cross product
         // phi_0 = (x2 - x0)^(x - x0)/(2*area)
         // dphi_0/dx = (x2.y - x0.y)/(2*area)
         // dphi_0/dy = (x0.x - x2.x)/(2*area)
         // And so on for the other shape functions.

         // Get the horizontal coordinates of the nodes id2 and id3 in the local frame
         const Point2 p2 = elem.getLocalHorizontalNodePosition(id2);
         const Point2 p3 = elem.getLocalHorizontalNodePosition(id3);

         const double dPhi_j_dx = (p2.y - p3.y)/(2*area); 
         const double dPhi_j_dy = (p3.x - p2.x)/(2*area); 

         elem.gradPhi[j] = {dPhi_j_dx, dPhi_j_dy};

      }
      else  { // Quad, the derivatives are evaluated at the centroid

         // The derivation is complicated, but it starts by defining the shape functions in terms of the local coordinates (eta, xi).
         // phi_0 = (1 + eta*eta_n(0))*(1 + xi*xi_n(0))/4
         // ...
         // and then the derivatives are computed at the centroid (eta = 0, xi = 0)

         // Get the horizontal coordinates of the nodes id2 and id4 in the local frame
         const Point2 p2 = elem.getLocalHorizontalNodePosition(id2);
         const Point2 p4 = elem.getLocalHorizontalNodePosition(id4);

         const double dPhi_j_dx = (p2.y - p4.y)/(2*area); 
         const double dPhi_j_dy = (p4.x - p2.x)/(2*area); 

         elem.gradPhi[j] = {dPhi_j_dx, dPhi_j_dy};

      } // if Tri/Quad

   } // for j   
   
   

}


__device__ void Mesh::populateBoundaryNodes(int *nond, int *iond,
                                            double *emo, double *efa, 
                                            double *umo, double *ufa, 
                                            double *vmo, double *vfa,
                                            int *ifltype, int *iettype) {

   NODE_BND(nodeBnd);
   const int idx = idx_node_bnd;    

   nodeBnd->idx = idx;

   // Need to find the segment index and the index of the node on the segment
   int idxNode = -1;
   int iseg, i;
   for (iseg = 0; iseg < nope; iseg++) {
      // for each segment...
      for (i = 0; i < nond[iseg]; i++) {
         idxNode = iond[iseg + i*nope];
         if (idxNode == nodeBnd->node->idx) goto found; // Cleanest way to exit the double loop
         
      }
   }

   assert(iseg < nope && "iseg is out of bounds");

found: // Jumped here when the node was found

   nodeBnd->isEssentialBC = true;
   nodeBnd->elbc = 0.0;

   nodeBnd->segment->ifltype = ifltype[iseg];
   nodeBnd->segment->iettype = iettype[iseg];

   for (int jfr = 0; jfr < nbfr; jfr++) {

      const int idx2 = iseg + i*nope + jfr*nope*mnond; 

      nodeBnd->efa[jfr] = efa[idx2];
      nodeBnd->emo[jfr] = emo[idx2];
      nodeBnd->umo[jfr] = {umo[idx2], vmo[idx2]};
      nodeBnd->ufa[jfr] = {ufa[idx2], vfa[idx2]};

   }
        

}

__device__ void Mesh::populateBoundarySideLinks(void) {

   SIDE_BND(sideBnd);

   const Side *side = sideBnd->side;

   // We assume that the leftSide is the side s that has s->n2 == side->n1
   // and the rightSide is the side s that has s->n1 == side->n2
   for (int i = 0; i < mesh->nbs; i++) {

      Side *s = &mesh->localSideData[i];

      if (s == side) continue;  // Do not link the side to itself
   
      if (s->n1 == side->n2) sideBnd->leftSide  = s;
      if (s->n2 == side->n1) sideBnd->rightSide = s;
   }
   
}

__device__ void Mesh::populateMesh(int *i34, int *elnode, int *elside,
                double *node_coords, double *node_coords_ll, bool use_ll,
                double *wind_rotate_angle, int *sharedIdxNodeBoundary, int *sharedIdxSideBoundary,
                int *isdel, int *isidenode, int *isbnd,
                int *nond, int *nlnd, int *iond, int *ilnd, 
                int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                int *localNodePart,     int *localElemPart,     int *localSidePart,
                double *drag, double *roughness,
                double *emo, double *efa, 
                double *umo, double *ufa, 
                double *vmo, double *vfa) {

   printf("Populating mesh for partition %d\n", partNum);

   // XXXXXX THIS function does almost nothing. 

   nbs = *sharedIdxSideBoundary + 1;
   nbp = *sharedIdxNodeBoundary + 1;

   nEntries = (nodeNeighNodesStart[nbp] - nodeNeighNodesStart[0]) + nbp;

}


extern "C" __global__ void setupCoefficientMatrixAndRHS(double time) {
   // Setup coefficient matrix, sparsem, for the shallow water wave equation
   // No elevation essential b.c. are imposed yet but other b.c. is imposed


   NODE(node1);

   node1->clearMatrixRow();
   node1->qel = 0.0;

   // THIS IS FOR STRUCTURES XXX
   // if ((ihydraulics!=0 && nhtblocks>0) && isblock_el(ie)>0) cycle //active block

   
   for (int j=0; j<node1->numElements; j++) {
      // For each element connected to the node...

      Point2 x0, x1, x2, x3; // Local horizontal node positions

      Element *elem = node1->elements[j];

      // Compute the average gravity
      const double grav3 = elem->mean(&Node::grav2);

      // Look at the documentation for the description of the terms in the following integrals

      // Integral I_1 
      // \bar{\breve{H}} 
      const double hhatb = elem->mean(&Side::hhat);        

      assert((hhatb >= 0.0 || ihhat != 1) && "Invalid hhat");

      const int id1 = elem->nodeLocalIdx(node1);  // Index of node1 in element->nodes[]
      const double area = elem->area;

      // Integral I_4
      // I_4 = \int_\Omega \left[ \phi_i\eta^n + \theta\Delta t \nabla\phi_i\cdot\pmb{E} + (1-\theta)\Delta t \nabla\phi_i \cdot\pmb{U}^n \right]d\Omega
      // Where
      // U is pmb{U}^{n+1} = \pmb{E} - g\theta\breve{H}\Delta t\nabla\eta^{n+1}

      // Since \nabla\phi_i is constant over the element, for the integral of U we can compute the mean value
      // of U over the element.
      const Vector2 U = elem->mean(&Side::U);

      const Vector2 gradPhi_1 = elem->gradPhi[id1];

      // Integral of \nabla\phi_i \cdot \pmb{U}^n
      const double I_4_1 = dot(gradPhi_1, U);


      // Integral of \nabla\phi_i \cdot \pmb{E}
      const double I_4_2 = dot(gradPhi_1, elem->ghat1);

      node1->qel += dt*area*((1 - thetai)*I_4_1 + thetai*I_4_2);      


      // Normal component of the bed velocity
      double ubedNormalComp = 0.0;
      if constexpr(imm == 2)  { // Pre-compute vnorm on the element
         Vector3 ubed = {0.0, 0.0, 0.0};
         // update_bdef(time, elem->center, elem->depth, ubed); XXX
         assert(false && "computeVerticalVelocity: not implemented for imm==2");

         // Computes a unit vector perpendicular to the bottom face. It is computed as the unit vector parallel to (grad(h), 1)
         const Vector2 gradDepth = elem->grad2(&Node::depth);
         const Vector3 N = normalize(toVector3(gradDepth, 1.0));

         // Computes the normal component of the bed velocity
         ubedNormalComp = dot(ubed, N);
      } // imm==2

      const double fricCoeff = grav3*thetai*(1.0 - thetai)*dt*dt*hhatb;

      if (elem->isTri()) {

         // To clarify the notation in this section: 
         //   In formulas, we use the following notation:
         //       i stands for the node under consideration (node1)
         //       j (node2) and k (node3) are the other nodes of the triangle.
         //   However, in the code, we use the following notation:
         //       id1 is the node under consideration (node1), id2, id3 are the local indices of the other nodes (node2, node3)
         const int id2 = elem->nxq(1 + id1);
         const int id3 = elem->nxq(2 + id1);
         Node *node2 = elem->nodes[id2];
         Node *node3 = elem->nodes[id3];

         elem->getLocalHorizontalNodePositionsTri(x0, x1, x2, id1);

         // integral of phi_i*phi_j
         const double ipp_ii = (1.0/6.0 )*area;
         const double ipp_ij = (1.0/12.0)*area;
         const double ipp_ik = (1.0/12.0)*area;

         // integral of [(grad phi_i) * (grad phi_j)]
         const double igpgp_ii = ((1.0/4.0)/area)*dot(x2-x1, x2-x1);
         const double igpgp_ij = ((1.0/4.0)/area)*dot(x2-x1, x0-x2);
         const double igpgp_ik = - igpgp_ii - igpgp_ij; // This comes from the fact that if a,b,c are vectors forming a triangle (in order), then a+b+c = 0, hence: -a.c = a.(a+b+c-c) = a.a + a.b

         // Computes the different components of integral I_1 of each node with respect to node i.
         const double I1_ii = ipp_ii + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ii;
         const double I1_ij = ipp_ij + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ij;
         const double I1_ik = ipp_ik + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ik;

         node1->addToSparseMatrix(node1, I1_ii);
         node1->addToSparseMatrix(node2, I1_ij); 
         node1->addToSparseMatrix(node3, I1_ik);

         #pragma unroll
         for (int l=0; l<3; l++) {

            Node *node = elem->nodes[l];

            // Additional variation in eta due to bed deformation or vertical velocity
            double extraEta = 0.0;
            if      constexpr(imm == 2) extraEta = dt*ubedNormalComp;
            else if constexpr(imm == 1) extraEta = node->bdef2 - node->bdef1;


            const double ipp   = elem->intPhiPhiTri(id1, l);
            const double igpgp = elem->intGradPhiGradPhiTri(id1, l);

            const double isElemWet = elem->isDry ? 0.0 : 1.0;  // XXXXX I added this to replicate Fortran's code, but Joseph is going to check if it is needed.

            node1->qel += (node->eta2 + extraEta)*ipp - node->eta2*fricCoeff*igpgp*isElemWet;
           
         }


      }
      else { // Quad
         // To clarify the notation in this section: 
         //   In formulas, we use the following notation:
         //       i stands for the node under consideration (node1)
         //       j (node2), k (node3), l (node4) are the other nodes of the quad.
         //   However, in the code, we use the following notation:
         //       id1 is the node under consideration (node1), id2, id3, id4 are the local indices of the other nodes (node2, node3, node4)
         const int id2 = elem->nxq(1 + id1);
         const int id3 = elem->nxq(2 + id1);
         const int id4 = elem->nxq(3 + id1);
         Node *node2 = elem->nodes[id2];
         Node *node3 = elem->nodes[id3];
         Node *node4 = elem->nodes[id4];   
         
         elem->getLocalHorizontalNodePositionsQuad(x0, x1, x2, x3);

         // integrals of phi_i*phi_j, phi_i*phi_k, phi_i*phi_l
         const double ipp_ii = elem->intPhiPhiQuad(id1, id1, x0, x1, x2, x3);
         const double ipp_ij = elem->intPhiPhiQuad(id1, id2, x0, x1, x2, x3);
         const double ipp_ik = elem->intPhiPhiQuad(id1, id3, x0, x1, x2, x3);
         const double ipp_il = elem->intPhiPhiQuad(id1, id4, x0, x1, x2, x3);

         // integrals of (grad phi_i) * (grad phi_i), (grad phi_i) * (grad phi_j), (grad phi_i) * (grad phi_k), (grad phi_i) * (grad phi_l)
         const double igpgp_ii = elem->intGradPhiGradPhiQuad(id1, id1, x0, x1, x2, x3);
         const double igpgp_ij = elem->intGradPhiGradPhiQuad(id1, id2, x0, x1, x2, x3);
         const double igpgp_ik = elem->intGradPhiGradPhiQuad(id1, id3, x0, x1, x2, x3);
         const double igpgp_il = elem->intGradPhiGradPhiQuad(id1, id4, x0, x1, x2, x3);


         // Computes the different components of integral I_1 of each node with respect to node i.
         const double I1_ii = ipp_ii + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ii;
         const double I1_ij = ipp_ij + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ij;
         const double I1_ik = ipp_ik + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ik;
         const double I1_il = ipp_il + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_il;

         node1->addToSparseMatrix(node1, I1_ii);
         node1->addToSparseMatrix(node2, I1_ij); 
         node1->addToSparseMatrix(node3, I1_ik);
         node1->addToSparseMatrix(node4, I1_il);
  
         #pragma unroll
         for (int l=0; l<4; l++) {

            const Node *node = elem->nodes[l];

            // Additional variation in eta due to bed deformation or vertical velocity
            double extraEta = 0.0;
            if      constexpr(imm == 2) extraEta = dt*ubedNormalComp;
            else if constexpr(imm == 1) extraEta = node->bdef2 - node->bdef1; 

            const double ipp = elem->intPhiPhiQuad(id1, l, x0, x1, x2, x3);
            const double igpgp = elem->intGradPhiGradPhiQuad(id1, l, x0, x1, x2, x3);
            node1->qel += (node->eta2 + extraEta)*ipp - node->eta2*fricCoeff*igpgp;

         }         
      } 


     node1->qel += (dt/elem->i34)*elem->vsource;
      

   } // for all elements connected to the node
   

   // Boundary integrals I_{2,3,5,6}; they all vanish at land boundaries
   // I_2 and I_6 are not needed if essential b.c. are enforced by elminating rows and columns
   NodeBoundary *nodeBnd = node1->nodeBnd;
   if (nodeBnd != nullptr && !nodeBnd->isLand) {
   // If it is an open boundary node...

   // XXXXX DO THIS MORE NATURALLY IN A SideBnd kernel?

      // At least one open bnd side, maybe two
      for (int l=0; l<2; l++) { 

         // Gets each of the sides sharing the node. 
         SideBoundary *sideBnd = nodeBnd->sideBnd(l);

    

         // XXXX CAN WE DO THIS MORE NATURALLY IN A SideBnd Loop?
         if (sideBnd != nullptr) {
            // If the side is on a boundary...

            const Node *node2 = sideBnd->side->otherNode(node1); // Other node of the side

            // Get the boundary segment
            const BoundarySegment *sideSegment = sideBnd->segment;
            const Side *side = sideBnd->side;

            // Integral I_3
            if (sideSegment->ifltype != 0) {
            // Natural or Flather b.c.

               // Calculate I_3 even if i is on essential b.c. so as to check symmetry later
               // especially for Flather b.c.
               double ri3 = 0.0;
               if (!side->isDry) {
                  // If the side is not dry
                  
                  // Computes the side flux of (uth, vth)
                  double Vn = 0.0;
                  for (int k = side->kb; k < nvrt - 1; k++) {
                     const double vn1 = dot(sideBnd->ubc[k  ], side->n); // outer normal
                     const double vn2 = dot(sideBnd->ubc[k+1], side->n);

                     Vn += 0.5*(vn1 + vn2)*(side->zc[k+1] - side->zc[k]);
                  }

                  ri3 = 0.5*side->len*Vn;
               }

               if (sideSegment->ifltype == -1) { // Flather 1

                  assert(node1->eta2 >= -98.0 && node2->eta2 >= -98.0 && "Mismatch 1") ;
                  assert(side->depth() > 0.0 && "Negative depth at Flather boundary");

                  const double con0 = (1.0/6.0)*side->len*sqrt(grav*side->depth()); // For coefficient matrix
                  ri3 -= con0*(2.0*node1->eta2 + node2->eta2);
                  
                  node1->addToSparseMatrix(node1, (thetai*dt*con0)*2.0);
                  node1->addToSparseMatrix(node2, (thetai*dt*con0));

               } 
               else if (sideSegment->ifltype == -2) { // Discharge relation (outgoing only)
                  // Reset ri3
                  ri3 = 0.0;

                  assert(polyDischarge != nullptr && "polyDischarge is not initialized. This part might not even be implemented.");
                  assert(sideSegment->clen != 0.0 && "clen is 0.0, this means that the segment is completely wetted");

                  const double etam = 0.5*(node1->eta2 + node2->eta2);

                  const double tmp = polyDischarge->eval(etam)/sideSegment->clen; // XXXX IMPLEMENT Poly. It may need the grade (4 in this case)
                  assert(tmp >= 0.0 && "bar{f} < 0 at discharge boundary");

                  const double con0 = (1.0/6.0)*side->len*tmp*thetai*dt;
                  node1->addToSparseMatrix(node1, con0*2.0);
                  node1->addToSparseMatrix(node2, con0);

               } // if Flather or Discharge

               node1->qel -= thetai*dt*ri3;  


            } // I_3
          

            // Integral I_5
            if(!side->isDry) {

               const double Unbar = dot(side->U, side->n);
               node1->qel -= 0.5*(1.0 - thetai)*dt*side->len*Unbar; 

            } // I_5

         } // if is boundary side

      } // for l

   } // if isbnd


}



extern "C" {

   __global__ void gatherBFromMesh(double* b) {
      // b[local] <- qel  (local = idx_node on this device)
      NODE(node);
      b[idx_node] = node->qel;
   }

   
   __global__ void scatterXToMesh(double* x) {
      // mX <- x[local]
      NODE(node);
      node->mX = x[idx_node];
   }

   __global__ void setupCGSolver(double* data_d, double* b_d, 
                                 int* indices, int* indptr) {
      // XXXXXXXXXXXX REMOVE!!!
      assert(false && "setupCGSolver has beed deprecated");

      // NODE(node);
      
      // const int idx = node->idx;
      
      // // Computes the index in the nzz array      
      // const int nzzIdx = (node->nodes - mesh->nodeNeighNodes) + idx; // idx is there to account for the diagonal element that is not stored in the nodeNeighNodes array


      // indptr[idx] = nzzIdx;
      // indices[nzzIdx] = idx;

      // data_d[nzzIdx]    = node->coeffDiag;

      // b_d[idx] = node->qel;

      // for (int k = 0; k < node->numNodes; ++k) {
      //    const Node* nb = node->nodes[k];

      //    const int j = nb->idx;

      //    data_d[nzzIdx + k + 1]    = node->coeffMatrix[k];

      //    indices[nzzIdx + k + 1] = j;
      // }

   }

   __global__ void setupMultiGPUCGSolver(double* data_d, double* b_d, 
                                         int *indices, int *indptr) {
   // Builds all the A_ij matrices and b_i and x_i vectors for the solver.

      NODE(node);

      // Global and local indices of the node.
      // In this function, I stands for the global index of the current node,
      // and J stands for the global index of the neighbor node under consideration.
      // i and j are the local indices of the current node and the neighbor node under consideration.
      // Notice that i and j may be local to different GPUs.
      // idxPart is the number of the partition of the *neighbor* node. There is
      // a biunivocal correspondence between the partitions and the GPUs.
      const int I = node->idx;      
      const int i = idx_node;       
      const int partI = node->part; // Partition of the current node.
      

      setBValue(partI, i, node->qel);

   
 
   }
}



   __global__ void solverCopyResult(double* x) {
      // Copies a solution vector to the member variable mX. XXXX It might be possible to remove all the unnecessary references to mX, mR, mP, mAp, etc.
      NODE(node);
      const int idx = node->idx;

      node->mX = x[idx];
   }


   __global__ void rollUpElevation(void) {

      NODE(node);

      node->eta1 = node->eta2;

   }


   __global__ void backtrackAllNodes(void) {

      // max/min # of sub-divisions
      const int ndelt_max = (int)max(1.0, dt/dtb_min);
      const int ndelt_min = (int)max(1.0, dt/dtb_max); 

      SIDE(side);

      const int kb = side->kb;


      // Initialize btrack velocities to the actual velocities at the side 
      for (int k = 0; k < nvrt; k++) side->btu[k] = side->u2[k];
  
      if (side->isDry) return;
     
// addDebugPlot(side->center, side->center, 0xFF001F, 'x');    

      // Checks that at least one element is not dry
      const Element *elem1 = side->elems[0];
      const Element *elem2 = side->elems[1];

// if (!((elem1 != nullptr && !(elem1->isDry)) || (elem2 != nullptr && !(elem2->isDry)))) return;      

      // assert ((elem1 != nullptr && !(elem1->isDry)) || (elem2 != nullptr && !(elem2->isDry)) && "backtrackAllNodes: non dry side without two dry elements.");
      
      // XXXXXX
      // !swild_tmp to store global coord. of the starting side
      // swild_tmp(1)=xcj(isd0); swild_tmp(2)=ycj(isd0); swild_tmp(3)=zcj(isd0)
      // !swild10_tmp to store frame at starting pt for ics=2 (not used for ics=1)
      // !Use sframe2 
      // swild10_tmp(1:3,1:3)=sframe2(:,:,i) !pframe(:,:,isidenode(1,isd0)) 


      for (int k = kb; k < nvrt; k++) {
      //  For ics=2, the coordinates are all in the side's local frame

         const int iadvf = side->iadv();  // Advection flag   
         assert(nadv == 1 && "Advection flag is only implemented for nadv = 1");

         const double vmag2 = norm2(side->u2[k]);

         if (vmag2 <= velmin_btrack_sq) {  
            // No activity 

            side->btu[k] = side->u2[k];  // I don't think this is necessary because it has already been set to this value XXXX
            
         } else {  // do btrack

            // Compute number of sub-divisions based on local gradients of the velocity field

            // XXX IS THIS NECESSARY?
            // The gradient at a side is the average of the gradients at the elements sharing the side.
            double suma = 0.0;
            int icount = 0;
            #pragma unroll
            for (int l=0; l<2; l++) {
            // For each element sharing the side...

               const Element *elem = side->elems[l];
                 
               if (elem == nullptr) continue;

               icount++;

               // Original Fortran code comment: "Not strictly along z; in local frame for ics=2"
               // I understand it to mean that u2[k] is assumed to be horizontal-only when transformed to
               // the local frame. This is not strictly true, but it is a reasonable approximation.

               double utmp[4];
               double vtmp[4];
               #pragma unroll
               for (int j = 0; j < elem->i34; j++) {

                  const Node *node = elem->nodes[j];

                  if constexpr(ics==1) {
                     utmp[j] = node->u2[k].x;
                     vtmp[j] = node->u2[k].y;  

                  } else {
                     // 'Project' the velocity from the node frame to the element frame
                     const Vector2 tmp_u2 = node->frame.projHorVec(node->u2[k], elem->frame);
                     utmp[j] = tmp_u2.x;
                     vtmp[j] = tmp_u2.y;

                  }
               } // for j

               const Vector2 grad_u = elem->grad2(utmp); 
               const Vector2 grad_v = elem->grad2(vtmp);
   
               suma += dt*sqrt(dot(grad_u, grad_u) + dot(grad_v, grad_v));
            } // for l

            assert(icount > 0 && "Impossible 77");

            // '4' is somewhat arbitrary
            const int ndelt = clamp(int(suma/icount)*4, ndelt_min, ndelt_max); 

            const double dtbk = dt/ndelt; // Target btrack step; may be smaller sometimes
            Vector3 u_int = {0.0, 0.0, 0.0};

            side->btrack(k, dt, dtbk, u_int);


            side->btu[k] = {u_int.x, u_int.y};
               
         } // if vmag (do backtrack)


      } // for k=kb,nvrt-1

   }

} // extern "C"

template<int nc>
__device__ void tridag_sch(int n, double *a, double *b, double *c, double *rhs, double *u, double *gam) {

// Solves a tridiagonal system. It was adapted from "Numerical Recipes in FORTRAN (pp.43 ).
//
// a,b,c,rhs: input vectors and are not modified by this program.
// b is the main diagonal, a below and c above. a(1) and c(n) are not used.
// rhs is the r.h.s.
// n: actual rank of the system.
// rows in u() and rhs();
// u: output with nc columns (transposed to rows)
// gam: a working array.
//
// NOTES:
// This is a template function that allows for different numbers of columns in the rhs vector.
// [n, nc] is the shape of Fortran-type arrays rhs and u.

   assert(n >= 1       && "TRIDAG: n must be >= 1");
   assert(b[0] != 0.0  && "TRIDAG: b[0]=0");

   double bet = b[0];
   if constexpr (nc == 1) {
      u[0] = rhs[0]/bet;
   } else {
      for (int l = 0; l < nc; l++) u[l] = rhs[l]/bet;
   }

   for (int j = 1; j < n; j++) {
      gam[j] = c[j - 1]/bet;
      bet = b[j] - a[j]*gam[j];

      assert (bet != 0.0 && "TRIDIAGONAL SOLVER: bet=0");

      if constexpr (nc == 1) {
         u[j] = (rhs[j] - a[j]*u[j - 1])/bet;
      } else {
         for (int l = 0; l < nc; l++) {
            u[j*nc + l] = (rhs[j*nc + l] - a[j]*u[(j - 1)*nc + l])/bet;
         }
      }
   }

   // Back-substitution
   for (int j = n - 2; j >= 0; j--) {
      if constexpr (nc == 1) {
         u[j] = u[j] - gam[j + 1]*u[j + 1];
      } else {
         for (int l = 0; l < nc; l++) {
            u[j*nc + l] = u[j*nc + l] - gam[j + 1]*u[(j + 1)*nc + l];
         }
      }
   }
}



// Helper function that contains the core algorithm
template<typename XAccessor, typename YAccessor>
__device__ void initCubicSplineImpl(int k1, int k2, XAccessor getX, YAccessor getY, 
                                    double  yp_top, double  yp_bot, 
                                    CubicSpline &spline) {
   // Compute the 1st and 2nd derivatives of y = f(x) using continuous cubic splines.
   // 
   // This function uses the General expression for a C2 interpolating cubic spline (search for "C2 interpolating cubic spline" in Wikipedia)
   // This method uses the 2nd derivatives and the y values at the endpoints to generate the polynomial coefficients for each level.
   //
   // NOTE: (x,y) here are just any function y=f(x). In fact we use these splines for vertical interpolation.
   // 
   // Inputs: 
   //            n: dimension for x etc; # of data pts (>=2);
   //            getX, getY: functions to access x,y coordinates
   //            yp_top and yp_bot: 1st derivatives at x[0] and x[n-1];
   //     Output: 
   //            spline (CubicSpline): spline structure with x,y,y',y''

   const int n = k2 - k1 + 1;
   assert(n >= 2 && "CUBIC SPLINE: n < 2");

   // Allocate memory for the splines
   assert(spline.x == nullptr && spline.y == nullptr && spline.yp == nullptr && spline.ypp == nullptr && "CUBIC SPLINE: x,y,yp,ypp already allocated");
   spline.x   = new double[n];
   spline.y   = new double[n];
   spline.yp  = new double[n];
   spline.ypp = new double[n];

   for (int k=k1; k<=k2; k++) {
      spline.x[k] = getX(k);
      spline.y[k] = getY(k);
   }

   // This is to simplify expressions later
   double *x = spline.x;
   double *y = spline.y;

   double *a = nullptr;
   double *b = nullptr;
   double *c = nullptr;
   double *gam = nullptr;
   double *rhs = nullptr;
   

      a = new double[n];   // Lower diagonal
      b = new double[n];   // Diagonal
      c = new double[n];   // Upper diagonal
      gam = new double[n]; // Gamma, used by the tridiagonal solver
      rhs = new double[n]; // Right-hand side

   for (int k=k1; k<=k2; k++) {
      if (k == k1) {
         // Bottom level
         assert(x[k+1] - x[k] != 0.0 && "CUBIC SPLINE: bottom problem");
         b[k] = (x[k+1] - x[k])/3.0;
         c[k] = b[k]/2.0;
         rhs[k] = (y[k+1] - y[k])/(x[k+1] - x[k]) - yp_bot;
      } else if (k == n - 1) {
         // Top level (surface)
         assert(x[k] - x[k-1] != 0.0 && "CUBIC SPLINE: surface problem");
         b[k] = (x[k] - x[k-1])/3.0;
         a[k] = b[k]/2.0;
         rhs[k] = -(y[k] - y[k-1])/(x[k] - x[k-1]) + yp_top;
      } else {
         // Intermediate levels
         assert(x[k] - x[k-1] != 0.0 && x[k+1] - x[k] != 0.0 && "CUBIC SPLINE: middle problem");
         b[k] = (x[k+1] - x[k-1])/3.0;
         a[k] = (x[k  ] - x[k-1])/6.0;
         c[k] = (x[k+1] - x[k  ])/6.0;
         rhs[k] = (y[k+1] - y[k  ])/(x[k+1] - x[k  ]) - 
                  (y[k  ] - y[k-1])/(x[k  ] - x[k-1]);
      }
   }

   // Solve the tridiagonal system
   tridag_sch<1>(n, a, b, c, rhs, spline.ypp, gam);
   
   spline.yp[k1] = yp_bot;
   spline.yp[k2] = yp_top;

   // Compute rest of 1st derivatives using the splines
   for (int k=k1+1; k<k2; k++) {
      spline.yp[k] = (y[k+1] - y[k])/(x[k+1] - x[k]) - 
                     (x[k+1] - x[k])/6.0*(2.0*spline.ypp[k] + spline.ypp[k+1]);
   }
   
   delete[] a;
   delete[] b;
   delete[] c;
   delete[] gam;
   delete[] rhs;

}

// Array version
__device__ void initCubicSpline(int k1, int k2, double *x, double *y, 
                                double y_top, double y_bot, 
                                CubicSpline &spline) {
   auto getX = [x] (int k) { return x[k]; };
   auto getY = [y] (int k) { return y[k]; };
   initCubicSplineImpl(k1, k2, getX, getY, y_top, y_bot, spline);
}

// Lambda version - allows custom access patterns
template<typename XFunc, typename YFunc>
__device__ void initCubicSpline(int k1, int k2, XFunc getX, YFunc getY,
                                double y_top, double y_bot,
                                CubicSpline &spline) {
   initCubicSplineImpl(k1, k2, getX, getY, y_top, y_bot, spline);
}


__device__ double evalCubicSpline(CubicSpline spline, double z, bool extrapolate = false)
{
   assert("Not implemented");

   return 0.0;
}

// __device__ void vertLinearInterp(int nc, double zt, int k1, int k2, double *za, double *s_in, double *s_out, int *ibelow) {
// //     Routine to do vertical linear interpolation in z
// //     Inputs:
// //
// //       nc: actual # of variables 
// //       k1,k2: lower and upper limits for za, s_in 
// //       zt: desired interpolation level
// //       za[0:nvrt-1]: z-coordinates for s_in (must be in ascending order)
// //       s_in[0:nc-1,0:nvrt-1]: values to be interpolated from; dimensions must match driving program

// //     Outputs:
// //       s_out[0:nc-1]: interpolated value @ z=zt (bottom value if ibelow=1). Constant extrapolation
// //                   is used below bottom or above surface.
// //       ibelow: flag indicating if zt is below za(k1)
// //

//    assert(k1 >= 0 && k2 < nvrt && "vertLinearInterp: k1 or k2 out of bounds");
//    assert(k1 <= k2 && "vertLinearInterp: k1 > k2");
   
//    ibelow = 0;

//    if (zt < za[k1]) {
//       // Below bottom
//       ibelow = 1;
//       for (int i = 0; i < nc; i++) s_out[i] = s_in[i + k1*nc];

//    } else if(zt>=za[k2]) {
//       // Above surface
//       for (int i = 0; i < nc; i++) s_out[i] = s_in[i + k2*nc];

//    } else {
//       // Normal case (between bottom and surface)
      
//       // Finds the level index using binary search
//       int l1 = k1;
//       int l2 = k2;
//       int k = (l1 + l2)/2;
//       while (l2 - l1 > 1) {  // Stop when we have adjacent levels
//           int k = (l1 + l2)/2;
//           if (zt >= za[k]) l1 = k;
//           else             l2 = k;
//       }
//       // Now l1 and l2 are adjacent levels and zt is between them

//       assert(za[l2] - za[l1] > 0.0 && "vertLinearInterp: za[l2] - za[l1] <= 0");

//       const double zrat = (zt - za[l1])/(za[l2] - za[l1]);
//       for (int i = 0; i < nc; i++) s_out[i] = s_in[i + l1*nc]*(1.0 - zrat) + s_in[i + l2*nc]*zrat;
//    }
// }


// __device__ void vertLinearInterp(double zt, int k1, int k2, double *za, double *s_in, double *s_out, int *ibelow) {
// //     Routine to do vertical linear interpolation in z
// //     Inputs:
// //
// //       k1,k2: lower and upper limits for za, s_in 
// //       zt: desired interpolation level
// //       za[0:nvrt-1]: z-coordinates for s_in (must be in ascending order)
// //       s_in[0:nvrt-1]: values to be interpolated from; dimensions must match driving program
// //     Outputs:
// //       s_out: interpolated value @ z=zt (bottom value if ibelow=1). Constant extrapolation
// //                   is used below bottom or above surface.
// //       ibelow: flag indicating if zt is below za(k1)
// //

//    assert(k1 >= 0 && k2 < nvrt && "vertLinearInterp: k1 or k2 out of bounds");
//    assert(k1 <= k2 && "vertLinearInterp: k1 > k2");
   
//    ibelow = 0;

//    if (zt < za[k1]) {
//       // Below bottom
//       ibelow = 1;
//       s_out = s_in[k1];

//    } else if(zt>=za[k2]) {
//       // Above surface
//       s_out = s_in[k2];

//    } else {
//       // Normal case (between bottom and surface)
      
//       // Finds the level index using binary search
//       int l1 = k1;
//       int l2 = k2;
//       int k = (l1 + l2)/2;
//       while (l2 - l1 > 1) {  // Stop when we have adjacent levels
//           int k = (l1 + l2)/2;
//           if (zt >= za[k]) l1 = k;
//           else             l2 = k;
//       }
//       // Now l1 and l2 are adjacent levels and zt is between them

//       assert(za[l2] - za[l1] > 0.0 && "vertLinearInterp: za[l2] - za[l1] <= 0");

//       const double zrat = (zt - za[l1])/(za[l2] - za[l1]);
//       s_out = s_in[l1]*(1.0 - zrat) + s_in[l2]*zrat;
//    }
// }

__device__ __forceinline__ int localIdx(Node* node) {
   const int part = node->part;
   const int localIdx = node - mesh->nodes[part]; // The local index of the node in the partition is computed as the difference between the node and the start of the node memory block in the partition
   return localIdx;
}

extern "C" {

   __global__ void scaleGravity(double loadtide_coef) {


      NODE(node);

      // SAL option: scale gravity
      if constexpr(iloadtide==2) {  // simple const
         node->grav2 = grav*(1.0 - loadtide_coef);

      } else if constexpr(iloadtide==3) {  // Stepanov & Hughes (2004)
        
          const double depth = max(node->depth, 0.0);
          const double tmp1 = sqrt(depth);
          double beta_bar = -9.8169e-3 + 1.8289e-3*tmp1 + 4.3787e-4*depth - 2.9042e-5*depth*tmp1 + 
                             6.6038e-7*depth*depth - 4.7393e-9*depth*depth*tmp1 - 1.9354e-11*depth*depth*depth + 
                             2.6969e-13*depth*depth*depth*tmp1;
          
          beta_bar=max(0.0, min(loadtide_coef, beta_bar));
          


          node->grav2 = grav*(1.0 - beta_bar);

      } else {  // iloadtide=0,1
          node->grav2 = grav;

      }

   }


   __global__ void setWinds(double *xWind, double *yWind, double *pressAir) {
      NODE(node);

      const int idx = node->idx;
      Vector2 wind = {xWind[idx], yWind[idx]};

      if constexpr(nws != 0 && iwindoff != 0) {
         assert(false && "setWinds: windfactor not implemented (nws != 0 && iwindoff != 0)");
      }

      // Sets the wind speed (rotated to the mesh angle) and surface stress vector due to wind.
      node->setWind(wind);

      node->pr = pressAir[idx];
      node->p_air = pressAir[idx]; // XXXXXX is this the same as pr?

   }

   __global__ void computeTidalPotential(double time, double* tfreq, double* tamp, double* tnf, double* tear, int* jspc) {

      NODE(node);

      node->ETP = 0.0;

      if constexpr(ntip == 0) return;

      // For all tidal potential components
      for (int j=0; j<ntip; j++) {

         const double nodeETP = node->computeETP(j, time, tfreq, tamp, tnf, tear, jspc);

         node->ETP += nodeETP;

         if constexpr(iloadtide==1) { // loading tide
            assert(false && "Loading tide not implemented");
            // node->ETP += node->rloadtide[0,j]*cos(tfreq[j]*time - node->rloadtide[1,j]);
         }
         
      }

     
   }

   __global__ void getBoundaryConditions(double time) {
      // Get boundary conditions
      // XXX IMPLEMENT THIS
      // assert(true && "getBoundaryConditions not implemented");


      NODE_BND(nodeBnd);

      assert(nfltype <= 0 && nettype <= 0 && "getBoundaryConditions: Not implemented");
      return;

      const Node *node = nodeBnd->node;

    

   //    if(nettype>0) then
   //      if(time>th_time(1,2,1)) then !not '>=' to avoid last step
   //        ath(:,1,1,1)=ath(:,1,2,1)
   //        read(50,*) tmp,ath(1:nettype,1,2,1)
   //        th_time(1,1,1)=th_time(1,2,1)
   //        th_time(1,2,1)=th_time(1,2,1)+th_dt(1,1)
   //      endif !time

      
   //      rat=(time-th_time(1,1,1))/th_dt(1,1)
   //      if(rat<-small1.or.rat>1.d0+small1) then
   //        write(errmsg,*) 'STEP: rat out in elev.th:',rat,time,th_time(1,1:2,1),th_dt(1,1)
   //        call parallel_abort(errmsg)
   //      endif

   //      icount=0
   //      do k=1,nope_global
   //        if(iettype(k)==1) then
   //          icount=icount+1
   //          if(icount>nettype) call parallel_abort('Wrong counting 1')
   //          eth(1,k)=(1-rat)*ath(icount,1,1,1)+rat*ath(icount,1,2,1)
   //        endif
   //      enddo 
   //    endif !nettype

   //    if(nfltype>0) then
   //      if(time>th_time(1,2,2)) then
   //        ath(:,1,1,2)=ath(:,1,2,2)
   //        read(51,*) tmp,ath(1:nfltype,1,2,2)
   //        th_time(1,1,2)=th_time(1,2,2)
   //        th_time(1,2,2)=th_time(1,2,2)+th_dt(1,2)
   //      endif !time


   //      rat=(time - th_time(1,1,2))/th_dt(1,2)
   //      if(rat<-small1.or.rat>1.d0+small1) then
   //        write(errmsg,*) 'STEP: ratio out of range while interpolating &
   //   &flux.th. Probably times are not equally spaced or dt has changesd &
   //   &from a prior run (ratio, time, th times):',rat,time,th_time(1,1:2,2)
   //        call parallel_abort(errmsg)
   //      endif
   //      icount=0
   //      do k=1,nope_global
   //        if(ifltype(k)==1) then
   //          icount=icount+1
   //          if(icount>nfltype) call parallel_abort('STEP: wrong counting 2')
   //          qthcon(k)=(1.d0-rat)*ath(icount,1,1,2)+rat*ath(icount,1,2,2)
   //        endif
   //      enddo !k
   //    endif !nfltype

 

   //    if(nettype2>0) then
   //      if(time>th_time2(2,1)) then        
   //        ath2(:,:,:,1,1)=ath2(:,:,:,2,1)
   //        icount3=time/th_dt2(1)+2
   //        j=nf90_inq_varid(ncid_elev2D, "time_series",mm)
   //        if(j/=NF90_NOERR) call parallel_abort('step: time_series in elev2D.th.nc')
   //        j=nf90_get_var(ncid_elev2D,mm,ath2(1,1,1:nnode_et,2,1), &
   //  &(/1,1,1,icount3/),(/1,1,nnode_et,1/))
   //        if(j/=NF90_NOERR) call parallel_abort('step: time_series in elev2D.th.nc (2)')

   //        th_time2(1,1)=th_time2(2,1)
   //        th_time2(2,1)=th_time2(2,1)+th_dt2(1)
   //      endif !time


   //      rat=(time-th_time2(1,1))/th_dt2(1)
   //      if(rat< -small1.or.rat>1.d0+small1) then
   //        write(errmsg,*) 'STEP: rat out in elev2D.th:',rat,time,th_time2(1:2,1),th_dt2(1)
   //        call parallel_abort(errmsg)
   //      endif
   //      icount=0
   //      icount2=0
   //      do k=1,nope_global
   //        if(iettype(k)>=4) then
   //          icount=icount+1
   //          if(icount>nettype2) call parallel_abort('STEP: wrong counting 7')
   //          do j=1,nond_global(k)

   //            icount2=icount2+1
   //            if(icount2>nnode_et) call parallel_abort('STEP: wrong counting nodes')

   //            eth(j,k)=(1.d0-rat)*ath2(1,1,icount2,1,1)+rat*ath2(1,1,icount2,2,1)
   //          enddo !j
   //        endif
   //      enddo !k
   //    endif !nettype2

   //    if(nfltype2>0) then
   //      if(time>th_time2(2,2)) then
   //        ath2(:,:,:,1,2)=ath2(:,:,:,2,2)


   //        icount3=time/th_dt2(2)+2
   //        j=nf90_inq_varid(ncid_uv3D, "time_series",mm)
   //        if(j/=NF90_NOERR) call parallel_abort('step: time_series in uv3D.th.nc')
   //        j=nf90_get_var(ncid_uv3D,mm,ath2(1:2,1:nvrt,1:nnode_fl,2,2), &
   //  &(/1,1,1,icount3/),(/2,nvrt,nnode_fl,1/))

   //        if(j/=NF90_NOERR) call parallel_abort('step: time_series in uv3D.th.nc')



   //        th_time2(1,2)=th_time2(2,2)
   //        th_time2(2,2)=th_time2(2,2)+th_dt2(2)
   //      endif !time


   //      rat=(time-th_time2(1,2))/th_dt2(2)
   //      if(rat<-small1.or.rat>1.d0+small1) then
   //        write(errmsg,*) 'STEP: rat out in uv3D.th:',rat,time,th_time2(1:2,2),th_dt2(2)
   //        call parallel_abort(errmsg)
   //      endif
   //      icount=0
   //      icount2=0
   //      do k=1,nope_global
   //        if(iabs(ifltype(k))>=4) then
   //          icount=icount+1
   //          if(icount>nfltype2) call parallel_abort('STEP: wrong counting 6')
   //          do j=1,nond_global(k)
   //            icount2=icount2+1
   //            if(icount2>nnode_fl) call parallel_abort('STEP: wrong counting vel')

   //            uthnd(1:nvrt,j,k)=(1.d0-rat)*ath2(1,1:nvrt,icount2,1,2)+rat*ath2(1,1:nvrt,icount2,2,2) !ll frame if ics=2
   //            vthnd(1:nvrt,j,k)=(1.d0-rat)*ath2(2,1:nvrt,icount2,1,2)+rat*ath2(2,1:nvrt,icount2,2,2)
   //          enddo !j
   //        endif
   //      enddo !k
   //    endif !nfltype2
   }

   __global__ void resetBoundaryAreas(void) {
      // Reset boundary areas and lengths for each segment.

      // This should be done per segment, not per boundary side (each segment has multiple sides),
      // However we don't have a way to access the segments by themselves and this is simple and fast enough.
      SIDE_BND(sideBnd);

      BoundarySegment *segment = sideBnd->segment;

      segment->carea = 0.0;
      segment->clen = 0.0;
 
   }

   __global__ void computeBoundaryAreas(void) {
   // Calculation of cross-section areas and lengths for flow b.c.

      SIDE_BND(sideBnd);

      const Side *side = sideBnd->side;
      BoundarySegment *segment = sideBnd->segment;

      if (segment->ifltype != 0) {

         const Node *n1 = side->n1;
         const Node *n2 = side->n2;

         const double htot = side->htot2();

         if (!side->isDry) {
            atomicAdd(&segment->carea, htot*side->len);
            atomicAdd(&segment->clen, side->len);
         }
      }

   }


   __global__ void computeFlowBcVelocities(double time) {

      SIDE_BND(sideBnd);

      const Side *side = sideBnd->side;
 
      const bool isDry = side->isDry;
      const bool isLand = sideBnd->segment->isLand;

   
      for (int k = 0; k < nvrt; k++) sideBnd->ubc[k] = {0.0, 0.0};

      if (isDry || isLand) return;

      // For wet sides in open boundaries...
      const Vector2 n = side->n;
      const BoundarySegment *segment = sideBnd->segment;

      assert(segment != nullptr && "computeFlowBcVelocities: segment is nullptr");

      const int ifltype = segment->ifltype;

      const double r_uth0 = ramp*segment->qthcon/segment->carea;  // r_ stands for ramp-up
      const Vector2 r_nuth0 = r_uth0*n;


      if (ifltype == 1 || ifltype == 2) {
      
         assert(segment->carea > 0.0 && "computeFlowBcVelocities: dry bnd side");

         const double unth0 = segment->qthcon*ramp/segment->carea;
         
         for (int k = 0; k < nvrt; k++) sideBnd->ubc[k] += unth0*n;  // This case has no ramp-up
         
      } else if (ifltype == -1) { // Flather 1

         for (int k = 0; k < nvrt; k++) {

            // uth is the normal vel.; no ramp up
            const double uth = sideBnd->nuth(k);
            assert(uth >= -98.0 && uth <=  98.0 && "computeFlowBcVelocities: wrong velocity b.c. at side");

            sideBnd->ubc[k] += uth*n;
         }

      } else if (ifltype == 3 || ifltype == 5) {
         // vnth0=0 !normal vel.

         Vector2 r_u1_bar = {0.0, 0.0};
         // Computes the contribution to the velocity from all tide frequencies
         for (int jfr = 0; jfr < nbfr; jfr++) { // For all tide frequencies

            const Vector2 ufa  = sideBnd->ufa(jfr); 
            const Vector2 umo  = sideBnd->umo(jfr); 

            const int ncyc = int(amig[jfr]*time/(2.0*pi));
            const double argx = amig[jfr]*time - (2.0*pi)*ncyc + face[jfr] - ufa.x;
            const double argy = amig[jfr]*time - (2.0*pi)*ncyc + face[jfr] - ufa.y;
            r_u1_bar += ramp*ff[jfr]*Vector2{umo.x*cos(argx), umo.y*cos(argy)};

         }


         for (int k = 0; k < nvrt; k++) sideBnd->ubc[k] += r_u1_bar;

      } // if (ifltype = 3 || ifltype == 5



   }

   __global__ void computeElevationBoundaryConditions(double time) {
      // Compute elevation essential boundary conditions 

      NODE_BND(nodeBnd); // Each boundary node.

      assert(nodeBnd->segment != nullptr && "computeElevationBoundaryConditions: nodeBnd->segment is nullptr");

      const Node *node = nodeBnd->node;
      const BoundarySegment *segment = nodeBnd->segment;
      const int iettype = segment->iettype;
      

      // Elevation initial condition at boundary nodes. It is ramped down to zero
      double eic = 0.0;
      if (nramp_elev == 1) eic = nodeBnd->etaic;  // XXXX is etaic ever set?  
      else                 eic = 0.0;

      const double r_eic = (1.0 - ramp)*eic;   // r_ stands for ramp
      const double r_eth = ramp*nodeBnd->eth;
      double r_eta1_bar = 0.0;

      // Compute tide elevation
      if (iettype == 3 || iettype == 5) {   

         if constexpr(iloadtide==1) { // loading tide
            assert(false && "Loading tide not implemented");
            // node->ETP += node->rloadtide[0,j]*cos(tfreq[j]*time - node->rloadtide[1,j]);
         }
         
         for (int jfr = 0; jfr < nbfr; jfr++) { // For all tide frequencies

            const double efa  = nodeBnd->efa[jfr]; 
            const double emo  = nodeBnd->emo[jfr]; 

            const int ncyc = int(amig[jfr]*time/(2.0*pi));
            const double arg = amig[jfr]*time - (2.0*pi)*ncyc + face[jfr] - efa;
            r_eta1_bar += ramp*ff[jfr]*emo*cos(arg);


         } // for jfr

      } // if iettype

      switch (iettype) {
         case 1:
         case 2:
         case 4:
            nodeBnd->elbc = r_eic + r_eth;
            break;
         case 3:
            nodeBnd->elbc = r_eic + r_eta1_bar;
            break;
         case 5:
            nodeBnd->elbc = r_eic + r_eta1_bar + r_eth;
            break;
      } // switch


      //Add inverse barometric effects
      if (iettype != 0 && inv_atm_bnd == 1) {
         nodeBnd->elbc += ramp*(prmsl_ref - node->pr)/grav/rho0;
      }

   }


   __global__ void computeHorizontalViscosity(void) { 

      assert(ihorcon == 0 || ics == 1 && "computeHorizontalViscosity: implemented only for ics=1");

      SIDE(side);

      if (side->isDry) return;

      for (int k = 0; k < nvrt; k++) {
         side->nuGradU[k] = {0.0, 0.0};
      }

      if constexpr(ihorcon == 1) {
         for (int k = side->kb; k < nvrt - 1; k++) {
            // Computes the Laplacian by using the simplified formula presented in the paper:
            // SELFE: A semi-implicit Eulerian-Lagrangian finite-element model for cross-scale ocean circulation
            // By Zhang and Baptista.

            // sum(nu*(u_i - u_c))/area, where u_i is the velocity at the i-th surrounding side and u_c is the velocity at the side under consideration.
            // This is a standard discretization of the Laplacian operator.
            // Following the paper, nu is estimated as: gamma*area/dt, where gamma is a non-dimensional parameter.
            // The area is the area of the cell we use for this calculation (it's not the area of an element). However, since it is in the numerator and the denominator, it cancels out.
            // Gamma is hvis_coef0 in the code.
            // The value of gamma needs to be < 0.5.

            // For the computation at the side center, we use the velocities at the two sides that are connected to the 
            // side under consideration in each of the two possible elements that share the side.
            Vector2 sumU = {0.0, 0.0};
            for (int l = 0; l < 2; l++) {
               const Element *elem = side->elems[l];

               if (elem == nullptr) continue;  // This happens at the domain boundaries

               // Find the local index of the side in the element
               const int idxSide = elem->sideLocalIdx(side);

               // Find the two sides that are connected to the side under consideration
               const Side *s1 = elem->sides[elem->nxq(idxSide + 1)];
               const Side *s2 = elem->sides[elem->nxq(idxSide - 1)];

               sumU += s1->u2[k] + s2->u2[k] - 2*side->u2[k];   // XXXX Is this correct?

            }

            side->nuGradU[k] = hvis_coef0*sumU/dt;

            assert(false && "computeHorizontalViscosity: Check if the above is correct");
   
         }
      }
      else if constexpr(ihorcon == 2) {
         assert(false && "computeHorizontalViscosity: implemented only for ihorcon == 1");
      }


   }


   __global__ void setEssentialBC(void) {  // XXXXXXX CONVERT TO NODE_BND????????

      // Set essential boundary conditions 
      NODE(node); // Each boundary node
      NodeBoundary *nodeBnd = node->nodeBnd;
      
 
      // Check if this node has essential BC
      if (nodeBnd && nodeBnd->isEssentialBC) {

         // Verify BC value is valid
         assert(nodeBnd->elbc >= -9998.0 && "setEssentialBC: wrong b.c. at node");

       
         // Set diagonal to 1 and zero off-diagonals
         for (int j = 0; j <= node->numNodes; j++) {
               node->coeffRow[j] = 0.0;
         }

         *(node->coeffDiag) = 1.0;


         // Set solution and RHS to BC value
         node->mX  = nodeBnd->elbc;
         node->qel = nodeBnd->elbc;

      } else {
         // Not an essential BC node. Still need to modify RHS and columns that belong to neighbors with essential BCs

         for (int j = 0; j < node->numNodes; j++) {
               const Node* nbNode = node->nodes[j];
   	         assert(nbNode != nullptr && "setEssentialBC: nbNode is nullptr");

               const NodeBoundary* nbNodeBnd = nbNode->nodeBnd;
               assert(nbNodeBnd != nullptr && "setEssentialBC: nbNodeBnd is nullptr");

               if (nbNodeBnd && nbNodeBnd->isEssentialBC) {
                  // Verify neighbor's BC value
                  assert(nbNodeBnd->elbc >= -9998.0 && "setEssentialBC: wrong b.c. at node");

                  // Subtract contribution to RHS
                  node->qel -= node->coeffRow[node->csrSlot[j]]*nbNodeBnd->elbc;
                  node->coeffRow[node->csrSlot[j]] = 0.0;  
               }
         }
// printf("setEssentialBC: node->qel = nodeBnd->elbc - node->coeffMatrix[j]*nbNodeBnd->elbc\n");
      }
   }


   __global__ void updateElevation(void) {
      NODE(node);


      node->eta2 = node->mX;

      node->updateZCoordData();
   }


   __global__ void computeVerticallyIntegratedVelocity(void) {

      SIDE(side);

      side->U = {0.0, 0.0}; 
         
      if(side->isDry) return;

      // Vertically integrated velocity U
      for (int k = side->kb; k < nvrt - 1; k++) { 
         const Vector2 uc = 0.5*(side->u2[k] + side->u2[k + 1]);
         side->U += side->Dz(k + 1)*uc;  
      }

   }

   __global__ void computeBottomDragVegetation(void) {
      // UNFINISHED. CHECK THIS:
      // How many kernels are needed? Why so many for loops?
      // Is the loop on sides needed?

      NODE(node);


      const int    kb      = node->kb;
      const double htot    = node->htot2(); 
      
      if (node->isDry) {
         // Default values for dry nodes
         node->veg_alpha_vert_mean = node->veg_alpha0;
         if constexpr(iveg != 0) for (int k = kb; k < nvrt; k++) node->veg_alpha3D[k] = node->veg_alpha0;

         return;
       
      }

      // Computes bottom drag Cd based on the selected nchi
      node->computeBottomDragCd();
      


      //   Bottom stress in m^2/s^2
      if (node->rho[kb + 1] >= -98.0)
      {
         const double mag_u = length(node->u2[kb + 1]);
         node->tau_bot = node->rho[kb + 1]*node->Cd*mag_u*node->u2[kb + 1]; // units: kg/m/s^2 (Pa)

         // XXX In the original code there is an additional component of tau_bot, used in icm_sfm.F90 with this comment:
         // !sediment erosion: this part needs more documentation
      }
    


      // Compute vegetation parameters
      //--------------------------------

      // SAV constants
      constexpr double veg_cfk   = 0.07; //Shimizu & Tsujimoto (1994)
      constexpr double veg_cfpsi = 0.16;

      // Add vertical variation to veg_alpha and compute vertical mean etc.
      if constexpr(iveg == 1) { //  Vertical variation specified by functions

         // Compute the vegetation "density" \alpha = D*N*Cd/2
         // D, N and Cd are given as piecewise linear functions: veg_vert_scale_cd, veg_vert_scale_N, veg_vert_scale_D
         int kbin = 0; // bin index
         for (int k = kb; k < nvrt; k++) {
            // Computes node height above the bottom level
            const double tmp_h = node->z[k] - node->z[kb]; // >=0
            assert(false && "Vertical variation of vegetation not implemented");

            // if (tmp_h >= veg_vert_z[nbins_veg_vert]) {
            //    cff1 = veg_vert_scale_cd[nbins_veg_vert];
            //    cff2 = veg_vert_scale_N [nbins_veg_vert];
            //    cff3 = veg_vert_scale_D [nbins_veg_vert];
            // } else {
            //    // Find the bin at which the current vertical level is
            //    while (kbin<nbins_veg_vert) {
            //       const double zn = node->z[kbin];
            //       if (tmp_h >= veg_vert_z[kbin] && tmp_h <= veg_vert_z[kbin+1]) {
            //          break;
            //       } // if
            //       kbin++;
            //    } // while
            //    assert (kbin < nbins_veg_vert && "computeBottomDragVegetation: veg vert failed");

            //    // for linear interpolation, compute the ratio of the height above the bottom level to the bin width
            //    const double zrat = (tmp_h - veg_vert_z[kbin])/(veg_vert_z[kbin+1] - veg_vert_z[kbin]);
               
            //    // Linear interpolation of the vegetation parameters
            //    cff1 = veg_vert_scale_cd[kbin]*(1.0 - zrat) + veg_vert_scale_cd[kbin+1]*zrat;
            //    cff2 = veg_vert_scale_N [kbin]*(1.0 - zrat) + veg_vert_scale_N [kbin+1]*zrat;
            //    cff3 = veg_vert_scale_D [kbin]*(1.0 - zrat) + veg_vert_scale_D [kbin+1]*zrat;
            // } // if rl10

            // node->veg_alpha3D[k] = node->veg_alpha0*cff1*cff2*cff3;

         } // for k


      } else { // Ganthy (iveg==2)
      //  Modify canopy height, diameter, and density

         // Make sure there is vegetation on this node
         if (node->veg_h_unbent > 0.0) {
            const double wtmp2 = length(node->dav)*100.0; // cm/s

            node->veg_h = 0.72e-2*(4.657 - 0.158*wtmp2 + 0.262*veg_lai - 0.011*veg_lai*wtmp2 + 
                                   0.0022*wtmp2*wtmp2 + 0.048*veg_lai*veg_lai) - 0.00784;
            node->veg_h = max(node->veg_h, 1.0e-2);  //impose min of 1 cm

            const double phi_b = veg_cw*node->veg_di_unbent*node->veg_h_unbent/node->veg_h; // \phi_b
            node->veg_di = (node->veg_di_unbent + phi_b)*0.5;
            const double surpHeight = max(0.0, node->veg_h_unbent - node->veg_h); //surplus height
            node->veg_nv = node->veg_nv_unbent*(1.0 + surpHeight/node->veg_h);
         } // if veg_h_unbent

         for (int k = 0; k < nvrt; k++) {
            node->veg_alpha3D[k] = 0.5*node->veg_di*node->veg_nv*node->veg_Cd;
         } // for k
      

      } // if iveg
     

      // Compute vertical mean of veg_alpha3D
      if (iveg != 0) {
         const double ztop = node->zTopCanopy(); // Top of the canopy
         double sum1=0.0;
         double sum2=0.0;
         // Integrate over each level under the canopy...
         for (int k = kb; (k < nvrt - 1) || (node->z[k + 1] <= ztop); k++) { /// XXX CHECK INDICES
            const double Dz = node->Dz(k);

            sum1 += 0.5*(node->veg_alpha3D[k] + node->veg_alpha3D[k + 1])*Dz;
            sum2 += Dz;

         } // for k
         if (sum2 == 0.0) node->veg_alpha_vert_mean = node->veg_alpha3D[kb];
         else             node->veg_alpha_vert_mean = sum1/sum2;

      } // if iveg

   }

   __global__ void computeVegetationAndDragParameters(void) {
      // Computes the vegetation and drag model parameters:
      // chi, \tilde\chi, beta, \breve{H}
      
      SIDE(side);

      side->hhat = 0.0;
      side->chi2 = 0.0;
      side->chi = 0.0;
      side->veg_c2 = 0.0;
      side->veg_beta = 0.0;

      if (side->isDry) return;

      const double htot = side->htot2();
      const int    kb   = side->kb;

      assert(htot > 0.0 && "computeVegetationModelParameters: htot <= 0");

      const double veg_h     = side->veg_h();
      const double veg_nv    = side->veg_nv();
      const double veg_c     = side->veg_c();
      const double veg_beta  = side->veg_beta;
      const double veg_alpha = side->veg_alpha_vert_mean();

      // Compute the drag coefficients chi2 and chi (\tilde\chi in notes) // XXX recommended to change names, chi->chitilde, chi2->chi
      // Velocity magnitude at the bottom level
      // XXXX Check if we have added -dt*(bcc(1,k,n1)+bcc(1,k,n2))/2.d0 - dt*swild2(k,1) to u2 (I think this only pertains to iupwind_mom \= 0)
      const double vmag1 = length(side->btu[kb + 1]); //Backtracked velocity at the bottom level
      side->chi2 = side->Cd()*vmag1;
      side->chi  = side->chi2;

      if constexpr(iveg != 0) {
         const Node *n1 = side->n1;
         const Node *n2 = side->n2;
         const double veg_alpha_bot = 0.5*(n1->veg_alpha3D[n1->kb] + n2->veg_alpha3D[n2->kb]);
         side->chi = side->chi/(1.0 + veg_alpha_bot*vmag1*dt);
      }

      // Calc consts in SAV model; make sure veg_[c2,beta]=0 at 2D, dry side, and emergent side

      // Compute the vertically averaged velocity magnitude Ua up to the top of the canopy.
      // \int_{-h}^{z_v} |u|dz / H^\alpha = \bar{|u|}^\alpha
      double magUa = 0.0;
      if constexpr(iveg != 0) {
         if (side->is3D() && veg_h > 0.0) {
            // 3D wet side with veg...

            double magUa = 0.0;
            for (int k = kb; k < nvrt - 1; k++) { // XXX CHECK INDICES
               if (side->hasVegetation(k + 1)) { //partial or full SAV layer

                  // Velocity magnitude must be computed at the center of each level to perform the integration.
                  // The tricky part is to compute that velocity at the level that may be partially submerged.
                  const double zrat = side->zLevelVegRatio(k + 1); // in (0,1]

                  // Velocity magnitude at the top (h) of the level (may be partially submerged)
                  const double uh = length((1.0 - zrat)*side->u2[k] + zrat*side->u2[k + 1]);

                  // Velocity magnitude at the bottom (l) of the level
                  const double ul = length(side->u2[k]);

                  const double Dz = zrat*side->Dz(k);  // Takes into account the partially submerged level
                  magUa += Dz*0.5*(uh + ul);

               } // if HasVegetation
            } // for k
            const double ua = magUa/min(htot, veg_h); // Vertical average

            // Computes beta, used for computing the stress term at the top of the canopy.      
            side->veg_c2 = veg_alpha*dt*ua; // >=0
            if (side->isSubmerged()) { // 3D wet submergent side with vegetation
               const double zCpy = side->zTopCanopy();
               const double I = max(1.0e-5, side->chi2*vmag1/(grav*htot)); // Energy gradient
               const double Dz_veg = veg_h - side->Dz(kb + 1); // Distance from the top of the canopy to the top of the bottom level (equal to side->topOfCanopy() - side->zc[kb + 1])
               const double beta2 = sqrt(sqrt(veg_nv)/veg_h)*(-0.32 - 0.85*log10((htot - veg_h)/veg_h*I));
               const double beta = exp(beta2*Dz_veg) - 1.0; /// XXX CHECK INDICES
               side->veg_beta = max(0.0, min(10.0, beta));  // Clamp to [0,10]
            } // if isSubmerged
         } // if is3D
      } // if iveg   

      // Computes the friction modified depth, \breve{H} in notes
      if (side->is2D()) { // 2D

         double fricH = htot + side->chi2*dt; // Depth enhanced by the bottom friction and form drag, \tilde{H} in notes
         if constexpr(iveg != 0) fricH += veg_alpha*vmag1*htot*dt;

         assert(fricH > 0.0 && "COMPUTE CONTINUITY EQUATION: Impossible dry 53");

         side->hhat = htot*htot/fricH; 

      } else { // 3D
         side->hhat = htot - side->chi*dt; // Friction modified depth, \hat{H} in notes  

         if constexpr(iveg != 0) {
            
            if (side->isSubmerged()) { 
               const double Hhat_alpha = veg_h + veg_beta*side->chi*dt; 
               side->hhat -= veg_c*Hhat_alpha;

            } else { // Emergent
               side->hhat = side->hhat/(1.0 + side->veg_c2); // veg_c2 = alpha*|u|*dt
            }
         } // if iveg

         // Enforce positivity for 3D model
         if constexpr(ihhat == 1) side->hhat = max(0.0, side->hhat);
      } // if 3D

   }



   __global__ void getRhoDiffAndVec(const Element *e1, const Element *e2, const int j, const int k, const double rho2, Vector2 &vec, double &diffRho) {
   // Computes the density gradient between two elements, e1 and e2.
   // If the neighbor element (e2doesn't exist (boundary or dry element), the density gradient is assumed to be
   // zero in the direction normal to the side. 
   // 
   // NOTE: In fact it doesn't compute the gradient, but the difference in de-meaned density between e1 and e2 and
   // the vector from e1 to e2 (or normal to side if neighbor is dry or not present
   //
   // Input:
   //    e1: Element
   //    e2: Neighbor element (or nullptr if boundary)
   //    j: Index of the neighbor element (and the side of e1 that is shared with e2)
   //    k: Index of the half-level at which the density is interpolated
   //    rho2: De-meaned density at e2 interpolated at the k-th half-level of e1 (or zero if neighbor is dry or not present)
   // 
   // Output:
   //    vec: Vector from e1 to e2 (or normal to side if neighbor is dry or not present)
   //    diffRho: Difference in de-meaned density between e1 and e2 (or zero if neighbor is dry or not present)

      if (e2 == nullptr || e2->isDry) {
         // Neighbor doesn't exist or is dry, we assume the density gradient is zero in the direction normal to the side.

         Side *side = e1->sides[j];
         assert((side->otherElem(e1) == nullptr || side->otherElem(e1)->isDry) && "getRhoDiffAndVec: For boundary elements, side->otherElem(elem) must be nullptr or dry");
                 
         if constexpr(ics == 1) {
            vec = side->n; // Normal vector to the side
         } else {
            assert(false && "getRhoDiffAndVec: lat/lon not implemented");
            // IMPLEMENT THIS XXXXXX
            // !replace with xel, yel?
            //                   call project_pt('g2l',xnd(n1),ynd(n1),znd(n1),(/xctr(i),yctr(i),zctr(i)/), &
            //      &eframe(:,:,i),xn1,yn1,tmp)
            //                   call project_pt('g2l',xnd(n2),ynd(n2),znd(n2),(/xctr(i),yctr(i),zctr(i)/), &
            //      &eframe(:,:,i),xn2,yn2,tmp)
         } //ics
                   
         diffRho = 0.0;

      } else { //internal and wet
         if constexpr(ics==1) {
            assert(e2 != nullptr || e2->isDry && "getRhoDiffAndVec: neighbor element is nullptr or dry");

            // Vector from e1 to e2
            vec = e2->center - e1->center;

         } else { 
            assert(false && "getRhoDiffAndVec: lat/lon not implemented");
            // Remember that elem->center is in global coordinates
            // IMPLEMENT THIS XXXXXX
            //                   call project_pt('g2l',xctr(ie),yctr(ie),zctr(ie),(/xctr(i),yctr(i),zctr(i)/), &
            //      &eframe(:,:,i),xctr2,yctr2,tmp)
            //                   alow(j)=xctr2
            //                   bdia(j)=yctr2
         } // ics

         diffRho = rho2 - e1->demeanedDensity(k);
      }

   }

   __global__ void computeDensityCubicSplines(void) {
      // Computes cubic splines of de-meaned density (rho - rho_mean) for all elements

      ELEM(elem);

      if (elem->isDry) return;

      const int kb = elem->kb;

      // Computes the cubic splines of de-meaned density (rho - rho_mean) for all elements
      if (elem->is2D()) {
         // In the 2D case, the density is constant at the prism level. 
         // In this case, rho_z doesn't contain the z derivatives, but the actual de-meaned density value.
         elem->rho_z[nvrt-1] = elem->demeanedDensity(nvrt-1);
      } else { //3D

         // Lambda functions to access the function x,y points (x = z at element prism and y = de-meaned density)
         auto getX = [elem] (int k) { return elem->z_prism(k); };
         auto getY = [elem] (int k) { return elem->demeanedDensity(k); };
         ///// XXXX CHECK INDICES!!!!!
         initCubicSpline(kb, nvrt, getX, getY, 0.0, 0.0, elem->rho_spl);
      } // 2D/3D

   }
} // extern "C"


__device__ double interpolateDensityAtHalfLevel(const Element *elem, const Element *elem2, int k) {
// Interpolates the de-meaned density of elem2 at the same z coordinate as the k half-level of elem.

   if (elem2 == nullptr || elem2->isDry) return 0.0;

   const int kb  = elem->kb;
   const int kb2 = elem2->kb;

   // z-coordinate of prism of elem
   const double zp = elem->z_prism(k);

   double rho = 0.0;

   if (elem2->is2D()) {
      return elem->rho_z[nvrt - 1]; // Remember that we stored the de-meaned density in rho_z for the 2D case.

   } else {  // 3D

      // eta_min may be < zmax
      const double eta_min = min(elem->z[nvrt - 1], elem2->z[nvrt - 1]);
      const double zmax    = max(elem->z[kb   + 1], elem2->z[kb   + 1]);

      // (Extrapolate density below bottom, we will fix it later)
      // NOTE: interpolate densities at elem2 at coordinates from elem.
      rho = evalCubicSpline(elem->rho_spl, zp, true);
      // evalCubicSpline(nvrt-kbe(ie),swild2(kbe(ie)+1:nvrt,12), &
      //      &erho(kbe(ie)+1:nvrt,ie)-rho_mean(kbe(ie)+1:nvrt,ie), &
      //      &hp_int(kbe(ie)+1:nvrt,ie,1),nvrt-kbe(i),swild(kbe(i)+1:nvrt),ibot_fl,zmax,eta_min,swild2(kbe(i)+1:nvrt,j))
   } // 2D/3D
     

   // For the case where the elem is deeper than elem2, we extrapolated the density in evalCubicSpline.
   // Now we need to fix it.
   if (elem->z[k] < elem2->z[kb2]) {
      // We are at a z level that doesn't exist at elem2, so we need to adjust the density 

      // There are two methods to adjust the density at the bottom of the element.
      // iunder_deep=0: 
      const double depth1 = abs(elem->z[kb]); // larger depth

      if (iunder_deep==0) {

         const double weight = clamp((depth1 - h1_bcc)/(h2_bcc - h1_bcc), 0.0, 1.0); // weight

         rho = (1.0 - weight)*rho + weight*elem->demeanedDensity(k);

      } else { // iunder_deep != 0

         const double diffDepth = abs(elem->z[kb] - elem2->z[kb2]); // Change??? XXX Was like this in the original code

         if (depth1 >= hw_depth && diffDepth >= hw_ratio*depth1) {

         //   Demeaned density at the same z coordinate as
               rho = elem->demeanedDensity(k); // would make cupp() below=0
         }

      } // if iunder_deep
   } // if elem->z[kb] < elem2->z[kb]

   return rho;
}


extern "C" {
   __global__ void computeDensityGradient(void) {

      assert(ibc == 1 && "baroclinic formulation is not implemented yet: ibc != 0");
      if constexpr(ibc != 0) return;

         // Compute the horizontal density gradient at each element.
      // It works in two steps:
      // 1. Interpolate the density at the neighboring elements at the same zcor as the prism of element i.
      // 2. Compute the density gradient at elem by averaging the density gradients computed between elem and pairs of neighboring elements.
      //    To compute one of these gradients (that it later averages), it solves a 2nd rank linear system like:
      //             gradRho.r1 = rho1 - rho0
      //             gradRho.r2 = rho2 - rho0
      //    where r1, r2 are the vectors from the prism center to the neighboring prism centers.
      // 
      //  NOTE: We use de-meaned densities.

      ELEM(elem);

      if (elem->isDry) return;

      const int kb  = elem->kb;
      const int i34 = elem->i34;

      for (int k = kb + 1; k < nvrt; k++) {
         // Computes the density gradient at elem as average of the density gradients computed between elem and pairs of neighboring elements.
         //             gradRho.r1 = rho1 - rho0
         //             gradRho.r2 = rho2 - rho0
         int icount = 0;
         double sumRx = 0.0;
         double sumRy = 0.0;
         for (int j1=0; j1<i34; j1++) { 

            const int j2 = elem->nxq(j1 + 1);

            // Pairs of consecutive neighbor elements of elem
            const Element *e1 = elem->nghbElems[j1];
            const Element *e2 = elem->nghbElems[j2];

            // Using cubic splines, interpolate demeaned density difference at the same z half-levels as elem.
            const double rho1 = interpolateDensityAtHalfLevel(elem, e1, k);
            const double rho2 = interpolateDensityAtHalfLevel(elem, e2, k);

            // Compute the density difference and the vector from elem to e1 and e2
            // getRhoDiffAndVec takes into account boundaries/dry elements
            double diffRho1, diffRho2;
            Vector2 vec1, vec2;
            getRhoDiffAndVec(elem, e1, j1, k, rho1, vec1, diffRho1); // We could remove one of these calls, but it would make the code uglier.
            getRhoDiffAndVec(elem, e2, j2, k, rho2, vec2, diffRho2);

            const double rl1 = length(vec1);
            const double rl2 = length(vec2);
            assert(rl1 > 0.0 && rl2 > 0.0 && "computeDensityInterpolation: one of rl1 or rl2 is 0");

            const double delta = cross(vec1, vec2);
            const double sintheta = abs(delta)/(rl1*rl2);

            // To avoid problems, avoid elements with small angle between vectors
            constexpr double sin_one_deg = 0.01745240643728351; // sin(pi/180.0). We use 1 degree as threshold for well-conditioned system
            const bool validEq = (sintheta > sin_one_deg);

            if (validEq) {
               // Solve the system of equations to get the density gradient at elem
               icount++;
               const double Rx = (vec2.y*diffRho1 - vec1.y*diffRho2)/delta;
               const double Ry = (vec1.x*diffRho2 - vec2.x*diffRho1)/delta;
               sumRx += Rx;
               sumRy += Ry;
            }
         } // for j1

         // Compute the average density gradient at elem
         if (icount==0) {
           // Bad mesh quality; bail out with 0 b-cc
           elem->gradRho[k] = {0.0, 0.0};
         } else {
           elem->gradRho[k] = {sumRx/icount, sumRy/icount};
         }

      } // for k

   }

   __global__ void computeBaroclinicTerm(void) {
   // Computes the baroclinic term at sides (at half levels).
   // For this it has to compute the density gradient at the half-sides (face centers) and then integrate it vertically.
   // The baroclinic term is -(g/rho_0)*int_z^\eta grad(rho) dz

      if (ibc == 1) return;

      SIDE(side);

      if (side->isDry) return;

      const Element *e1 = side->elems[0];
      const Element *e2 = side->elems[1]; 
      const int kb = side->kb;
      
      assert(e1 != nullptr || e2 != nullptr && "computeSideDensityGradient: side needs to have at least one element");

      // interpolates vertically the density gradient of the two elements at the half level/

      int k1 = max(kb, nvrt / 2);  // Start the search for the half level at the bottom of the side
      int k2 = max(kb, nvrt / 2);

      // Starts the integration at the surface and goes down to the bottom
      side->bcc[nvrt-1] = {0.0, 0.0};
      for (int k = nvrt - 2; k >= kb; k--) {
      // For each prism level (from surface to bottom)...

         // For both elements, interpolate density gradient at the side's half level.
         // The gradient at the side is the average of the gradients at the elements.
         // w1, w2 are used to take into account boundary sides, with just one element.
         double w1 = (e1!=nullptr) ? 1.0 : 0.0;
         double w2 = (e2!=nullptr) ? 1.0 : 0.0;
         const double wsum = w1 + w2;
         w1 = w1/wsum;
         w2 = w2/wsum;
         assert(w1 + w2 > 0.99 && w1 + w2 < 1.01 && "computeSideDensityGradient: w1 + w2 is not 1");

         // Linear interpolation of both elements' density gradient at the shared side's level.
         Vector2 gradRho = {0.0, 0.0};
         if (e1 != nullptr && k > kb) {
            // XXXX WHAT HAPPENS IF k1 < kb??? (check ibelow)
            
            k1 = e1->findLevel(side->zc[k], k1); // finds the index at the half level starting the search at k1
            const double zrat = (side->zc[k] - e1->z[k1])/(e1->z[k1 + 1] - e1->z[k1]);
            assert(zrat >= 0.0 && zrat <= 1.0 && "computeSideDensityGradient: zrat is not in [0, 1]");

            // Grad(rho) vertically interpolated at the side's half level
            Vector2 gradRho1 = e1->gradRho[k1]*(1.0 - zrat) + e1->gradRho[k1 + 1]*zrat;

            // Projects the element density gradient to the side frame if needed.
            if constexpr(ics==2) gradRho1 = e1->frame.projHorVec(gradRho1, side->frame);

            // Adds the contribution of the density gradient at the side's half level
            gradRho += gradRho1*w1;
         }      
         if (e2 != nullptr && k > kb) {
            k2 = e2->findLevel(side->zc[k], k2); // finds the index at the half level starting the search at k2
            const double zrat = (side->zc[k] - e2->z[k2])/(e2->z[k2 + 1] - e2->z[k2]);
            assert(zrat >= 0.0 && zrat <= 1.0 && "computeSideDensityGradient: zrat is not in [0, 1]");
            
            // Grad(rho) vertically interpolated at the side's half level
            Vector2 gradRho2 = e2->gradRho[k2]*(1.0 - zrat) + e2->gradRho[k2 + 1]*zrat;

            // Projects the element density gradient to the side frame if needed.
            if constexpr(ics==2) gradRho2 = e2->frame.projHorVec(gradRho2, side->frame);

            // Adds the contribution of the density gradient at the side's half level
            gradRho += gradRho2*w2;
         }

         
         // Integrate (using a central Gauss point) the density gradient from the surface to the bottom 
         // multiplied by -rampbc*(g/rho0) 
         const double grav3 = side->grav3();
         side->bcc[k] = side->bcc[k + 1] - (rampbc*grav3/rho0)*side->Dz(k + 1)*gradRho;  

      } // for

   }



   __global__ void computeHorGradientAtSide(void) {

      // Compute the horizontal density gradient at each element.
      // It works in two steps:
      // 1. Interpolate the density at the neighboring elements at the same zcor as the prism of element i.
      // 2. Compute the density gradient at elem by averaging the density gradients computed between elem and pairs of neighboring elements.
      //    To compute one of these gradients (that it later averages), it solves a 2nd rank linear system like:
      //             gradRho.r1 = rho1 - rho0
      //             gradRho.r2 = rho2 - rho0
      //    where r1, r2 are the vectors from the prism center to the neighboring prism centers.
      // 
      //  NOTE: We use de-meaned densities.

      SIDE(side);


      if (side->isDry) return;

      // Computes the gradient of the variables at the side center.
      auto isWet1 = [](const Element* e) { return !(e->isDry); };
      auto isWet2 = [](const Element* e) { return !(e->isDry2()); };
      auto isWet3 = [](const Element* e) { return (e->depth >= tip_dp) && !(e->isDry); };

      side->gradEta1 = side->grad2(&Node::eta1, isWet1);
      side->gradEta2 = side->grad2(&Node::eta2, isWet2);
      side->gradETP  = side->grad2(&Node::ETP,  isWet3);
      side->gradPr   = side->grad2(&Node::pr,   isWet1);

   }


   __global__ void solveMomentumEquations(double time) {

      double alow[nvrt], cupp[nvrt], bdia[nvrt];
      Vector2 rrhs[nvrt], u[nvrt], soln[nvrt], gam[nvrt];

      SIDE(side);

      if (side->isDry) 
      {
         for (int k=0; k<nvrt; k++) side->u2[k] = {0.0, 0.0};       // WHY {0,0}? Shouldn't be just the normal component?  The FORTRAN code also uses {0,0}  in line ~6462   

         return;
      }

      const Node *n1 = side->n1;
      const Node *n2 = side->n2;

      const int kb = side->kb;

      const double veg_h = side->veg_h();
      const double grav3 = side->grav3();

      if (side->is2D()) {
         // Remember: Do not use htot2, which is updated.
         const double htot = side->htot1();
         const double hhat = side->hhat;

         assert(hhat > 0.0 && htot > 0.0 && "solveMomentumEquations: Impossible");

         // tau() takes care of the projection to the side frame (ics==1 or ics==2)
         const Vector2 tau2 = side->tau();

         // hat_gam has units of m/s

         const Vector2 coriolisTerm = side->coriolis*rot90(side->u2[nvrt-1]);
         Vector2 hat_gam = side->btu[nvrt-1] + dt*(coriolisTerm - side->gradPr/rho0 + grav3*side->gradETP +
                           side->bcc[kb] + tau2/htot - grav3*((1 - thetai)*side->gradEta1 + thetai*side->gradEta2));

         hat_gam += dt*side->nuGradU[nvrt-1];  // XXXX compute nuGradU


         const Vector2 vel = clampVelocity(hat_gam*hhat/htot);

         for (int k = 0; k < nvrt; k++) side->u2[k] = vel;

      }
      else { // 3D

         // For 3D it is necessary to solve an implicit tridiagonal system.
         // 
         // Each line of the tridiagonal system include variables for layers k-1, k, k+1.
         // However, for the vegetation terms, the canopy complicates things, because the vertical 
         // integration is only done up to the canopy. In this case we use the notation "t" (top)
         // to indicate the "top" level, which may not coincide with a layer. "t" can stand for both
         // k+1 and k.


         const int ndim = nvrt - kb;
         const double zTopCanopy = side->zTopCanopy();

         // Initialize the tridiagonal system
         for (int k = kb + 1; k < nvrt; k++) {
            
            const int kin = k - (kb + 1);

            // Vertical viscosity at the face center
            const double dfz = 0.5*(side->dfv(k - 1) + side->dfv(k)); // Vertical average.

 
            alow[kin] = 0.0;
            cupp[kin] = 0.0;
            bdia[kin] = 0.0;

            double cff1 = 1.0; //init \bar{c^{k+1}}              
            // Vegetation term: alpha*|u^n|*dt*H(dv-d)*u^{n+1}

            // We divide the code into two parts, one (L) that deals with layers k, k+1, and another (R) that deals with k-1, k.
            // There is a special case when k is the bottom layer, where the bottom drag is applied.
            // Notice that the top boundary is not included implicitly, because it depends on tau, but not on u.

            // First part: layers k, k+1
            if (k < nvrt - 1) {
                             
               if constexpr(iveg != 0) {

                  // If there is vegetation at this level...
                  if (zTopCanopy > side->zc[k]) {
                     // If the canopy is above the k layer...
                     
                     // Computes the velocity at the level above (typically layer k+1, unless the canopy is below it)
                     const double zt = min(zTopCanopy, side->zc[k + 1]);
                     const double zrat = (zt - side->zc[k])/side->Dz(k + 1); // value in (0,1], 1 means that the canopy is above zt.
                     const Vector2 ut2 = (1.0 - zrat)*side->u2[k] + zrat*side->u2[k + 1];

                     // Computes the velocity magnitude at the face center (the section of the face covered by the vegetation)
                     const double vmag2 = length(ut2);
                     const double vmag1 = length(side->u2[k]);
                     const double vmag = 0.5*(vmag1 + vmag2);
                     
                     // Computes the vegetation alpha at the face center
                     const double veg_alpha_sd = 0.5*(side->veg_alpha3D(k) + side->veg_alpha3D(k+1));

                     cff1 += veg_alpha_sd*dt*vmag;

                  } // if (zTopCanopy
               } // if iveg

               const double dzz = side->Dz(k + 1);
               const double tmp = dt*dfz/dzz;
               assert(dzz > 0.0 && "STEP: dzz<=0 in momentum");

               cupp[kin] += cff1*dzz/6.0 - tmp;
               bdia[kin] += cff1*dzz/3.0 + tmp;
            }


            // Second part: layers k-1, k
            if (k > kb + 1) {

               double cff1 = 1.0; //init \bar{c^k}

               if constexpr(iveg != 0) {
                  // If there is vegetation at this level...
                  if (zTopCanopy > side->zc[k-1]) {
                     const double zz1 = min(zTopCanopy, side->zc[k]);
                     const double zrat = (zz1 - side->zc[k-1])/side->Dz(k); // value in (0,1]

                     const Vector2 ut2 = (1.0 - zrat)*side->u2[k-1] + zrat*side->u2[k];

                     const double vmag2 = length(ut2);
                     const double vmag1 = length(side->u2[k-1]);
                     const double vmag = 0.5*(vmag1 + vmag2);

                     // Computes the vegetation alpha at the face center
                     const double veg_alpha_sd = 0.5*(side->veg_alpha3D(k) + side->veg_alpha3D(k-1));
                     cff1 += veg_alpha_sd*dt*vmag;
                  } // if zTopCanopy

               } // if iveg

               const double dzz = side->Dz(k);
               const double tmp = dt*dfz/dzz;
               assert(dzz > 0.0 && "STEP: dzz<=0 in momentum");

               alow[kin] += cff1*dzz/6.0 - tmp;
               bdia[kin] += cff1*dzz/3.0 + tmp;

            } else { // Bottom boundary
               // This occurs at the bottom, where the bottom drag is applied.

               bdia[kin] += dt*side->chi2; // XXXX

            } // if (k > kb + 1)


         } // for (int k



         // Computes the RHS.
         // -----------------

         //	BC to be imposed at the end.

         const Vector2 gradETP  = side->gradETP;
         const Vector2 gradPr   = side->gradPr;
         const Vector2 gradEta1 = side->gradEta1;
         const Vector2 gradEta2 = side->gradEta2;

         const Vector2 *btu = side->btu;
         const Vector2 *u2  = side->u2;
         const double coriolis = side->coriolis;
         const Vector2 *nuGradU = side->nuGradU;

         const Vector2 grad2D = grav3*(thetai*gradEta2 + (1.0-thetai)*gradEta1 - gradETP) + gradPr/rho0;
         for (int k = kb + 1; k < nvrt; k++) {
            const int kin = k - kb - 1;

            const double dz_0 = side->Dz(k);
            const double dz_1 = (k < nvrt - 1) ? side->Dz(k + 1) : 0.0;  // k always < nvrt - 1 when used.

            rrhs[kin] = {0.0, 0.0};
            double dz = 0.0;

            // dz is the distance between the center of the k-th layer and the center of the (k+1)-th layer
            // In the boundary layers, it is the distance between the center of the layer and the boundary.
            if (k < nvrt  ) dz += dz_1;
            if (k > kb + 1) dz += dz_0;

            // Elevation gradient, atmospheric pressure and tidal potential
            rrhs[kin] = - (0.5*dz)*dt*grad2D;

            // Coriolis, advection, wind stress, and horizontal viscosity
            if (k < nvrt - 1) {


               rrhs[kin] += (dz_1/6.0)*(2.0*btu[k] + btu[k + 1] + 
                                        dt*(coriolis*rot90(2.0*u2[k] + u2[k + 1]) + 
                                            2.0*nuGradU[k] + nuGradU[k + 1]));
                                                

               //diff stress tensors 1006
               if constexpr(itur==5) {  // 1018:itur==5 XXX What is itur?
                  // Does it make sense to have dz_1 multiplying and dividing? XXX
                  assert(false && "Not implemented: itur == 5");
                  // rrhs[kin] += 0.5*dt*dz_1*(side->TDz[k+1] - side->TDz[k])/dz_1
               }
               // else assert(false && "itur # Not implemented");

            } else { // k == nvrt - 1 SURFACE
               // tau() takes care of the projection to the side frame (ics==1 or ics==2)
               const Vector2 tau2 = side->tau();

               rrhs[kin] += dt*tau2;          

            } // if k<nvrt


            if (k > kb + 1) {
               rrhs[kin] += (dz_0/6.0)*(2.0*btu[k] + btu[k - 1] + 
                                       dt*(coriolis*rot90(2.0*u2[k] + u2[k - 1]) + 
                                           2.0*nuGradU[k] + nuGradU[k - 1]));

               // diff stress tensors 1006
               if constexpr(itur==5) { // 1018:itur==5 XXX What is itur?
                  assert(false && "Not implemented: itur == 5");
                  // Does it make sense to have dz_0 multiplying and dividing? XXXX
                  // rrhs[kin] += 0.5*dt*dz_0*(side->TDz[k] - side->TDz[k-1])/dz_0;
               }
               // else assert(false && "itur # Not implemented");


            } // if k>kbs + 1

            // Baroclinic
            Vector2 *bcc = side->bcc;
            if constexpr(ibc == 0) {
               if (k < nvrt  ) rrhs[kin] += (dz_1*dt/6.0)*(2.0*bcc[k] + bcc[k+1]);
               if (k > kb + 1) rrhs[kin] += (dz_0*dt/6.0)*(2.0*bcc[k] + bcc[k-1]);

            } // if ibc==0

         } // k=kb + 1,nvrt


         tridag_sch<2>(nvrt - kb - 1, alow, bdia, cupp, (double *)rrhs, (double *)soln, (double *)gam);

         for (int k = kb; k < nvrt; k++) {  // XXXXX YYYYYYYYY can be kb 
            const int kin = k - kb - 1;
   
            // Velocity is always clamped to maximum values
            side->u2[k] = clampVelocity(soln[kin]);
         } // k

      } // if side->is2D()
      

      // Impose velocity boundary condition at the bottom
      if constexpr(imm == 2) { // 3D bed deformation model
         double dep;
         Vector2 vel;

         assert(false && "Not implemented: imm == 2");
         //update_bdef(time, x0, y0, dep, vel);
         // side->u2[kb] = vel;

      } else { // imm = 0 or 1
         if (side->Cd() == 0.0) {
            side->u2[kb] = side->u2[kb + 1];

         } else if (side->is3D()) { // 3D no slip bottom   
            side->u2[kb] = {0.0, 0.0};

         }
          
      } // if imm

      // The velocities below the bottom level are set to zero
      for (int k = 0; k < kb; k++) side->u2[k] = {0.0, 0.0};


   }


   __global__ void applyMomentumBC(void) {

      SIDE_BND(sideBnd);

      Side *side = sideBnd->side;
      const BoundarySegment *segment = sideBnd->segment;
      const Node *n1 = side->n1;
      const Node *n2 = side->n2;

      assert(segment != nullptr && "applyMomentumBC: segment is nullptr");

      const bool isLand = segment->isLand;
      const int ifltype = segment->ifltype;

      const int kb = side->kb;


      if (!isLand && ifltype != 0) {
         const Vector2 *ubc = sideBnd->ubc;
  
         for (int k = kb; k < nvrt; k++) { // YYYYYYYYY can be kb + 1  XXXX
            
            assert(ifltype == -2 || ((ubc[k].x >= -98 && ubc[k].y >= -98) 
                   && (ubc[k].x <= 98 && ubc[k].y <= 98)) && "Momentum BC: Side has no velocity boundary conditions");

            if (ifltype == -1) {
               // Flather 1

               const double depth = side->depth();

               assert(depth > 0.0 && "Flather bnd has 0 or negative depth");

               const double vnorm = sqrt(grav/depth)*(side->eta2() - sideBnd->etaMeanBnd()) + dot(ubc[k], side->n);

               side->u2[k] = vnorm*side->n;

            } else if (ifltype == -2) { // Discharge

               const double eta1 = side->eta1();
               const double eta2 = side->eta2();
               const double htot = side->htot2();

               // Discharge is a polynomial. See disch_coef in readTides. We must implement this.
               assert(false && "Not implemented: Discharge bc");
               // const double dischrg = discharge(eta1);
               const double dischrg = 0.0; // XXX REMOVE
               

               assert(htot > 0.0 && "Discharge bc has 0 or negative htot");

               const double vnorm = dischrg*eta2/htot;

               side->u2[k] = vnorm*side->n;

            }
            else if (ifltype == -4 || ifltype == -5) { // 3D radiation
               const double vnorm = dot(side->u2[k], side->n);

               double rat;
               if (vnorm <= 0.0) rat = segment->vobc1; // incoming
               else              rat = segment->vobc2; // outgoing

               side->u2[k] = (1.0 - rat)*side->u2[k] + rat*ubc[k];

            }
            else { // not Flather or 3D radiation
               side->u2[k] = ubc[k];

            } // if (segment->ifltype...

         }

      }
      else if (isLand) { // if it is land boundary side

         if constexpr(islip == 0) { // free slip

            for (int k = kb; k < nvrt; k++) {

               const Vector2 t = side->t();
               const Vector2 n = side->n;

               double vnorm = 0.0; // for most cases
               const double vtan = dot(side->u2[k], t);
               
               side->u2[k] = vtan*t - vnorm*n;

            }

         } else { // no slip
            for (int k=kb; k < nvrt; k++) side->u2[k] = {0.0, 0.0};

         }
      } // if (side segment type...

   }


   __global__ void shapiroFilter1(void) {
      // First step of the Shapiro filter (normally used if indvel<=0)
      // Computes the filtered values and stores them in bcc (as a temporary array).

      SIDE(side);

      // Only internal wet sides are processed
      if (side->isDry || !(side->isInternal())) return;

      const int kb = side->kb;

      for (int mm = 0; mm < niter_shap; mm++) {

         for (int k = kb + 1; k < nvrt; k++) {
            // WARNING: reusing bcc as a temporary array to store the suru values.
            Vector2 suru = {0.0, 0.0};

            for (int j = 0; j < 4; j++) { // 4 neighbors per side

               const Side *nSide = side->neigh(j); // Each internal side has 4 neighbors. They are the closest sides to each of the two elements that share this side

               int kin = k;
               if (!(nSide->isDry)) kin = max(k, nSide->kb + 1);
               
               // Get the velocity at the neighbor's side.
               Vector2 nSideU2 = nSide->u2[kin];

               // Projects the velocity to the current side frame if needed.
               if constexpr(ics==2) nSideU2 = nSide->frame.projHorVec(nSideU2, side->frame);

               suru += nSideU2;
            }

            // Reusing bcc as a temporary array to store the suru values.
            side->bcc[k] = side->u2[k] + side->shapiro*(0.25*suru - side->u2[k]); // All in the side frame if ics=2

         }
      }
   }

   __global__ void shapiroFilter2(void) {
      // Second step of the Shapiro filter. Assigns to u2 the values stored in bcc.

      SIDE(side);


      if (side->isDry || !(side->isInternal())) return;

      const int kb = side->kb;

      for (int k = 0; k < kb - 1; k++) side->u2[k] = {0.0, 0.0};

      if (side->is2D()) {
         side->u2[kb + 1] = side->bcc[kb + 1];
         side->u2[kb    ] = side->u2 [nvrt - 1];  // XXXXX isn't side->u2 [nvrt - 1] == side->bcc[kb + 1]?

      } else {
         // 3D
         for (int k = kb + 1; k < nvrt; k++) side->u2[k] = side->bcc[k];
         
      }

   }


   __global__ void ul0UpdateZCoords(void) {
   // For wet nodes updates the nodes z-coordinates and bottom index, kb.
      
      NODE(node);      


      if (false) { //node->isDry) { XXXXXXX

         assert((ivcor != 2) || (node->depth < mesh->h_s) && "Deep depth dry");

         if constexpr(ivcor != 1) node->kb = 0;

      } else {

         // Compute the updated node's kb index. This must be done before compZCoord() is called.
         node->updateKb();

         for (int k = node->kb; k < nvrt; k++) node->z[k] = node->compZCoord(k); 

      }

      // We use the temp.b field to store the previous state of isDry.
      node->tmp.b = node->isDry;

      // In the next step we will compute which nodes are wet, so we first assume that all nodes are dry.
      node->isDry = true;

   }

   __global__ void ul0ComputeDryElements(void)
   // Updates isDry based on eta2 for the elements. Also computes kb and z.
   {
      ELEM(elem);

      elem->isDry = elem->isDry2(); // A element is dry if any of its nodes are dry (based on current eta2, not on isDry)

      elem->kb = 0;

      if (elem->isDry) return;


   //    if(elem->isDry2()) { // XXXX
   //        write(errmsg,*)'level0: Element-node inconsistency (0):',ielg(i),idry_e2(i), &
   //   &iplg(elnode(1:i34(i),i)),idry2(elnode(1:i34(i),i)),idry(elnode(1:i34(i),i))
   //        call parallel_abort(errmsg)
   //    }

      const Node *n1 = elem->nodes[0];
      const Node *n2 = elem->nodes[1];
      const Node *n3 = elem->nodes[2];  
      const Node *n4 = elem->nodes[3]; // nullptr for triangles      

      // Computes element's kb as the minimum of the nodes' kb.
      const int kb1 = n1->kb;
      const int kb2 = n2->kb;
      const int kb3 = n3->kb;
      const int kb4 = (elem->isQuad()) ? n4->kb : INT_MAX;  // INT_MAX is never used (see below).


      elem->kb = min(min(kb1, kb2), min(kb3, kb4));


      // Computes the element z as the mean of the nodes z.
      for (int k = elem->kb; k < nvrt; k++) {

         // Compute the element z as the mean of the nodes z.
         elem->z[k] = n1->z[max(k, kb1)] + n2->z[max(k, kb2)] + n3->z[max(k, kb3)];
         if (elem->isQuad()) elem->z[k] += n4->z[max(k, kb4)];          
         elem->z[k] /= elem->i34;

         assert(k == elem->kb || (elem->z[k] > elem->z[k-1]) && "ivcor = 1: z must increase monotonically");
      }
                  
   }

   __global__ void ul0ComputeDryNodes(void)
   // Updates node's isDry based on the elements' isDry. 
   {
      NODE(node);

      node->isDry = true;

      // A node is wet if at least one of the elements sharing it is wet.
      for (int j = 0; j < node->numElements; j++) {
         const Element *elem = node->elements[j];
         if (!elem->isDry) {
            node->isDry = false;
            break;
         }
      }

       
   }


   __global__ void ul0UpdateRewettedNodes(int it, int iths)
   // Obtains values of velocity, S, T for newly wet nodes. It uses the average of the previously wet surrounding nodes.
   // node->tmp.b is true if the node was previously dry.
   // This function avoids race conditions by only modifying nodes that just became wet and only using the values of the previously wet surrounding nodes.
   {
      NODE(node);
      // Compute velocity, S, T for re-wetted nodes (q2 and xl are fine)

      // Avoid this the first time step.
      if (it == iths) return;  

      const bool wasDry = node->tmp.b;

//  if (node->isDry && !wasDry) node->eta2 = max(node->eta2, -node->depth); // XXXX KKKKKK NOT IN THE ORIGINAL CODE

      if ((!node->isDry) && wasDry) {  // rewetted, remember that tmp.b is true if the node was previously dry.

         for (int k=0; k < nvrt; k++) {

            Vector3 utmp={0.0, 0.0, 0.0};

            double ttmp=0.0;
            double stmp=0.0;
            int icount=0;
            for (int j = 0; j < node->numNodes; j++) {

               const Node *nbrNode = node->nodes[j];
               
               const bool wasNbrDry = nbrNode->tmp.b;

               // Wet neighbors not affected by this part and so each sub-domain should use same values
               if (!wasNbrDry) { // Incluide only previously wet surrounding nodes

                  icount++;

                  // Assume small element size in wet/dry zone so pframes are close to each other
                  utmp += nbrNode->uvw2(k);
                  
                  ttmp += nbrNode->tracers[k].tracer[0]; 
                  stmp += nbrNode->tracers[k].tracer[1]; 
               }
            }

            if (icount > 0) {
               const double inv_icount = 1.0/icount;
               utmp = inv_icount*utmp;
               node->u2[k]     = {utmp.x, utmp.y};
               node->w2[k]     = utmp.z;
               node->tracers[k].tracer[0]   = inv_icount*ttmp;
               node->tracers[k].tracer[1]   = inv_icount*stmp;
            } else {
               node->u2[k] = {0.0, 0.0};
               node->w2[k] = 0.0;
               // Initialize tracers to 0.0??? XXXXX
            }
         }
      }
   }

  
   


   __global__ void ul0UpdateSides(int it, int iths) {
   // Obtains values of velocity for newly wet sides. It uses the average of the previously wet surrounding sides.
   // Notice that, since we are only updating newly wetted side from sides that were previously wet, there is no race condition here.
   // The dryness of the sides is updated in ul0UpdateSides2.

      SIDE(side);

      // A side is wet if and only if at least one of its elements is wet.
      const bool isWet = !side->elems[0]->isDry || (side->isInternal() && !side->elems[1]->isDry);
 
      const Node *n1 = side->n1;
      const Node *n2 = side->n2;
     
      const int kb1 = n1->kb;
      const int kb2 = n2->kb;

      side->kb = 0; // default when dry

      if (isWet) {

         assert(~(n1->isDry || n2->isDry) && "Side-node inconsistency. At least one node of a wet side is dry.");

         assert(side->htot2() > h0 && "Weird side. Total depth is smaller than h0.");

         side->kb = min(kb1, kb2);

         for (int k = 0; k < side->kb; k++) side->zc[k] = 0.0;  // XXXX NECESSARY???

         for (int k = side->kb; k < nvrt; k++) {
            side->zc[k] = 0.5*(n1->z[max(k, kb1)] + n2->z[max(k, kb2)]);

            assert(((k == side->kb ) || (side->zc[k] > side->zc[k-1])) && "Weird side (3)");

            // Compute velocity for re-wetted sides 
            if ((it != iths) && (side->isDry && isWet)) {

               Vector2 utmp={0.0, 0.0};

               int icount=0;

               // For all the sides on the two elements that share this side
               for (int j = 0; j < 2; j++) {
                  const Element *elem = side->elems[j];
                  if (elem == nullptr) continue;

                  for (int jj=0; jj < elem->i34; jj++) {
                     const Side *nSide = elem->sides[jj];
                     if (nSide->isDry) continue;

                     icount++;
                     utmp += nSide->u2[k];  
                     //  XXXXXX DO THE SAME FOR TRACERS????                   
                  }
               }
                  
               // assert(icount > 0 && "No wet sides around this side");    XXXX KKKKK
               if (icount > 0)  side->u2[k] = utmp/icount;
               else             side->u2[k] = {0.0, 0.0};
                  
            } // if ((it != iths) && (side->isDry && isWet))

         } // for k

      } // if isWet 
      


   }

   __global__ void ul0UpdateSides2(void) {
      // Obtains values of velocity for newly wet sides. It uses the average of the previously wet surrounding nodes.
   
         SIDE(side);
   
         const Element *e1 = side->elems[0];
         const Element *e2 = side->elems[1];

         // A side is wet if and only if at least one of its elements is wet.
         const bool isWet = (!e1->isDry) || (side->isInternal() && (!e2->isDry));

         // Update the side's dryness
         side->isDry = !isWet;

   }


   __global__ void zeroNodalVelocities(void) {
      NODE(node);

      for (int k = 0; k < nvrt; k++) {
         node->u2[k]   = {0.0, 0.0};
         node->w2[k]   = 0.0;
         node->utmp[k] = {0.0, 0.0};         
      }
      node->tmp.d = 0.0;  // reset the temporary real value. It will be used to compute the area-averaged vertical velocity
   }
 

   __global__ void computeNodalVerticalVelocities1(void) {
   // Compute node's vertical velocities from element's vertical velocities at WHOLE levels.
   // Adds the element's vertical velocity times the element's area.
      ELEM(elem);

      for (int j = 0; j < elem->i34; j++) {
         Node *node = elem->nodes[j];

         for (int k = node->kb; k < nvrt; k++) {
            atomicAdd(&(node->w2[k]), elem->w[k]*elem->area);
         }
         atomicAdd(&(node->tmp.d), elem->area);
      }
   }

   __global__ void computeNodalVerticalVelocities2(void) {
   // Compute node's vertical velocities from element's vertical velocities at WHOLE levels.
   // Divides by the total area of the elements sharing the    node.
      NODE(node);


      const double invArea = 1.0 / node->tmp.d;
      for (int k = node->kb; k < nvrt; k++) node->w2[k] = invArea*node->w2[k];  // node->tmp.d contains the total area of the elements sharing the node

      node->tmp.i = 0; // Reset the temporary integer value. It will be used to compute the number of elements sharing the node.
   }


   __global__ void computeNodalHorizontalVelocities1(void) {

      assert(indvel == 0 && "computeNodalVelocities only works for indvel <= 0. IMPLEMENT THIS FOR INDVEL = 1");

      ELEM(elem);

      if (elem->isDry) return;

      const int kb = elem->kb;
      const bool isTri = elem->isTri();

      if constexpr(indvel <= 0) {
         // Assuming discontinous hvel.

         if (isTri) { // For triangles
         // Interpolates the side velocities at the nodes. A magnitude defined at the sides and interpolated (extrapolated) at one node can be
         // computed as the sum of the adjacent side magnitudes minus the side magnitude that is not adjacent to the node.
            #pragma unroll
            for (int j = 0; j < 3; j++) {
               // For every node of the triangle...
               Node *node = elem->nodes[j];

               if (node->isDry) continue;

               const Side *s1 = elem->sides[j];  // Side opposite to node j.
               const Side *s2 = elem->sides[(j + 1) % 3];
               const Side *s3 = elem->sides[(j + 2) % 3];

               assert(s1->n1 != node && s1->n2 != node && "s1 must not be connected to node");

               atomicAdd(&(node->tmp.i), 1);

               for (int k = kb; k < nvrt; k++)  {
                  Vector2 s1u2 = s1->u2[k];
                  Vector2 s2u2 = s2->u2[k];
                  Vector2 s3u2 = s3->u2[k];

                  // If ics == 2, project the velocity 2D vectors from the side frame to the element frame
                  if constexpr(ics == 2) {
                     s1u2 = s1->frame.projHorVec(s1u2, elem->frame);
                     s2u2 = s2->frame.projHorVec(s2u2, elem->frame);
                     s3u2 = s3->frame.projHorVec(s3u2, elem->frame);
                  }

                  Vector2 u = s2u2 + s3u2 - s1u2; // u extrapolated linearly at node j.   
                  u = clampVelocity(u);
                  atomicAdd(&(node->utmp[k].x), u.x);
                  atomicAdd(&(node->utmp[k].y), u.y);

               }

            }
         }
         else { // For quads


            // Compute the bilinear coordinates of the center of each diagonal.
            double xi, eta;
            double shp13[4], shp24[4];
            Point2 x0, x1, x2, x3;
            elem->getLocalHorizontalNodePositionsQuad(x0, x1, x2, x3);
            const Point2 p13 = 0.5 * (x0 + x2);   // Center of sides 1 & 3
            const Point2 p24 = 0.5 * (x1 + x3);   // Center of sides 2 & 4

            elem->inverseBilinearQuad(p13, xi, eta, shp13);
            elem->inverseBilinearQuad(p24, xi, eta, shp24);




            for (int k = kb; k < nvrt; k++) {

               //Each quad is split into two triangles. For each triangle we extrapolate as we did for triangles.
               Vector2 sideU[4];
               #pragma unroll
               for (int j = 0; j < 4; j++) {
                  Vector2 u = elem->sides[j]->u2[k];
                  if constexpr(ics == 2) u = elem->sides[j]->frame.projHorVec(u, elem->frame);
                  sideU[j] = u;
               }

               Vector2 ud1 = shp13[0]*sideU[0] + shp13[1]*sideU[1] + shp13[2]*sideU[2] + shp13[3]*sideU[3];
               Vector2 ud2 = shp24[0]*sideU[0] + shp24[1]*sideU[1] + shp24[2]*sideU[2] + shp24[3]*sideU[3];

               Vector2 ufg0 = clampVelocity(sideU[2] + sideU[3] - ud2);  // node 1: sides 3,4
               Vector2 ufg1 = clampVelocity(sideU[0] + sideU[3] - ud1);  // node 2: sides 1,4
               Vector2 ufg2 = clampVelocity(sideU[0] + sideU[1] - ud2);  // node 3: sides 1,2
               Vector2 ufg3 = clampVelocity(sideU[1] + sideU[2] - ud1);  // node 4: sides 2,3
               const Vector2 ufg[4] = {ufg0, ufg1, ufg2, ufg3};
               #pragma unroll
               for (int j = 0; j < 4; j++) {
                  Node *node = elem->nodes[j];
                  if (node->isDry) continue;

                  if (k == kb) atomicAdd(&(node->tmp.i), 1);
                  atomicAdd(&(node->utmp[k].x), ufg[j].x);
                  atomicAdd(&(node->utmp[k].y), ufg[j].y);
               }
            }

         }

      }
      else assert(false && "computeNodalHorizontalVelocities1: not implemented for indvel > 0");
      
   }


   __global__ void computeNodalHorizontalVelocities2(void) {
      NODE(node);

      const int nodeCount = node->tmp.i;

      if (nodeCount == 0) 
      {
         for (int k = 0; k < nvrt; k++) node->u2[k] = {0.0, 0.0};
         return;
      }

      assert(!node->isDry && "computeNodalHorizontalVelocities2: node is dry");
      
      const int kb = node->kb;

      for (int k = 0;  k < kb;   k++) node->u2[k] = {0.0, 0.0};
      for (int k = kb; k < nvrt; k++) node->u2[k] = node->utmp[k]*(1.0/nodeCount);
      
      
   }   

 

   __global__ void resetVelocityOpenBndNodes(void) {
   // Set velocity to zero for open boundary nodes where the velocity is imposed.

      NODE_BND(nodebnd);
printf("DDDDDDDDDD\n");  // XXXXXXX
      if (nodebnd->segment != nullptr && nodebnd->segment->ifltype != 0) {
         for (int k = 0; k < nvrt; k++) { nodebnd->node->u2[k] = {0.0, 0.0}; nodebnd->node->w2[k] = 0.0; }  // XXXXXXX
      }

   }


   __global__ void setVelocityOpenBndNodes(void) {
   // For open boundary nodes sets their velocity to the imposed values.

      SIDE_BND(sidebnd);

      if (sidebnd->segment != nullptr && sidebnd->segment->ifltype != 0) {
         // If the sidebnd is an open boundary, set the velocity of the nodes that belong
         // to its segment to the average of the side velocity.

         for (int k = 0; k < nvrt; k++) {
            Node *n1 = sidebnd->side->n1;
            Node *n2 = sidebnd->side->n2;


            Vector2 n1_u2(n1->u2[k]);
            Vector2 n2_u2(n2->u2[k]);

            if (n1->nodeBnd->segment == sidebnd->segment) {
               atomicAdd(&(n1_u2.x), sidebnd->side->u2[k].x/n1->nodeBnd->numSidesForAverage());
               atomicAdd(&(n1_u2.y), sidebnd->side->u2[k].y/n1->nodeBnd->numSidesForAverage());
            }
            if (n2->nodeBnd->segment == sidebnd->segment) {
               atomicAdd(&(n2_u2.x), sidebnd->side->u2[k].x/n2->nodeBnd->numSidesForAverage());
               atomicAdd(&(n2_u2.y), sidebnd->side->u2[k].y/n2->nodeBnd->numSidesForAverage());
            }

         }

      }

   }


   __global__ void computeVerticalVelocity(void) {
      // Computes the vertical velocity at the (2D) element's centroid.
      // In this function we will be considering volumetric elements, with side faces and a 
      // bottom and a top face. We will call them SFs, BF and TF.
      // It computes the flux coming from BF and the SFs.
      // If the bed is moving, it takes that into account.
      // Since water is incompressible, the integral of the velocity times the normal vector
      // along the surface of the prism is zero. In other words, the flux coming from the SFs 
      // and the BF is equal to the flux going out from the TF.
      // TF = - BF - sum(SFs).
      // It then computes the horizontal velocity at the center of TF.
      // The vertical velocity W_TF at the center of TF is computed by knowing that:
      //    (U_TF, V_TF, W_TF) = TF*n_TF/area(TF)
      // where n_TF is the normal vector to TF (3D) and area(TF) is the area of TF (also 3D).

      ELEM(elem);

      for (int l = 0; l < nvrt; l++) elem->w[l] = 0.0; 

      if (elem->isDry) return;

      const int kb = elem->kb;
      const bool isTri = elem->isTri();

      assert(kb >= 0 && "computeVerticalVelocity: wrong kb");

      // Compute the bottom velocity
      if constexpr(imm != 0) {
         assert(false && "computeVerticalVelocity: not implemented yet for imm != 0");
      }


      double bottomArea3D;
      Vector3 elem_bottom_n = elem->n(kb, bottomArea3D); // normal vector to the bottom face of the "volumetric" element (3D)

      assert(elem_bottom_n.z > 0.0 && "computeVerticalVelocity: elem_bottom_n.z <= 0");

      double topArea3D;

      // Bottom face flux
      double bottomFlux = 0.0;
      Vector2 u2meanBottom = {0.0, 0.0};
      const Vector2 gradh = elem->grad2(&Node::depth);
      if constexpr(imm==2) {
         Vector3 ubed;

         assert(false && "computeVerticalVelocity: not implemented for imm==2");
         // updateBDEF(ubed);
     
         // const double bflux0 = dot(ubed, elem_bottom_n);

         // elem->w[kb] += ubed.z;
      }
      else {

         u2meanBottom = elem->mean(kb, &Side::u2Bottom);  // The mean computation includes the projection to the element frame if ics=2
  
         // Computes the vertical velocity such that a horizontal velocity of u2meanBottom has the same slope as the bed
         // in the direction of motion.
         const double bottomW = - dot(gradh, u2meanBottom);  
         // Error: /=0 for 2D (but OK b/cos fluxes are 0 below for transport)  XXXXX

         // VElocity associated to the bed deformation
         double bedVelocity = 0.0;
         if constexpr(imm == 1) assert (false && "computeVerticalVelocity: not implemented for imm==1"); //bedVelocity = (av_bdef2 - av_bdef1)/dt;

         elem->w[kb] = bedVelocity + bottomW; 

         bottomFlux = 0.0;  // No flux from the bed face. XXXX
         
      } // if imm

      // Compute w for kb+1 to nvrt-1
      for (int l = kb; l < nvrt - 1; l++) {

         // normal vector to the top face of the "volumetric" element (3D)
         Vector3 elem_top_n = elem->n(l + 1, topArea3D);  // Get the normal vector and the area of the top face (as a 3D surface)

         assert(topArea3D > 0.0 && "computeVerticalVelocity: topArea3D <= 0");
         assert(elem_top_n.z > 0.0 && "computeVerticalVelocity: elem_top_n.z >= 0");

         // Compute the total flux at the sides
         double sideFlux = 0.0;

         for (int j = 0; j < elem->i34; j++) {
            const Side *side = elem->sides[j];
            const int kbs = side->kb;
            const Vector2 u2_top = side->u2[max(l + 1, kbs)];
            const Vector2 u2_bot = side->u2[max(l,     kbs)];
            Vector2 u2c = 0.5*(u2_top + u2_bot); // velocity at the center of the side

            // The area is "clamped" meaning that it takes kb into account.
            const double sideArea = side->areaClamped(l);

            const double sign = elem->ssign[j]; 

            // sign is positive if the side is oriented in the direction of the normal vector to the element, negative otherwise
            sideFlux += sign*dot(u2c, side->n)*sideArea;
         }
         
         // Compute the horizontal velocity at the center of the top face as an average of the 
         // velocities at the centers of the sides (takes into account the projection to the element frame if ics=2)
         const Vector2 u2meanTop = elem->mean(l + 1, &Side::u2); 

         // The flux at the top ace can be computed from the other flows.
         // All fluxes are in directions into the prism.
         const double topFlux = - sideFlux - bottomFlux;

         // Flux at the top face due to the horizontal velocities only
         // topFlux = topHorizontalFlux + topVerticalFlux
         const double topHorizontalFlux = dot(u2meanTop, toVector2(elem_top_n))*topArea3D; // Part of the flux related to the horizontal velocity

         // From the flux balance, (topVerticalFlux + topHorizontalFlux) + sideFlux + bottomFlux = 0, so:
         const double topVerticalFlux = topFlux - topHorizontalFlux;  

         elem->w[l + 1] = topVerticalFlux/(elem_top_n.z*topArea3D);

         bottomFlux = -topFlux;
         
      }

   }  







// !===============================================================================
// !===============================================================================

} // extern "C"

#ifdef DONTCOMPILE
__device__ void computeVerticalAdvection(void) {
   // Compute the vertical advection term d(uw)/dz and d(vw)/dz, using cubic splines.

   sdbt(1:2,:,:)=0.d0

   const Side *side = FROM kernel;

   if(side->isDry) return;

   // If side is wet...


   // Computes the vertical velocity at the side using the average of the velocities at the element centers
   // XXX WHY DOES IT USE THIS METHOD? IN OTHER PLACES, IT INTERPOLATES FROM W AT THE NODES
   double ws[nvrt] = {} // w @ side
   int icount = 0;
   for (int j=0; j<2; j++) {
      const Element *elem = side->elements[j];
       
      if(elem != nullptr) {
         icount++;
         for (int k=0; k<nvrt; k++) {
            ws[k] = ws[k] + elem->w[k];
         }
      }
   }  // for j

   assert(icount > 0 && "STEP: Wet side has no wet elements");
   ws = ws/icount;

   // Computes the vertical advection term d(uw)/dz, using cubic splines.
   // Splines allow for a smooth transition between vertical levels.
   // Splines have continuity restrictions at the internal knots, but
   // we have to provide the 1st derivatives at the endpoints.
   const int kbs = side->kbs; 
   
   double top_duw_dz = (ws[nvrt - 1]*side->u2[nvrt] - ws[nvrt - 1]*side->u2[nvrt - 1])/(side->z[nvrt - 1] - side->z[nvrt - 2])
   double bot_duw_dz = (ws[kbs     ]*side->u2[kbs ] - ws[kbs     ]*side->u2[kbs  - 1])/(side->z[kbs ] - side->z[kbs  - 1])
   initCubicSpline(kbs, nvrt, side, ws, up)

   double top_dvw_dz = (ws[nvrt - 1]*side->v2[nvrt] - ws[nvrt - 1]*side->v2[nvrt - 1])/(side->z[nvrt - 1] - side->z[nvrt - 2])
   double bot_dvw_dz = (ws[kbs     ]*side->v2[kbs ] - ws[kbs     ]*side->v2[kbs  - 1])/(side->z[kbs ] - side->z[kbs  - 1])
   initCubicSpline(kbs, nvrt, side, ws, vp)

   // Total advection (average onto side)
   const Node *node1 = side->nodes[0];
   const Node *node2 = side->nodes[1];
   for (int k=kbs(i); k<nvrt; k++) {

      Vector2 bcc1 = node1->bcc[k];
      Vector2 bcc2 = node2->bcc[k];
      if constexpr(ics == 2) {
         // Project the velocities from the node frames to the side frame
         bcc1 = node1->frame.projHorVec(bcc1, side->frame);
         bcc2 = node2->frame.projHorVec(bcc2, side->frame);
      }

      sdbt(1,k,i) = side->u2[k] - dt*(bcc1.x + bcc2.x)/2.d0 - dt*side->u2p[k]
      sdbt(2,k,i) = side->v2[k] - dt*(bcc1.y + bcc2.y)/2.d0 - dt*side->v2p[k]

   }
}
#endif


extern "C" {
   __global__ void computeElemAverageGhat(void) {
      // Element average ghat1 (in eframe if ics=2; \breve{G} in notes). Dimensions: m^2/s.

      ELEM(elem);
      
      elem->ghat1= {0.0, 0.0};  

      if constexpr(ihydraulics != 0) {
         if (nhtblocks > 0 && elem->isblock) return;
      }

      
      // ghat1=elem average of \breve{G} includes: Coriolis, atmospheric pressure, tidal potential, horizontal diffusion, and baroclinic etc
      // cff1-3: coefficients for 2D/3D

      elem->gradEta = {0.0, 0.0};  // \nabla \eta
      elem->gradPr  = {0.0, 0.0};  // \nabla pressure
      elem->gradETP = {0.0, 0.0};  // \nabla Earth Tidal Potential

      if (elem->isDry) return;

      // Average of the coefficients
      double av_cff1 = 0.0;  // av. for cff1
      double av_cff2 = 0.0;  // av. for cff2
      double av_cff3 = 0.0;  // av. for cff3
      double av_df = 0.0;  // av. for H^\alpha / H
      double av_cff2_chi = 0.0;  // av. for cff2*chi
      double av_cff3_chi = 0.0;  // av. for cff3*chi

      const int i34 = elem->i34;
      const int kb  = elem->kb;

      for (int j=0; j<i34; j++) {
         // For each node and side in the element...

         const Node *node = elem->nodes[j];
         const Side *side = elem->sides[j];

         const double chi      = side->chi;
         const double veg_c2   = side->veg_c2;
         const double veg_beta = side->veg_beta;
         const double veg_h    = side->veg_h();
         const double veg_c    = veg_c2/(1.0 + veg_c2);  // veg_c2 = \alpha*|u|*dt

         const Vector2 gradPhi = elem->gradPhi[j];

         // Node - const in elem
         elem->gradEta += node->eta2*gradPhi;  // gradPhi is always in the element frame
         elem->gradPr  += node->pr  *gradPhi;

         if (elem->depth >= tip_dp) elem->gradETP += node->ETP*gradPhi;

        
         const double htot  = side->htot2();        // Total depth at the side
         const double zctr2 = side->zTopCanopy();  

         double cff1; 
         double cff2; 
         double cff3;   
         if (side->is2D()) {  // XXX Check indices, it used to have a +1
         // coefficients applied to 2D/3D case 
            cff1 = side->hhat/htot;
            cff2 = 0.0;
            cff3 = 0.0;
         } else { //3D
            //Init for no-veg case
            cff1 = 1.0;
            cff2 = 1.0;
            cff3 = 0.0;

            if constexpr(iveg != 0) {
               // Computes coefficients used with vegetation terms G1, G2
               // cff1, cff3 are used for the first two terms in G1, G2:
               // Emergent vegetation:   G1 = Ustar/(\alpha*|u|*dt)
               // Submergent vegetation: G2 = Ustar - c*Ustar_alpha + ..., with c = \alpha*|u|*dt/(1+\alpha*|u|*dt)
               // In general, cff1*Ustar - cff3*Ustar_alpha
               if(veg_h < 0.99*htot) { // submergent vegetation
                  cff1 = 1.0;
                  cff2 = 1.0 + veg_beta*veg_c;          // Term (1 + beta*c) in G2.
                  cff3 = veg_c;

               } else { // emergent vegetation
                  cff1 = 1.0/(1.0 + veg_c2);
                  cff2 = cff1;
                  cff3 = 0.0;

               }
            } // if iveg
         } //  2D/3D

         av_cff1 += cff1/i34;
         av_cff2 += cff2/i34;
         av_cff3 += cff3/i34;
         av_cff2_chi += cff2*chi/i34;
         av_cff3_chi += cff3*chi/i34;
         av_df += min(htot, veg_h)/(htot*i34);


         // btrack values: U^\star, U^{\star\alpha}
         // Computes the backtracked velocity integrated vertically (full depth and to the top of the canopy)
         Vector2 Ustar       = {0.0, 0.0}; 
         Vector2 Ustar_alpha = {0.0, 0.0}; 

         if (side->is2D()) {  //2D; no need to compute U^{\star\alpha}
               Ustar += side->btu[nvrt - 1]*side->Dz(nvrt - 1);
         } else { //3D
            for (int k = side->kb + 1; k < nvrt; k++) {
            
               Vector2 w2 = (side->btu[k] + side->btu[k-1])/2.0*side->Dz(k);
               Ustar += w2;

               if constexpr(iveg != 0) {
                   if (zctr2 > side->zc[k-1]) {  // XXX CHECK INDICES!!! 
                       const double zrat = clamp((zctr2 - side->zc[k-1])/side->Dz(k), 0.0, 1.0); // For the top of the canopy, ratio of the level that is below the canopy
                       Ustar_alpha += zrat*w2; // XXXX
                   }
               }  // if iveg
            }  // for k
         }  //2D/3D

         // If needed, project the side velocities to the element frame
         if constexpr(ics==2) {
            Ustar       = side->frame.projHorVec(Ustar,       elem->frame);
            Ustar_alpha = side->frame.projHorVec(Ustar_alpha, elem->frame);
         }
        

         // Add btrack & surface stress terms to \breve{G} first
         elem->ghat1 += cff1*Ustar - cff3*Ustar_alpha;


         // Both tau and ubstar are in the side frame (for both ics==1 and ics==2)
         Vector2 tau    = side->tau();  
         Vector2 ubstar = side->btu[side->kb + 1];  

         // If needed, project tau and ubstar to the element frame
         if constexpr(ics==2) {
            tau    = side->frame.projHorVec(tau,    elem->frame);  // Here we are doing two projections: from the node frame to the side frame, and from the side frame to the element frame. This is different from the original code where they only did one projection.
            ubstar = side->frame.projHorVec(ubstar, elem->frame);
         }

         // Term for surface and bottom stress
         elem->ghat1 -= chi*dt*cff2*ubstar + dt*cff1*tau;


         // All terms in F,F^\alpha,f_b except baroclinic
         // F, f_b: all explicit terms excluding vegetation 
         // F^\alpha=\int_{-h}^z_v f dz (f ncludes all explicit terms except for vegetation)
         // Only need to calculate F^\alpha for 3D sides (but init=0) 
         const double coriolis = side->coriolis;
         Vector2 F  =  coriolis*side->U;  // F_x
         Vector2 fb =  coriolis*side->u2[side->kb + 1];  // f_b_x
         Vector2 Fa = {0.0, 0.0};  // F^\alpha, vertically forcing terms baroclinic gradient, horizontal viscosity, Coriolis, earth tidal potential, atmospheric pressure, radiation stress)

         if (iveg != 0 && nvrt > side->kb) {  // 3D
            
            for (int k = side->kb; k < nvrt; k++) {
               if (zctr2 > side->zc[k-1]) {
                   const Vector2 w2 = side->Dz(k)*(side->u2[k] + side->u2[k-1])/2.0;
                   const double zrat = clamp((zctr2 - side->zc[k-1])/side->Dz(k), 0.0, 1.0);
                   Fa +=  coriolis*w2*zrat;
               } // if zctr2
            } // k
         } // iveg

         // Integrate horizontal viscosity force vertically
         for (int k=side->kb + 1; k<nvrt; k++) {
            const Vector2 levelHorVisc = side->Dz(k)*(side->nuGradU[k] + side->nuGradU[k-1])/2.0; 

            F += levelHorVisc;

            // Vegetation term
            if constexpr(iveg!=0) {
               if (nvrt > side->kb + 1 && zctr2 > side->zc[k-1]) { // 3D
                  const double zrat = clamp((zctr2 - side->zc[k-1])/side->Dz(k), 0.0, 1.0);
                  Fa += levelHorVisc*zrat;
               } 
            } // if iveg
         } // for k

         // Get nu*grad(U) at the side's bottom level
         Vector2 nuGradU = side->nuGradU[side->kb];

         fb += nuGradU;

         // If needed, project vectorsF, fb and Fa to the element frame
         if constexpr(ics==2) {
            F    = side->frame.projHorVec(F,    elem->frame);
            fb   = side->frame.projHorVec(fb,   elem->frame);
            Fa   = side->frame.projHorVec(Fa,   elem->frame);
         } // if ics

         elem->ghat1 += cff1*dt*F - cff2*chi*dt*dt*fb - cff3*dt*Fa;
      } // for j (nodes and sides)

      elem->ghat1 = (1.0/i34)*elem->ghat1;


      // Finish off terms in F, F^\alpha and f_b
      const double grav3 = elem->mean(&Node::grav2);
      const Vector2 fb = grav3*elem->gradETP - elem->gradPr/rho0; //const in each elem
      Vector2 tmp = {0.0, 0.0}; //elem average of all terms; into ghat1
      for (int j=0; j<i34; j++) { 
         // For each side in the element...
         const Side *side = elem->sides[j];
         const double chi = side->chi;

         const double htot = side->htot2();
         const double veg_h = side->veg_h();
         const Vector2 F  = htot*fb;
         const Vector2 Fa = min(htot, veg_h)*fb;

         tmp += dt*(av_cff1*F - av_cff2*chi*dt*fb - av_cff3*Fa);
      } // for j

      // tmp is already in elem frame
      elem->ghat1 += tmp/i34;

      // Baroclinic force
      if constexpr(ibc==0) {
         assert(elem->nodes[0]->rho[0] >= -98.0 && elem->nodes[1]->rho[0] >= -98.0 && elem->nodes[2]->rho[0] >= -98.0 && "computeElemAverageGhat1: Impossible dry");

         Vector2 Fa = {0.0, 0.0};
         Vector2 intGradRho = {0.0, 0.0};
         
         // Computes \sum_{l=k}^N dr*dz for the k level at trhe element frame
         for (int k=nvrt-2; k >= kb + 1; k--) { /// XXXX CHECK INDICES
            intGradRho += elem->gradRho[k+1]*elem->Dz(k+1); 

            // \in F [m^2/s/s]
            Fa -= grav3/rho0*elem->Dz(k)*(2.0*intGradRho + elem->gradRho[k]*elem->Dz(k));
         }

         const Vector2 fb = -grav3/rho0*intGradRho; //[m/s/s]; elem average
         const Vector2 F = av_df*Fa; //approx
         elem->ghat1 += rampbc*dt*(av_cff1*F - dt*av_cff2_chi*fb - av_cff3*Fa);

      } // if ibc



   }

} // extern "C"



#ifdef DONTCOMPILE
__global__ void checkMassConservation(Tracers* sumtr) {

   int i = blockIdx.x * blockDim.x + threadIdx.x;

   if (i >= mesh->ne) return;  // Guard against out-of-bounds access

   Element *elem = &(mesh->localElemData[i]);

   // Skip dry elements
   if (elem->isDry) return;

   // Local sum array for this thread
   double localSum[N_TRACERS] = {0.0};

   // Sum over vertical levels for this element
   const int kbe = elem->kbe;
   const double area = elem->area;

   for (int k = kbe; k < nvrt; k++) {
      const double vol = (elem->z[k] - elem->z[k-1]) * area;
      
      // Accumulate tracer mass
      for (int tr = 0; tr < N_TRACERS; tr++) {
         sumtr[i].tracer[tr] += vol * elem->tracers[k].tracer[tr];
      }
   }

   //  This function returns per element sums. A reduction must be done on the host side.
}


__global__ void transport(void) {

   int ie = blockIdx.x * blockDim.x + threadIdx.x;

   if (ie >= mesh->ne) return;  // Guard against out-of-bounds access

   Element *elem = &(mesh->localElemData[ie]);

   if (ip_weno==2 && nweno2(ie)>0 && isten_qual2(ie) && isbe(ie)==0)  // p2 weno method

               for (int m=0; m<ntr; m++) {
                  for (int k=kbe(ie); k<nvrt; k++) {
                     // Calculate each polynomial's weight 
                     wm=0.d0; sum1=0.d0; 
                     wm1=0.d0; wm2=0.d0; 
                     for (int i=0; i<elem->nweno2; i++) {
                  
                     const double b1 = dot_product(elem->wts2(:,1,i),elem->tr_el(m,k,elem->isten2(1:6,i)))
                     const double b2 = dot_product(elem->wts2(:,2,i),elem->tr_el(m,k,elem->isten2(1:6,i)))
                     const double b3 = dot_product(elem->wts2(:,3,i),elem->tr_el(m,k,elem->isten2(1:6,i)))
                     const double b4 = dot_product(elem->wts2(:,4,i),elem->tr_el(m,k,elem->isten2(1:6,i)))
                     const double b5 = dot_product(elem->wts2(:,5,i),elem->tr_el(m,k,elem->isten2(1:6,i)))
                     b1 = b1*b1+b2*b2+(4.d0*b1*b3+2*b2*b4)*fwts2(1,ie)+(4.d0*b2*b5+2*b1*b4)*fwts2(2,ie)
                     b2 = (4.d0*b3*b3+b4*b4)*fwts2(3,ie)+4.d0*b4*(b3+b5)*fwts2(4,ie)+(4.d0*b5*b5+b4*b4)*fwts2(5,ie)
                     b3 = (4.d0*b3*b3+b4*b4+4.d0*b5*b5)*area(ie)
                     // -----------------Hu and Shu (1999)'s method----------------
                     b4=b1+b2+b3

//Error: add PRODUCTION CPP
#ifdef DEBUG
                  if(b4<0.d0 && abs(b4)>1.0d-50) {
                    write(errmsg,*)'b4<0',b1,b2,b3,b4
                    call parallel_abort(errmsg)
                  }
#endif

                  wm(i) = 1.0d0/((elem->epsilon2+b4)*(elem->epsilon2+b4))
                  // !----test-------------- 
                  // !wm1(i)=1.0/((epsilon2_elem(ie)+b1+b2)*(epsilon2_elem(ie)+b1+b2)); ! sum1=sum1+wm1(i)
                  // !wm2(i)=1.0/((epsilon3+b3)*(epsilon3+b3)); 
                  // !wm(i)=wm1(i)*wm2(i)
                  // !---------------------- 
                  sum1=sum1+wm(i)
                } // end for i

                wm=wm/sum1

                for (int j=0; j<elem->i34; j++) {

                  const Side *side  = elem->sides[j];
                  const Node *node1 = side->n1;
                  const Node *node2 = side->n2;

                  if (elem->ssign[j]*side->flux_adv_hface[k] >= 0.d0) { // outflow face or bnd face
                     trsd=0.d0
                     for (int kk=0; kk<nquad; kk++) {
                        for (int i=0; i<elem->nweno2; i++) {
                           trsd(kk)=trsd(kk)+wm(i)*dot_product(elem->wmat2(1:6,i,kk,j),elem->tr_el(m,k,elem->isten2(1:6,i)))
                        } // end for i   
                     } // end for kk
                     trsd_tmp(m,k,jsj)=sum(trsd(1:nquad))/real(nquad,rkind) //mean value of the two quadrature points
                  } // outward flux
               } // end for j=1,i34(ie) ! side
            }
         }
              

//             elseif((ip_weno==1.or.ip_weno==2).and.nweno1(ie)>0.and.isbe(ie)==0) then !p1 weno method
//             //   calculate each polynomial's weight 
//               do m=1,ntr; do k=kbe(ie)+1,nvrt!!!possible optimization
//                 wm=0.d0; sum1=0.d0
//                 do i=1,nweno1(ie) !for each polynomial
//                   b1=dot_product(wts1(:,1,i,ie),tr_el(m,k,isten1(1:3,i,ie)))
//                   b2=dot_product(wts1(:,2,i,ie),tr_el(m,k,isten1(1:3,i,ie)))
//                   b1=b1*b1+b2*b2 !smoothness indicator
//                   b2=1.d0/((epsilon1+b1)*(epsilon1+b1))
//                   wm(i)=b2
//                   sum1=sum1+b2
//                   !wm(i)=1.d0
//                   !sum1=sum1+1.d0

// // Error: PRODUCTION
// #ifdef DEBUG
//                   if(.not.(wm(i)>=0.d0.or.wm(i)<0.d0)) then
//                     write(errmsg,*)'wm(i)=nan',wts1(:,:,i,ie),tr_el(m,k,isten1(1:3,i,ie))
//                     call parallel_abort(errmsg)
//                   endif
// #endif

//                 enddo !i
//                 wm=wm/sum1 !final weight

                for (int j=0; j<elem->i34; j++) { // side
                  const Side *side = elem->sides[j];
                  const Node *node1 = side->n1;
                  const Node *node2 = side->n2;

                  if ((elem->ssign[j]*side->flux_adv_hface[k]>=0.d0) .or. side->isbs != 0) {  outflow face or bnd face
                    trsd=0.d0
                    for (int kk=0; kk<nquad; kk++) { //1 or 2 quadrature points
                      for (int i=0; i<elem->nweno1; i++) {
                        trsd(kk)=trsd(kk)+wm(i)*dot_product(elem->wmat1(1:3,i,kk,j),elem->tr_el(m,k,elem->isten1(1:3,i)))
                      }
                    enddo //kk
                    trsd_tmp(m,k,jsj)=sum(trsd(1:nquad))/real(nquad,rkind) //mean value of the two quadrature points
                  } // outward flux
                } // end for j=1,i34(ie) ! side
              } // end for k; ntr 

            else !simple upwind method
              do j=1,i34(ie)
                jsj=elside(j,ie)
                n1=isidenode(1,jsj);   n2=isidenode(2,jsj)
                do k=kbe(ie)+1,nvrt !!!possible optimization
                  if(ssign(j,ie)*flux_adv_hface(k,jsj)>=0.d0) then !outflow face
                    trsd_tmp(:,k,jsj)=tr_el(:,k,ie) 
                  endif
                enddo  !vertical layers
              enddo !j=1,i34(ie) ! side
            endif !order of accuracy

         //  enddo !ie
// !$OMP     end do

// !$OMP     master
//           !message passing
//           if (nproc>1) then
//             !make a copy of side concentrations
//             trsd_tmp0=trsd_tmp 

//             !this exchange routine swaps interface concentration between two partitions
//             !trsd_tmp0 receives neighboring info, while trsd_tmp remains the same
//             call exchange_s3d_tr3(trsd_tmp0,trsd_tmp) 

//             !Two sets of concentrations: before message-passing (trsd_tmp) and after (trsd_tmp0).
//             !(1) trsd_tmp0 is the same as trsd_tmp everywhere except on interface sides;
//             !(2) For an interface side, non-0 concentration can not be present in both trsd_tmp and trsd_tmp0;


//             !This obtains final values for interface sides
//             trsd_tmp(:,:,iside_table)= trsd_tmp(:,:,iside_table)+trsd_tmp0(:,:,iside_table)
//       :    endif !nproc>1
//           !message passing end
// !$OMP     end master
// !$OMP     barrier


// !         Do advection; conc @ dry elem will not be changed
// !         Follows the style of explicit TVD, but using WENO correction instead of TVD correction
// !$OMP     do reduction(+: h_mass_in)
//           do i=1,ne
//             if(idry_e(i)==1) cycle


// !           Wet elements with 3 wet nodes
//             do k=kbe(i)+1,nvrt


//               bigv=area(i)*(ze(k,i)-ze(k-1,i)) !volume
//               dtb_by_bigv = dtb/bigv
    
// !             Advective flux
//               if (ntd_weno==1 .or. itd_weno==1) then
//                 adv_tr(1:ntr)=trel_tmp(1:ntr,k,i) 
//               else
//                 adv_tr(1:ntr)=(1.d0-rk_coef(itd_weno))*trel_tmp0(1:ntr,k,i)+rk_coef(itd_weno)*trel_tmp(1:ntr,k,i) 
//               endif

//               !psum=0.d0 !sum of modified fluxes for all inflow bnds
//               do j=1,i34(i)
//                 jsj=elside(j,i) !side
//                 iweno=.true. !use weno reconstruction as default
//                 iel=ic3(j,i) !neighboring elem

//                 if(iel/=0) then
//                   if(idry_e(iel)==1) cycle
//                   trel_tmp_outside(:)=trel_tmp(:,k,iel)
//                   !reset to upwind; if Element i uses ELM, it will be overwritten by ELM value eventually
//                   if(iupwind_e(i)+iupwind_e(iel)>0) then
//                     iweno=.false.
//                   endif
//                 else !land or open bnd side
//                   !skip land bnd and outflowing open bnd, 
//                   !in this case iweno is not changed from initial value (.true.), but won't be used
//                   if(isbs(jsj)<=0.or.k>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(k,jsj)>=0.d0) then
//                     !For outflow open bnd side, estimate mass in (open
//                     !side cannot be interface side)
//                     if(max_iadjust_mass_consv>0.and.isbs(jsj)>0) then !outflow @open bnd
//                       h_mass_in(:)=h_mass_in(:)-trel_tmp(1:ntr,k,i)*flux_adv_hface(k,jsj)*dtb
//                     endif

//                     cycle
//                   endif

//                   !Open bnd side with inflow (must be wet elem) or k<kbs(jsj)+1
//                   iweno=.false.  !reset to upwind
//                   ibnd=isbs(jsj) !global bnd #
//                   !Find node indices on bnd segment for the 2 nodes (for type 4 b.c.)
//                   nwild(1:2)=0
//                   do ll=1,2 !nodes
//                     ndo=isidenode(ll,jsj)
//                     do lll=1,2 !2 possible bnds
//                       if(isbnd(lll,ndo)==ibnd) then
//                         nwild(ll)=isbnd(-lll,ndo) !global index
//                         exit
//                       endif
//                     enddo !lll
//                   enddo !ll
//                   ind1=nwild(1); ind2=nwild(2);

//                   do jj=1,natrm
//                     if(ntrs(jj)<=0) cycle

//                     if(itrtype(jj,ibnd)==0) then !set to be same as interior (so cancel out below)
//                       do ll=irange_tr(1,jj),irange_tr(2,jj)
//                         trel_tmp_outside(ll)=trel_tmp(ll,k,i)
//                       enddo !ll
//                     else if(itrtype(jj,ibnd)==1.or.itrtype(jj,ibnd)==2) then
//                       do ll=irange_tr(1,jj),irange_tr(2,jj)
//                         trel_tmp_outside(ll)=trobc(jj,ibnd)*trth(ll,1,1,ibnd)+(1-trobc(jj,ibnd))*trel_tmp(ll,k,i)
//                       enddo !ll
//                     else if(itrtype(jj,ibnd)==3) then
//                       do ll=irange_tr(1,jj),irange_tr(2,jj)
//                         tmp=sum(tr_nd0(ll,k,elnode(1:i34(i),i))+tr_nd0(ll,k-1,elnode(1:i34(i),i)))/2.d0/real(i34(i),rkind)
//                         trel_tmp_outside(ll)=trobc(jj,ibnd)*tmp+(1.d0-trobc(jj,ibnd))*trel_tmp(ll,k,i)
//                       enddo !ll
//                     else if(itrtype(jj,ibnd)==4) then
//                       do ll=irange_tr(1,jj),irange_tr(2,jj)
//                         trel_tmp_outside(ll)=trobc(jj,ibnd)* &
//        &(trth(ll,k,ind1,ibnd)+trth(ll,k,ind2,ibnd)+trth(ll,k-1,ind1,ibnd)+trth(ll,k-1,ind2,ibnd))/4.d0+ &
//        &(1.d0-trobc(jj,ibnd))*trel_tmp(ll,k,i)
//                       enddo !ll
//                     else
//                       write(errmsg,*)'TRASNPORT: INVALID VALUE FOR ITRTYPE:',jj,ibnd
//   !'
//                       call parallel_abort(errmsg)
//                     endif !itrtype

//                   enddo !jj

//                   !Tally horizontal mass in for inflow case
//                   if(max_iadjust_mass_consv>0) h_mass_in(:)=h_mass_in(:)-trel_tmp_outside(1:ntr)*flux_adv_hface(k,jsj)*dtb
//                 endif !iel

//                 !!upwind part
//                 !if(k>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(k,jsj)<0) then !inflow
//                 !  do jj=1,ntr
//                 !    adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp_outside(jj)-trel_tmp(jj,k,i))
//                 !  enddo
//                 !endif
//                 !!weno correction
//                 !if (iweno) then
//                 !  if(k>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(k,jsj)<0) then !inflow
//                 !    do jj=1,ntr
//                 !      adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trsd_tmp(jj,k,jsj)-trel_tmp_outside(jj))
//                 !    enddo
//                 !  else if(k>=kbs(jsj)+1.and.ssign(j,i)*flux_adv_hface(k,jsj)>=0) then !outflow
//                 !    do jj=1,ntr
//                 !      adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trsd_tmp(jj,k,jsj) )
//                 !    enddo
//                 !  endif !inflow or outflow
//                 !endif !if use weno

//                 !alternative formulation
//                 if (k>=kbs(jsj)+1) then
//                   if (iweno) then !weno
//                     do jj=1,ntr
//                       tmp=rk_coef(itd_weno)*dtb_by_bigv*( ssign(j,i)*(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trsd_tmp(jj,k,jsj)) )
//                       adv_tr(jj)=adv_tr(jj)+tmp
//                       if (adv_tr(jj)<tr_min_max(1,jj) .or. adv_tr(jj)>tr_min_max(2,jj) ) then !reset to upwind
//                         adv_tr(jj)=adv_tr(jj)-tmp !reset to previous value
//                         if(ssign(j,i)*flux_adv_hface(k,jsj)<0.d0) then  !inflow 
//                           adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp_outside(jj)-trel_tmp(jj,k,i))
//                         endif
//                       endif
//                       !when ntd_weno=1, reduces to: adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*( ssign(j,i)*(flux_adv_hface(k,jsj))*(trel_tmp(jj,k,i)-trsd_tmp(jj,k,jsj)) )
//                     enddo
//                   else !upwind
//                     if(ssign(j,i)*flux_adv_hface(k,jsj)<0) then  !inflow 
//                       do jj=1,ntr
//                         adv_tr(jj)=adv_tr(jj)+dtb_by_bigv*abs(flux_adv_hface(k,jsj))*(trel_tmp_outside(jj)-trel_tmp(jj,k,i))
//                       enddo
//                     endif
//                   endif !weno or upwind
//                 endif

//               enddo !j

//               tr_el(1:ntr,k,i)=adv_tr(1:ntr) 
//             enddo !k=kbe(i)+1,nvrt

//             !Overwrite with ELM value if enabled.
//             if(ielem_elm(i)/=0) then
//               rat=time_r/dt !time ratio
//               rat=max(0.d0,min(1.d0,rat))
//               do k=kbe(i)+1,nvrt
//                 do jj=1,ntr
//                   psumtr(jj)=0.d0
//                   do j=1,i34(i)
//                     jsj=elside(j,i)
//                     n1=isidenode(1,jsj); n2=isidenode(2,jsj)
//                     swild4(1,1)=0.5d0*(tr_nd(jj,k,n1)+tr_nd(jj,k,n2))
//                     swild4(2,1)=0.5d0*(tr_nd(jj,k-1,n1)+tr_nd(jj,k-1,n2))
  
//                     swild4(1,1)=swild4(1,1)*rat+sdbt(2+jj,max(k,kbs(jsj)),jsj)*(1-rat)
//                     swild4(2,1)=swild4(2,1)*rat+sdbt(2+jj,max(k-1,kbs(jsj)),jsj)*(1-rat)
//                     psumtr(jj)=psumtr(jj)+0.5d0*(swild4(1,1)+swild4(2,1)) 
//                   enddo !j
//                 enddo !jj
//                 tr_el(1:ntr,k,i)=psumtr(:)/dble(i34(i))
//               enddo !k
//             endif !ielem_elm(i)/=0


// !           Extend
//             do k=1,kbe(i)
//               tr_el(1:ntr,k,i)=tr_el(1:ntr,kbe(i)+1,i)
//             enddo !k
//           enddo !i=1,ne
// !$OMP     end do 

// !$OMP     master
// !         Update ghosts
// #ifdef INCLUDE_TIMING
//           cwtmp=mpi_wtime()
//           timer_ns(1)=timer_ns(1)+cwtmp-cwtmp2
// #endif
//           call exchange_e3d_2t_tr(tr_el)

// #ifdef INCLUDE_TIMING
//           cwtmp2=mpi_wtime()
//           wtimer(9,2)=wtimer(9,2)+cwtmp2-cwtmp
// #endif

#endif


extern "C" 

__global__ void createMesh(int partNum, 
                           int npgb, int negb, int nsgb, int nseg,
                           int np, int ne, int ns, 
                           int max_np, int max_ne, int max_ns,
                           NodeBoundary *plandNodeBndMem, BoundarySegment *plandSegMem,
                           Node **nodesPtr, Element **elemsPtr, Side **sidesPtr,
                           Node *nodes, Element *elems, Side *sides,
                           Node **partNodeNeighNodesMem, Element **partNodeNeighElemsMem,
                           int16_t *localNodePermutationSortMem,
                           double *partCoeffMatrixMem,
                           NodeBoundary *localNodeBoundaryMem, SideBoundary *localSideBoundaryMem, BoundarySegment *localSegmentMem,
                           int *kb, double h_c, double h_s, double theta_b, double theta_f, double *sigma, 
                           double *ztot, double *cs, double *dcs,
                           double *_debugPlot, int *_debugPlotIdx) {



   if (mesh != nullptr) 
   {
      printf("ERROR: Mesh is already initialized\n");
      return;
   };


   debugPlot = _debugPlot;
   debugPlotIdx = _debugPlotIdx;
   // return;   

   printf("partNum = %d, npgb = %d, negb = %d, nsgb = %d, nseg = %d\n", partNum, npgb, negb, nsgb, nseg);
   printf("np = %d, ne = %d, ns = %d\n", np, ne, ns);
   printf("max_np = %d, max_ne = %d, max_ns = %d\n", max_np, max_ne, max_ns);
   printf("nvrt = %d, nzlev = %d\n", nvrt, nzlev);

   printf("Creating mesh for partition %d, np = %d\n", partNum, np);
   printf("CREATEMESH: elemsPtr = %p, elems = %p\n", *elemsPtr, elems);
   
   mesh = new Mesh(partNum, npgb, negb, nsgb, 
                   np, ne, ns,
                   max_np, max_ne, max_ns, nseg, 
                   plandNodeBndMem, plandSegMem,
                   nodesPtr, elemsPtr, sidesPtr, 
                   nodes, elems, sides,
                   partNodeNeighNodesMem, partNodeNeighElemsMem, 
                   localNodePermutationSortMem,
                   partCoeffMatrixMem,
                   localNodeBoundaryMem, localSideBoundaryMem, localSegmentMem,
                   kb, h_c, h_s, theta_b, theta_f, sigma, ztot, cs, dcs);

   // XXX read from file
   polyDischarge = new Polynomial<4>(0.0, 0.0, 0.0, 0.0, 0.0);                   
   
}

extern "C"
                                 
  __global__ void populateMesh(int *i34, int *elnode, int *elside,
                   double *node_coords, double *node_coords_ll, bool use_ll,
                   double *wind_rotate_angle,  int *sharedIdxNodeBoundary, int *sharedIdxSideBoundary,
                   int *isdel, int *isidenode, int *isbnd,
                   int *nond, int *nlnd, int *iond, int *ilnd, 
                   int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                   int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                   int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                   int *localNodeIndices, int *localElemIndices, int *localSideIndices,
                   int *localNodePart,    int *localElemPart,    int *localSidePart,
                   double *drag, double *roughness,
                   double *emo, double *efa, 
                   double *umo, double *ufa, 
                   double *vmo, double *vfa) {

    mesh->populateMesh(i34, elnode, elside, node_coords, node_coords_ll, use_ll, wind_rotate_angle, 
                       sharedIdxNodeBoundary, sharedIdxSideBoundary,
                       isdel, isidenode, isbnd,
                       nond, nlnd, iond, ilnd, 
                       nodeNeighNodesStart, nodeNeighNodesIndices,
                       nodeNeighElemsStart, nodeNeighElemsIndices,
                       globalNodeIndices, globalElemIndices, globalSideIndices,
                       localNodeIndices, localElemIndices, localSideIndices,
                       localNodePart, localElemPart, localSidePart,
                       drag, roughness,
                       emo, efa, umo, ufa, vmo, vfa);  // XXXXX Remember to refactor this, removing the parameters that are not used.
  }

extern "C" __global__ void populateNodes(int *i34, int *isbnd, int *sharedIdxNodeBoundary,
   double *node_coords, double *node_coords_ll, bool use_ll, double *wind_rotate_angle, 
   int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
   int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
   int16_t *nodePermutationSort,
   int *localNodeNeighNodesStart, int *localNodeNeighElemsStart,
   int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
   int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
   int *localNodePart,     int *localElemPart,     int *localSidePart,
   double *drag, double *roughness, double *rmanning, int *kb, double *localSigma) {

   mesh->populateNodes(i34, isbnd, sharedIdxNodeBoundary, 
                       node_coords, node_coords_ll, use_ll, wind_rotate_angle, 
                       nodeNeighNodesStart, nodeNeighNodesIndices,
                       nodeNeighElemsStart, nodeNeighElemsIndices,
                       nodePermutationSort,
                       localNodeNeighNodesStart, localNodeNeighElemsStart,
                       globalNodeIndices, globalElemIndices, globalSideIndices,
                       localNodeIndices, localElemIndices, localSideIndices,
                       localNodePart, localElemPart, localSidePart,
                       drag, roughness, rmanning, kb, localSigma);
}  

extern "C" __global__ void populateNodes2(void) {

   mesh->populateNodes2();
}


extern "C" __global__ void populateSides(int *i34, int *sharedIdxSideBoundary,
                                         double *node_coords, double *node_coords_ll, bool use_ll,
                                         double *wind_rotate_angle,   
                                         int *isdel, int *isidenode, int *isbnd,
                                         int *nond, int *nlnd, int *iond, int *ilnd,                                  
                                         int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                                         int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                                         int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                                         int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                                         int *localNodePart,     int *localElemPart,     int *localSidePart,
                                         double *drag, double *roughness) {
   mesh->populateSides(i34, sharedIdxSideBoundary,
                       node_coords, node_coords_ll, use_ll,
                       wind_rotate_angle,   
                       isdel, isidenode, isbnd,
                       nond, nlnd, iond, ilnd,                                  
                       nodeNeighNodesStart, nodeNeighNodesIndices, nodeNeighElemsStart, nodeNeighElemsIndices, globalNodeIndices, globalElemIndices, globalSideIndices, localNodeIndices, localElemIndices, localSideIndices, localNodePart, localElemPart, localSidePart, drag, roughness);                          

}

extern "C" __global__ void populateElements(int *i34, int *elnode, int *elside,
   double *node_coords, double *node_coords_ll, bool use_ll,
   double *wind_rotate_angle,   
   int *isdel, int *isidenode, int *isbnd,
   int *nond, int *nlnd, int *iond, int *ilnd,                                  
   int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
   int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
   int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
   int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
   int *localNodePart,     int *localElemPart,     int *localSidePart,
   double *drag, double *roughness) {
   mesh->populateElements(i34, elnode, elside,
                          node_coords, node_coords_ll, use_ll,
                          wind_rotate_angle,   
                          isdel, isidenode, isbnd,
                          nond, nlnd, iond, ilnd,                                  
                          nodeNeighNodesStart, nodeNeighNodesIndices, nodeNeighElemsStart, nodeNeighElemsIndices, globalNodeIndices, globalElemIndices, globalSideIndices, localNodeIndices, localElemIndices, localSideIndices, localNodePart, localElemPart, localSidePart, drag, roughness);                          

} 

extern "C" __global__ void populateBoundaryNodes(int *nond, int *iond,
                                                double *emo, double *efa, 
                                                double *umo, double *ufa, 
                                                double *vmo, double *vfa,
                                                int *ifltype, int *iettype) {
   mesh->populateBoundaryNodes(nond, iond, emo, efa, umo, ufa, vmo, vfa, ifltype, iettype);

}


extern "C" __global__ void populateBoundarySideLinks(void) {
   // mesh->populateBoundarySideLinks();   // XXXXXX Not really used at the moment

}





extern "C"  __global__ void destroyMeshCPU(void) {
   delete mesh;
}




__global__ void getStructSizes(int* sizes)  {
   // Get the size of each struct. This is used to allocate memory for the structs in the host code.
      printf(">>>>>*** C++ 3 getStructSizes, sizes = %p\n", sizes);
      sizes[0] = (int)sizeof(Node);
      sizes[1] = (int)sizeof(Element);
      sizes[2] = (int)sizeof(Side);
      sizes[3] = (int)sizeof(NodeBoundary);
      sizes[4] = (int)sizeof(SideBoundary);
      sizes[5] = (int)sizeof(BoundarySegment);
}


__global__ void printMatrix(bool full) {
   // Debugging: print the matrix

   const int size = mesh->np;
   double sum = 0.0;
   for (int idx_node = 0; idx_node < size; idx_node++) {
      const Node *node = mesh->localNodeData + idx_node;
      if (full) printf("idx_node = %d [%d] diag=%lf, [", idx_node, node->numNodes, node->coeffDiag);
      for (int j = 0; j < node->numNodes; j++) {
         Node* nbNode = node->nodes[j];
         if (full && node->coeffMatrix[j] != 0.0) printf("%d:%f,  ", nbNode->idx, node->coeffMatrix[j]);
         sum += node->coeffMatrix[j];
      }
      if (full) printf("]\n");


   }
   if (!full) printf("sum=%lf\n", sum);
}




__global__ void returnMatrix(double *res) {
   // Debugging: returns the matrix as a CuPy array

   const int size = mesh->np;
   for (int idx_node = 0; idx_node < size; idx_node++) {
      const Node *node = mesh->localNodeData + idx_node;
      res[idx_node*size + idx_node] = *(node->coeffDiag);
      for (int j = 0; j < node->numNodes; j++) {
         Node* nbNode = node->nodes[j];
         if (node->coeffRow[j] != 0.0) res[nbNode->idx*size + idx_node] =  node->coeffRow[csrSlot[j]];
      }

   }
   
}

__global__ void returnX(double *res) {
   // Debugging: returns the matrix as a CuPy array

   const int size = mesh->np;
   for (int idx_node = 0; idx_node < size; idx_node++) {
      const Node *node = mesh->localNodeData + idx_node;

      res[idx_node] = node->mX;

   }
   
}


__global__ void returnB(double *res) {
   // Debugging: returns the matrix as a CuPy array

   const int size = mesh->np;
   for (int idx_node = 0; idx_node < size; idx_node++) {
      const Node *node = mesh->localNodeData + idx_node;

      res[idx_node] = node->qel;

   }
   
}



extern "C" __global__ void getElementVarFunc(int varCode, double *res) {

   ELEM(elem);

   switch (varCode) {
      case (int)VarCode::outputCoords:
         const Point2 coords = elem->outputCoords();
         res[elem->idx*2    ] = coords.x;
         res[elem->idx*2 + 1] = coords.y;
         break;
   }
}


extern "C" __global__ void getSideVarFunc(int varCode, double *res) {

   SIDE(side);

   switch (varCode) {
      case (int)VarCode::outputCoords:
         const Point2 coords = side->outputCoords();
         res[side->idx*2    ] = coords.x;
         res[side->idx*2 + 1] = coords.y;
         break;
   }
}


extern "C" __global__ void getNodeVarFunc(int varCode, double *res) {

   NODE(node);

   switch (varCode) {
      case (int)VarCode::outputCoords:
         const Point2 coords = node->outputCoords();
         res[node->idx*2    ] = coords.x;
         res[node->idx*2 + 1] = coords.y;
         break;
   }
}


extern "C" __global__ void getElementVarAddr(int varCode, void **addr, int *offset, int *sizeElement, int *sizeType) {
   Element *elems = mesh->localElemData;

   *addr = elems;
   *sizeElement = sizeof(Element);
   switch (varCode) {
      case (int)VarCode::z:
         *offset = (int)((char*)&elems[0].z - (char*)elems);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::x:
         *offset = (int)((char*)&elems[0].center.x - (char*)elems);
         *sizeType = sizeof(Vector2);
         break;
      case (int)VarCode::gradETP:
         *offset = (int)((char*)&elems[0].gradETP.x - (char*)elems);
         *sizeType = sizeof(Vector2);
         break;
      case (int)VarCode::area:
         *offset = (int)((char*)&elems[0].area - (char*)elems);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::depth:
         *offset = (int)((char*)&elems[0].depth - (char*)elems);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::radius:
         *offset = (int)((char*)&elems[0].radius - (char*)elems);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::ghat:
         *offset = (int)((char*)&elems[0].ghat1.x - (char*)elems);
         *sizeType = sizeof(Vector2);
         break;
      case (int)VarCode::kb:
         *offset = (int)((char*)&elems[0].kb - (char*)elems);
         *sizeType = sizeof(int);
         break;
      case (int)VarCode::isdry:
         *offset = (int)((char*)&elems[0].isDry - (char*)elems);
         *sizeType = sizeof(int);
         break;            
      case (int)VarCode::w:
         *offset = (int)((char*)&elems[0].w - (char*)elems);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::FaceIndices:
         // Notice that this returns pointers to the nodes of the element, some arithmetic is needed to get the actual indices
         *offset = (int)((char*)&elems[0].nodes[0] - (char*)elems);
         *sizeType = sizeof(Node);
         break;
      default:
         *offset = -1;
   }
}


extern "C" __global__ void getSideVarAddr(int varCode, void **addr, int *offset, int *sizeSide, int *sizeType) {
   Side *sides = mesh->localSideData;

   *addr = sides;
   *sizeSide = sizeof(Side);
   switch (varCode) {
      case (int)VarCode::z:
         *offset = (int)((char*)&sides[0].zc - (char*)sides);
         *sizeType = sizeof(double);
         break;
      // case (int)VarCode::depth:
      //    *offset = (int)((char*)&sides[0].depth - (char*)sides);
      //    *sizeType = sizeof(double);
      //    break;         
      case (int)VarCode::x:
         *offset = (int)((char*)&sides[0].center.x - (char*)sides);
         *sizeType = sizeof(Vector2);
         break;
      case (int)VarCode::len:
         *offset = (int)((char*)&sides[0].len - (char*)sides);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::kb:
         *offset = (int)((char*)&sides[0].kb - (char*)sides);
         *sizeType = sizeof(int);
         break;
      case (int)VarCode::isdry:
         *offset = (int)((char*)&sides[0].isDry - (char*)sides);
         *sizeType = sizeof(int);
         break;            
      // case (int)VarCode::w:
      //    *offset = (int)((char*)&sides[0].w2 - (char*)sides);
      //    *sizeType = sizeof(double);
      //    break;
      case (int)VarCode::u:
         *offset = (int)((char*)&sides[0].u2[0].x - (char*)sides);
         *sizeType = sizeof(Vector2);
         break;   
      case (int)VarCode::btu:
         *offset = (int)((char*)&sides[0].btu[0].x - (char*)sides);
         *sizeType = sizeof(Vector2);
         break;          
      case (int)VarCode::U:
         *offset = (int)((char*)&sides[0].U.x - (char*)sides);
         *sizeType = sizeof(Vector2);
         break;
      case (int)VarCode::gradETP:
         *offset = (int)((char*)&sides[0].gradETP.x - (char*)sides);
         *sizeType = sizeof(Vector2);
         break;
      case (int)VarCode::gradPr:
         *offset = (int)((char*)&sides[0].gradPr.x - (char*)sides);
         *sizeType = sizeof(Vector2);
         break;
      case (int)VarCode::gradEta1:
         *offset = (int)((char*)&sides[0].gradEta1.x - (char*)sides);
         *sizeType = sizeof(Vector2);
         break;
      case (int)VarCode::gradEta2:
         *offset = (int)((char*)&sides[0].gradEta2.x - (char*)sides);
         *sizeType = sizeof(Vector2);
         break;
      // case (int)VarCode::ghat:
      //    *offset = (int)((char*)&sides[0].ghat.x - (char*)sides);
      //    *sizeType = sizeof(Vector2);
      //    break;   
      case (int)VarCode::hhat:
         *offset = (int)((char*)&sides[0].hhat - (char*)sides);
         *sizeType = sizeof(Vector2);
         break;      
      case (int)VarCode::chi:
         *offset = (int)((char*)&sides[0].chi - (char*)sides);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::n1:
         *offset = (int)((char*)&sides[0].n1 - (char*)sides);
         *sizeType = sizeof(Node);
         break;
      case (int)VarCode::n2:
         *offset = (int)((char*)&sides[0].n2 - (char*)sides);
         *sizeType = sizeof(Node);
         break;                   
      default:
         *offset = -1;
   }
}


extern "C" __global__ void getNodeVarAddr(int varCode, void **addr, int *offset, int *sizeNode, int *sizeType) {
   Node *nodes = mesh->localNodeData;

   *addr = nodes;
   *sizeNode = sizeof(Node);

   switch (varCode) {
      case (int)VarCode::depth:
         *offset = (int)((char*)&nodes[0].depth - (char*)nodes);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::eta:
         *offset = (int)((char*)&nodes[0].eta2 - (char*)nodes);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::x:
         *offset = (int)((char*)&nodes[0].p.x - (char*)nodes);
         *sizeType = sizeof(Vector2);
         break;
      case (int)VarCode::z:
         *offset = (int)((char*)&nodes[0].z - (char*)nodes);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::kb:
         *offset = (int)((char*)&nodes[0].kb - (char*)nodes);
         *sizeType = sizeof(int);
         break;
      case (int)VarCode::isdry:
         *offset = (int)((char*)&nodes[0].isDry - (char*)nodes);
         *sizeType = sizeof(int);
         break;         
      case (int)VarCode::w:
         *offset = (int)((char*)&nodes[0].w2 - (char*)nodes);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::u:
         *offset = (int)((char*)&nodes[0].u2[0].x - (char*)nodes);
         *sizeType = sizeof(Vector2);
         break;  
      case (int)VarCode::u_air:
         *offset = (int)((char*)&nodes[0].u_air - (char*)nodes);
         *sizeType = sizeof(Vector2);
         break;
      default:
         *offset = -1;    
   }  
   
}




extern "C" __global__ void getSideBndVar(int varCode, double *var, int k = -1) {
   // Get the value of a variable from an side
   // varName is the name of the variable to get
   // var is the array to store the value
   // k is the vertical level index of the variable to get
   SIDE_BND(sideBnd);   

   int idx = sideBnd->side!=nullptr?sideBnd->side->idx:-1;
   if (sideBnd->segment->isLand) idx = -7;

   if      (varCode == (int)VarCode::idx     && k == -1)   var[idx_side_bnd] = idx;
   else if (varCode == (int)VarCode::ifltype && k == -1)   var[idx_side_bnd] = (sideBnd->segment!=nullptr)?sideBnd->segment->ifltype:-1.1;

}