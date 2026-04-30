// BEGIN INCLUDE: constants.h
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
constexpr double btrack_nudge = 0.009013;   // Nudging factor for backtracking

constexpr int nadv = 1; // Advection flag

constexpr double dfv0 = 0.01; 

constexpr double dzb_min = 0.5; // Min. bottom boundary layer thickness [m]

constexpr double dt = 120;  // Time step  XXX 

constexpr double h0 = 0.0001;  // Minimum water depth for wetting and drying [m]

enum class VarCode {x=0, z=1, u=2, w=3, U=4, nx=5, ny=6, hhat=7, ghat=8, chi=9, ssign=10, btu=11, bcc=12, uth=13, ubc=14, island=15, ifltype=16, elbc=17, isdry=18, eta=19, area=20, qel=21, depth=22, radius=23, len=24, lonlat=25, kb=26, n1=27, n2=28, idx=29, htot=30, gradPr=31, gradEta1=32, gradEta2=33, gradETP=34, FaceIndices=35, u_air=36, outputCoords=37}; // Var codes for getElementVar, getNodeVar...
// END INCLUDE: constants.h
// BEGIN INCLUDE: meshdata.h
#define USE_LOCALIZED_SIGMA
#define USE_FRAMES  // XXXXXXXXXXXXXXX
constexpr double small1 = 1.0e-7;       // Small non-negative number
constexpr double small2 = small1*100.0;  // Slightly larger than small1
constexpr double small = 1.0e-7;            // Small positive number or 0

constexpr int internalSegment = -1;  // It is not really a segment, but we keep the name.
constexpr int landSegment     = -2;

constexpr double pi = 3.14159265358979323846;

constexpr double eps = 1.019e-2; // Nudging factor XXX WHY THIS VALUE?

constexpr double grav        = 9.81;      // [m/s^2]
constexpr double omega_e     = 7.292e-5;  // [1/s]

constexpr double ixi_n[4] = {-1, 1, 1,-1}; // Local  xi-coordinate of the vertices of a quad 
constexpr double iet_n[4] = {-1,-1, 1, 1}; // Local eta-coordinate of the vertices of a quad

#define INT_MAX 2147483647

#ifdef DEBUG
  #define assertFmt(cond, fmt, ...)                                      \
    do {                                                                   \
      if (!(cond)) {                                                       \
        printf("ASSERT FAILED: " fmt "\n", ##__VA_ARGS__);                \
        asm("trap;");                                                      \
      }                                                                    \
    } while (0)
#else
  #define assertFmt(cond, fmt, ...) do { (void)sizeof(cond); } while (0)
#endif

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


__device__ __forceinline__ double clamp(double x, double min, double max) {
    return fmin(fmax(x, min), max);
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



struct Frame {
   // Column 0 (East):  (-sin_lon, cos_lon, 0)              — store 2 values, t20=0 implicit
   // Column 1 (North): Normalized oblate ellipsoid vector  — store 3 values
   // Column 2 (Up):    cross(col0, col1)                   — recomputed on demand
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
// printf("lat = %f, lon = %f, cos_lat = %f, sin_lat = %f, cos_lon = %f, sin_lon = %f\n", lat, lon, cos_lat, sin_lat, cos_lon, sin_lon);
// printf("t00 = %f, t10 = %f, t01 = %f, t11 = %f, t21 = %f\n", t00, t10, t01, t11, t21);

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


   __forceinline__ __device__ Vector2 projHorVec(Vector2 const &vec0, Frame const &frameout) const {
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

__device__ bool celestialWalk(Element *&curElem, Point3 x0, Point3 &xt, const Side *origSide, Side *& lastCrossSide);


struct CubicSpline { 
   double *x   = nullptr;
   double *y   = nullptr;
   double *yp  = nullptr;
   double *ypp = nullptr;
};


struct Element {

   // Geometry information
   int idx;                   // Index of the element
   int i34;                   // Number of sides for the element 
   Side *sides[4];            // Sides of the element
   Node *nodes[4];            // Nodes of the element
   Element *nghbElems[4];     // Neighbors of the element

   int kb;                    // Bottom index of the element
   double area;               // Area of the element
   double ssign[4];           // Sign associated with each side of an element
   double radius;             // Element equivalent radii

   #ifdef USE_FRAMES
   Frame frame;
   #endif

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
         const Point3 local = frame.G2L(p - c);
         return local;
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
   __device__ __forceinline__ void getLocalHorizontalNodePositionsTri (Point2 &x0, Point2 &x1, Point2 &x2, int idx = 0) const;
   __device__ __forceinline__ void getLocalHorizontalNodePositionsQuad(Point2 &x0, Point2 &x1, Point2 &x2, Point2 &x3) const;


   // (eta,xi) coordinates of the Gauss-Legendre quadrature points (idxNode = 0,1,2,3)
   __device__ __forceinline__ Point2 GaussPoint(int idxNode) const {
       constexpr double g = 0.5773502691896257; // 1/sqrt(3)
       assert(idxNode >= 0 && idxNode < 4 && "GaussPoint: idxNode out of range");
       const double eta = (idxNode & 1) ? +g : -g;
       const double xi  = (idxNode & 2) ? +g : -g;
       return {eta, xi};
   }  

   // (eta,xi) coordinates of the Gauss-Legendre quadrature points (idx = 0,1,2,3)
   __device__ __forceinline__ double etaSign(int idxNode) const {
      // bit1: 0 -> -1, 1 -> +1
      assert(idxNode >= 0 && idxNode < 4 && "etaSign: idxNode out of range");
      return (idxNode & 2) ? 1.0 : -1.0;     // [-1,-1,+1,+1]
   }

   __device__ __forceinline__ double xiSign(int idxNode) const {
      // +1 if (bit0 XOR bit1), else -1
      assert(idxNode >= 0 && idxNode < 4 && "xiSign: idxNode out of range");
      return ((idxNode & 1) ^ (idxNode >> 1)) ? 1.0 : -1.0;   // [-1,+1,+1,-1]

   }

   // Shape function Phi at a given point p = (eta,xi), idxNode = node index 0,1,2,3
   __device__ __forceinline__ double Phi(int idxNode, Point2 p) const {
      return 0.25*(1.0 + p.x*etaSign(idxNode))*(1.0 + p.y*xiSign(idxNode));
   }

   // Derivative of shape function Phi_idx/dxi at a given point p = (eta,xi)
   __device__ __forceinline__ double dPhi_dxi(int idxNode, Point2 p) const {
      return 0.25*xiSign(idxNode)*(1.0 + p.x*etaSign(idxNode));
   }   

   // Derivative of shape function Phi_idx/deta at a given point p = (eta,xi) 
   __device__ __forceinline__ double dPhi_deta(int idxNode, Point2 p) const {
      return 0.25*(1.0 + p.y*xiSign(idxNode))*etaSign(idxNode);
   }      
   
   __inline__ __device__ double detJ_eta(void) const;

   __inline__ __device__ double detJ_xi(void) const;

   __inline__ __device__ double jacobianDet(Point2 p) const { return area/4.0 + detJ_eta()*p.x + detJ_xi()*p.y; }

   // Check if the element is dry using the more up-to-date elevation.
   __inline__ __device__ bool isDry2(void) const;

   __inline__ __device__ bool isEdgeBoundary(Side *side) const;

   __device__ bool inverseBilinearQuad(Point2 p, double &xi, double &eta, double shp[4]) const;

   __inline__ __device__ double computeSignedArea(void) const;

   __inline__ __device__ void nudgeToBarycenter(Point2 &p, const Side *origSide = nullptr) const;

   enum EFixIf { fiNo, fiOnlyInside, fiAnywhere };

   __device__ bool computeShapeCoords(Point2 &xt, double shpCoords[4], EFixIf fixNegShapeCoords = fiNo, const Side *origSide = nullptr) const;
   
   __device__ void shapeCoordsToXY(double shpCoords[4], Vector2 &xt) const;

   __device__ double intPhiPhiQuad(int i1, int i2) const;

   __device__ double intGradPhiGradPhiQuad(int i1, int i2) const;

   __device__ Vector2 grad2(double *a) const
   {
      double gx = a[0]*gradPhi[0].x + a[1]*gradPhi[1].x + a[2]*gradPhi[2].x;
      double gy = a[0]*gradPhi[0].y + a[1]*gradPhi[1].y + a[2]*gradPhi[2].y;
      if (isQuad()) {
         gx += a[3]*gradPhi[3].x;
         gy += a[3]*gradPhi[3].y;
      }
      return {gx, gy};
   }

   __device__ Vector2 grad2(double Node::* field) const
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


   #ifdef USE_FRAMES
   __device__ void setFrame(void) {
      double lat, lon;

      const Point3 p = toPoint3(center, centerZ());

      // Compute the latitude and longitude of the element center
      computeLatLon(p, lat, lon);

      frame.update(lat, lon);        
   }
   #endif

   __device__ bool findNextIntersection(Point2 x1, Point2 x2, Side *&side, Element *&nextElement, double &t, 
                                        double& tIntersect, const Side *origSide = nullptr);

   __inline__ __device__ bool containsPoint(Vector2 xt, const Side *origSide) const;

   __inline__ __device__ bool containsPoint(Vector2 &xt, bool fixCoords, const Side *origSide) const;

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



   __device__ void trackGoesOutOfWater(const Side *nextSide, const Side *origSide, Element *curElem, double tt, int jlev,
                                       Point3 &xt, double zint, Vector3 &u_int);

   __device__ __forceinline__ Point2 getLocalNodePosition(int idx, const Side* origSide) const;

   __device__ __forceinline__ Node *getNodeInOrder(int idx) const;

   __device__ __forceinline__ void computePredicates(const Point2& xt, const Side* origSide,
                                                    Bool4& isLeft, Bool4& isObtuse, Bool4& isBisector) const;
    

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

   #ifdef USE_FRAMES
   Frame frame;
   #endif

   SideBoundary *sideBnd; // Pointer to additional boundary information. For non-boundary sides, 
                          // this is nullptr
   
   int kb;                // Bottom index of the side
   bool isDry;            // True if the side is dry
 
   Vector2 n;             // Side unit normal vector (for ics=2: these are done inside local lon/lat frame)

   Tracers tracers[nvrt];   // Tracer concentrations at nodes

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
   __inline__ __device__ Point3 p3_t(double t) const;

   #ifdef USE_FRAMES
   __device__ void setFrame(void) {
      double lat, lon;
      const Point3 p = toPoint3(center, centerZ());

      // Compute the latitude and longitude of the element center
      computeLatLon(p, lat, lon);

      frame.update(lat, lon);        
   }
   #endif

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


   __device__ __forceinline__ Vector2 projectToLocalFrame(Vector2 const& v, Frame const& v_frame) const {
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
 
   __inline__ __device__ bool btrack(int jlev, double dt, double dtbk, Vector3 &u_int) const;




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


};

struct Node {

   // Geometry information
   int idx;        // Index of the node
   Node **nodes;   // Array of pointers to the nodes that are neighbors of this node (there is no guarantee that they are connected)
   int numNodes;   // Number of nodes neighbors of this node

   Element **elements; // Array of pointers to the elements that share this node
   int numElements;    // Number of elements that share this node

   int numSides;       // Number of sides that are connected to this node

   Point2 p;           // Horizontal coordinates at the node. (if ics == 2, this are the x,y coordinates of the 3d global cartesian coordinates. Together with gz, they form the 3D global cartesian coordinates.)
   double gz;          // Global vertical coordinate at the node.
   Vector2 u2[nvrt];
   double  w2[nvrt];
   Vector2 utmp[nvrt];  // Temporary array to store velocities
   double z[nvrt];      // Local vertical coordinate at each vertical level. 
   double lon, lat;
   double fun_lat[ntip + 1]; // XXXXX CHECK IF THIS MAKES SENSE. Also, the +1 is only to avoid erros when ntip ==1
   double depth;

   #ifdef USE_FRAMES
   Frame frame;
   #endif

   #ifdef USE_LOCALIZED_SIGMA   
   double localSigma[nvrt];
   #endif

   int kb;       // Bottom index of the node

   bool isDry;     // True if the node is dry

   int iadv;     // Advection flag

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
   double bdef1; // Bottom deformation at previous timestep
   double bdef2; // Bottom deformation at current timestep

   double qel; // RHS
   double coeffDiag; // Diagonal coefficient of the matrix
   double *coeffMatrix; // Matrix coefficients for neighbors. [0..numNodes-1]

   // This was used by a previous version of the conjugate gradient method
   // double mP; // Search direction for conjugate gradient method
   // double mR; // Residual for conjugate gradient method
   // double mAp; // Product of matrix and search direction for conjugate gradient method

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

            const double cff = 0.001293*dragcoef*speed*rampwind;
            tau = cff*u_air;

         } else {
            tau = {0.0, 0.0};
         }
      }
   }

   #ifdef USE_FRAMES
   __device__ void setFrame(double lat, double lon) {

      frame.update(lat, lon);        
    }
   #endif

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

   __inline__ __device__ double Dz(int k) const { return z[k] - z[k-1]; }
   __inline__ __device__ double htot2(void) const { return eta2 + depth; }
   __inline__ __device__ double zTopCanopy(void) const { return z[kb] + veg_h; } 
   
   __inline__ __device__ Vector3 uvw2(int k) const { return {u2[k].x, u2[k].y, w2[k]}; }

   __inline__ __device__ bool isDry2(void) const { return htot2() <= h0; } // isDry, but updated to use eta2

   __inline__ __device__ bool isBnd(void) const { return nodeBnd != nullptr; } // True if the node is a boundary node
   // __inline__ __device__ int numElements(void) const { return isNbd() ? numElements - 1 : numElements; } // Number of elements that are connected to this node

   __device__ __forceinline__ void clearMatrixRow(void) {
      // Clear the matrix row for this node


      coeffDiag = 0.0;

      for (int i = 0; i < numNodes; i++) coeffMatrix[i] = 0.0;

   }

   __device__ void addToSparseMatrix(Node const *node, double value) {
   // Add the value to the sparse matrix for this node that corresponds to the file of node *this and the column of node *node
      
      if (node == this) {
         coeffDiag += value;
     
      } else {
         // Search (linearly)the index of the node in the array of nodes
         // XXX: This is not efficient, see if we can improve this
         int i;
         for (i = 0; i < numNodes; i++) {
            if (nodes[i] == node) {
               break;
            }
         }
         assert(i < numNodes && "addToSparseMatrix: node not found in array of nodes");
         coeffMatrix[i] += value;
      }
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
   Node    **nodeNeighNodes; // Array of pointers to the nodes neighbors of each node
   Element **nodeNeighElems; // Array of pointers to the elements neighbors of each node

   BoundarySegment landBoundarySegment;
   NodeBoundary landBoundaryNode;

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
                   Node **nodesPtr, Element **elemsPtr, Side **sidesPtr,
                   Node *nodesMem, Element *elemsMem, Side *sidesMem,
                   Node    **localNodeNeighNodesMem,
                   Element **localNodeNeighElemsMem,
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
                       localNodeBoundaryMem(localNodeBoundaryMem), localSideBoundaryMem(localSideBoundaryMem),
                       localSegmentMem(localSegmentMem),
                       h_c(h_c), h_s(h_s), theta_b(theta_b), theta_f(theta_f), sigma(sigma), ztot(ztot),
                       cs(cs), dcs(dcs),
                       coeffMatrixMem(partCoeffMatrixMem),
                       landBoundarySegment(true),
                       landBoundaryNode(true)
   {
      printf("Mesh constructor for partition %d\n", partNum);

      printf("npgb = %d, negb = %d, nsgb = %d, nseg = %d\n", npgb, negb, nsgb, nseg);
      printf("np = %d, ne = %d, ns = %d\n", np, ne, ns);
      printf("max_np = %d, max_ne = %d, max_ns = %d\n", max_np, max_ne, max_ns);
 
      printf("elemsPtr = %p, elems = %p\n", *elems, elemsMem);
      
      sinh_theta_f = sinh(theta_f);

   }

   __device__ void populateNodes(int *i34, int *isbnd, int *sharedIdxNodeBoundary,
                  double *node_coords, double *node_coords_ll, bool use_ll,
                  double *wind_rotate_angle, 
                  int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                  int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                  int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                  int *localNodePart,     int *localElemPart,     int *localSidePart,
                  double *drag, double *roughness, double *rmanning, int *kb, double *localSigma);

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

   __device__ void test(double x, double y, double z)
   {
      // XXXXX REMOVE THIS FUNCTION
      printf("Testing mesh for partition %d, x=%f, y=%f\n", partNum, x, y);
      for (int i = 0; i < ne; i++) {
         auto isInside = localElemData[i].containsPoint({x, y}, nullptr);

         if (isInside) {
            printf("Point %f %f %f is inside element %d\n", x, y, z, i);
            Element *elem = &localElemData[i];
            
            Vector3 u_int;

            // Choose the first side
            Side *side = elem->sides[0];
            int jlev = 5;  
            // RECORDAR:
            // - con y=0 no encuentra el triangulo
            // - jlev = 5, z = 0.1...
            side->btrack(5, 3000.0, 3.1, u_int);

         }
      }
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

#ifdef OLD_PCG_SOLVER
   __global__ void initializeScalars(double *scalars);
 
   __global__ void initializePreconditionedResidual(void);

   __global__ void computeMatrixProductAndDotProducts(void);

   __global__ void updateSolutionAndDotProducts(void);

   __global__ void updateSearchDirection(void);
   
   __global__ void computeHorGradientAtSide(void);  


   __global__ void setupCGSolvers(double* data_d, double* b_d, 
                                  int* indices, int* indptr);

#endif

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
 * Uses the Celestial Walk algorithm to find the element that contains the target point.
 */
int stackSize = 0;

__device__ bool celestialWalk(Element *&curElem, Point3 x0, Point3 &xt, const Side *origSide, Side *&lastCrossSide) {

   constexpr int MAX_ITERATIONS = 200;

   int iter = 0;
   while (iter < MAX_ITERATIONS) {
      iter++;

      Bool4 isLeft, isObtuse, isBisector;
      curElem->computePredicates(toPoint2(xt), origSide, isLeft, isObtuse, isBisector);
// showPoints(x0, xt, jlev, origSide, 0xFF0000u, 'x', 1, 0, 0);
      Bool4 isValid = !isLeft & (!isObtuse | isBisector); // bitwise operation are overloaded

      if (isValid.none()) 
      // Inside the element
      {
         // Just returns true and the element (curElem) is unchanged).
         return true;
      }
      else { // Outside the element

         // Outside this element. Pick the first valid exit side and cross.
         lastCrossSide = curElem->sides[isValid.firstOne()];
         curElem = lastCrossSide->otherElem(curElem);

         // If it crosses into a dry element or a boundary element, return false. This means that the track will now move 
         // along the boundary.
         if (curElem == nullptr || curElem->isDry) return false;
         
      }

   }

   // // First thing is to check if the target point is inside the element. If it is, there is nothing
   // // else to do. This means that the full step occurred within the current element.
   // assert(!isDry && "Starting element is dry");
   // assert((origSide == nullptr || ics == 2) && "origSide must be nullptr if ics == 1 and vice versa");

   // assert(containsPoint(toPoint2(x0), origSide) && "Starting point is not inside the element");

   // if (containsPoint(toPoint2(xt), true, origSide)) {  // The "true" means that, if the point is very close to the boundary, it is nudged to be definitely inside the element.

   //    assert(containsPoint(toPoint2(xt), origSide) && "Target point is not inside the element 1");
   //    assert(!isDry && "Target point is inside a dry element");

   //    u_int = interpolateVelocities(xt, origSide);

   //    nextElem = this;

   //    return true;
   // }


   

//    assert(norm2(toVector2(xt - x0)) > 0.0 && "Zero length track"); 
   
//    double t;
//    Side *nextSide = nullptr;
//    double tIntersect = 1e10;

// // showPoints(x0_local, xt, jlev, origSide, 0xFFFFFFu, '+', 2, 0);
//    // Finds the first intersection with the element sides. It is guaranteed to exist because
//    // (x0,y0) is inside the element but the target point is not.
//    bool found_intersection = findNextIntersection(toPoint2(x0), toPoint2(xt), nextSide, nextElem, t, tIntersect, origSide);


// showPoints(x0, xt, jlev, origSide, 0x0000FFu, 'o', 1, 5);

// const Node *n1 = nodes[0];
// const Node *n2 = nodes[1];
// const Node *n3 = nodes[2];

// Point3 p1 = toPoint3(n1->p, n1->gz);
// Point3 p2 = toPoint3(n2->p, n2->gz);
// Point3 p3 = toPoint3(n3->p, n3->gz);
// p1 = origSide->toLocalFrame(p1);
// p2 = origSide->toLocalFrame(p2);
// p3 = origSide->toLocalFrame(p3);
// int s1 = sides[0]->idx;
// int s2 = sides[1]->idx;
// int s3 = sides[2]->idx;

// if (!found_intersection) { return false;}
// // if (!found_intersection) {if (idx == 15876) printf("111112: No intersection found\n"); return false;}
//    assert(found_intersection && "No intersection found at %d");

//    Point2 xtmp = toPoint2(origSide->toLocalFrame(nextSide->p3_t(t))); // If needed, moves the target point to the intersection point

//    // If the next element is not dry, nudge it to the barycenter to ensure that the point is inside the element.
//    if (nextElem != nullptr && !nextElem->isDry) {
//       assert(nextElem->containsPoint(xtmp, origSide) && "Intersection point is not inside the next element");
// Vector3 xtmp_Before = origSide->toLocalFrame(nextSide->p3_t(t));
//       nextElem->nudgeToBarycenter(xtmp, origSide);
// Vector3 xtmp_After = {xtmp.x, xtmp.y, (1 - tIntersect)*x0.z + tIntersect*xt.z};
// showPoints(xtmp_Before, xtmp_After, jlev, origSide, 0x00FFFFu, '+', 2, 10);
//    }

// // printf("x0 = (%.17f, %.17f, %.17f), xt = (%.17f, %.17f, %.17f), tIntersect = %.17f\n", x0.x, x0.y, x0.z, xt.x, xt.y, xt.z, tIntersect);
//    const double zint = (1 - tIntersect)*x0.z + tIntersect*xt.z;
//    xt = {xtmp.x, xtmp.y, zint};

//    // if the tracking segment goes out of water (either dry or edge boundary)
//    if (isEdgeBoundary(nextSide) || (nextElem == nullptr || nextElem->isDry)) {

      
// Vector3 xt_Before = xt;
//       trackGoesOutOfWater(nextSide, origSide, this, t, jlev, xt, zint, u_int);
// showPoints(xt_Before, xt, jlev, origSide, 0x00FF00u, '+', 14, 5);

//       assert(this->containsPoint(toPoint2(xt), origSide) && "Target point is not inside the current element 2");
//       // When this happens, the backtracking stops.

//       dt *= tIntersect;

//       // In this case, the next element is again the current element.
//       nextElem = this;

// // printf("trackGoesOutOfWater: side->idx = %d, dt = %.17f, time_rm = %.17f, u_int = (%.17f, %.17f, %.17f)\n", nextSide->idx, dt, time_rm, u_int.x, u_int.y, u_int.z);

//       // Returns "early", meaning that we stopped before we were able to reach the target point.
//       // In this case, the caller (btrack) needs to take into account the remaining time.
//       return true;    
//    }


//    assert(nextElem->containsPoint(toPoint2(xt), origSide) && "Target point is not inside the next element 3");
//    // Interpolate velocities at the final position. This is a conceptually simple step,
//    // but the function is not simple, because it has to search for the vertical level containing zt.

   
//    dt *= tIntersect; // reduces the time step to the time it takes to reach the intersection

//    u_int = nextElem->interpolateVelocities(xt, origSide); 

//    assert(nextElem != this && "Next element is the same as the current element");

   return true;
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
   if (*debugPlotIdx >= 10000) ok = false;
   if (ok) {
      debugPlot[(*debugPlotIdx) + 0*10000] = x1.x;
      debugPlot[(*debugPlotIdx) + 1*10000] = x1.y;
      debugPlot[(*debugPlotIdx) + 2*10000] = x2.x;
      debugPlot[(*debugPlotIdx) + 3*10000] = x2.y;
      debugPlot[(*debugPlotIdx) + 4*10000] = (double)color;
      debugPlot[(*debugPlotIdx) + 5*10000] = (double)(((int)(unsigned char)marker & 0xFF) | ((dashLen & 0xFF) << 8) | ((width & 0xFF) << 16));
      debugPlot[(*debugPlotIdx) + 6*10000] = (double)label;
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
   if (origSide->idx == 64501 && jlev >= 28) {
      // Convert points back to original frame
      x1.z = 0.0;
      x2.z = 0.0;
      Point3 X1 = origSide->frame.L2G(x1) + toPoint3(origSide->center, origSide->centerZ());// + origSide->center3D(jlev);
      Point3 X2 = origSide->frame.L2G(x2) + toPoint3(origSide->center, origSide->centerZ());// + origSide->center3D(jlev);
// printf("wwwwww origSide->idx = %d, marker = %c, width = %d, dashLen = %d, x1 = (%.17f, %.17f, %.17f), x2 = (%.17f, %.17f, %.17f)\n", origSide->idx, marker, width, dashLen, x1.x, x1.y, x1.z, x2.x, x2.y, x2.z);

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

      assert(level >= 0 && level < nvrt && "Level is out of range");
      // Get the values of the field at the nodes of the element
      res[0] = (nodes[0]->*field)[max(level, nodes[0]->kb)];
      res[1] = (nodes[1]->*field)[max(level, nodes[1]->kb)];
      res[2] = (nodes[2]->*field)[max(level, nodes[2]->kb)];
      if (isQuad()) res[3] = (nodes[3]->*field)[max(level, nodes[3]->kb)];
      else          res[3] = T();
   }    

   template<typename T>
   __device__ __forceinline__ void Element::getNodalValues(int level, T (Node::*func)(int) const, T res[4]) const {

      assert(level >= 0 && level < nvrt && "Level is out of range");

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

      assert(level >= 0 && level < nvrt && "Level is out of range");

      // Get the values of the field at the nodes of the element
      res[0] = (sides[0]->*field)[max(level, sides[0]->kb)];
      res[1] = (sides[1]->*field)[max(level, sides[1]->kb)];
      res[2] = (sides[2]->*field)[max(level, sides[2]->kb)];
      if (isQuad()) res[3] = (sides[3]->*field)[max(level, sides[3]->kb)];
      else          res[3] = T();
   }    

   template<typename T>
   __device__ __forceinline__ void Element::getSideValues(int level, T (Side::*func)(int) const, T res[4]) const {

      assert(level >= 0 && level < nvrt && "Level is out of range");

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
printf("mean 2: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
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
printf("mean 2: values = (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", values[0].x, values[0].y, values[1].x, values[1].y, values[2].x, values[2].y, values[3].x, values[3].y);
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


__device__ __forceinline__ void Element::getLocalHorizontalNodePositionsTri(Point2 &x0, Point2 &x1, Point2 &x2, int idx) const {
// Performs "getLocalHorizontalNodePosition" for the three nodes of the triangle.

   assert(isTri() && "getLocalHorizontalNodePositions only works for triangles");

   const int id1 = idx;
   const int id2 = nxq(1 + id1);
   const int id3 = nxq(2 + id1);

   if constexpr(ics==1) {
      x0 = nodes[id1]->p;
      x1 = nodes[id2]->p;
      x2 = nodes[id3]->p;

   } else {
      const Vector3 c = {center.x, center.y, centerZ()}; // Element center

      const Point3 p0 = toPoint3(nodes[id1]->p, nodes[id1]->gz); // XXXXX use p3D or similar?
      const Point3 p1 = toPoint3(nodes[id2]->p, nodes[id2]->gz);
      const Point3 p2 = toPoint3(nodes[id3]->p, nodes[id3]->gz);

      x0 = toVector2(frame.G2L(p0 - c));  // Vectors are implicitly cast into points at the local frame centered at c
      x1 = toVector2(frame.G2L(p1 - c));
      x2 = toVector2(frame.G2L(p2 - c));

   }

}

__device__ __forceinline__ void Element::getLocalHorizontalNodePositionsQuad(Point2 &x0, Point2 &x1, Point2 &x2, Point2 &x3) const {
// Performs "getLocalHorizontalNodePosition" for the four nodes of the quad.

   assert(isQuad() && "getLocalHorizontalNodePositionsQuad only works for quads");

   if constexpr(ics==1) {
      x0 = nodes[0]->p;
      x1 = nodes[1]->p;
      x2 = nodes[2]->p;  
      x3 = nodes[3]->p;

   } else {
      const Vector3 c = {center.x, center.y, centerZ()}; // Element center

      const Point3 p0 = toPoint3(nodes[0]->p, nodes[0]->gz); // XXXXX use p3D or similar?
      const Point3 p1 = toPoint3(nodes[1]->p, nodes[1]->gz);
      const Point3 p2 = toPoint3(nodes[2]->p, nodes[2]->gz);
      const Point3 p3 = toPoint3(nodes[3]->p, nodes[3]->gz);

      x0 = toVector2(frame.G2L(p0 - c));  // Vectors are implicitly cast into points at the local frame centered at c
      x1 = toVector2(frame.G2L(p1 - c));
      x2 = toVector2(frame.G2L(p2 - c));
      x3 = toVector2(frame.G2L(p3 - c));
   }

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

      // Computes the cross product of vectors (n2 - n0) and (n1 - n0). 
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

/**
   For a point on a side of an element and a line segment starting at that point and ending at (x2,y2), 
   find the intersection point (given by parameter tIntersect) and the side (intersectingSide) that it intersects.
   If the segment is completely inside the element, then it returns true.

   This function assumes:
   - The initial point is on a side of the element. 
   - Its position is given by parameter t in the range [0,1].
   - The line segment starting at that point and ending at (x2,y2) is either completely inside the element or it 
      intersects the element (elements are convex).

   To find the intersection point, this function does the following:
   - it finds the intersection of the line segment with all the element side lines,
      except for the one that contains the initial point.
   - Because we do the interection with the lines (not the segments), the t parameter is not ncessearily in the range [0,1].
   - let's call t1i the t parameter for the track segment and t2i the t parameter for the element side line.
   - We want the closest intersection point, for which t1i > 0 and t2i in [0,1].
   - From construction, we know that this intersection point *must* exist, even in extreme cases, like when the
      track segment is almost parallel to the element side line.

   Inputs:
   @param t: t*side.n1 + (1-t)*side.n2 is the point on the side
   @param x2, y2: other end of the line segment

   Outputs:
   @param side:  side that intersects the line segment (input and output)
   @param nextElement (output): element that contains the intersection point
   @param tIntersect (output): New intersection point tIntersect*side.n1 + (1-tIntersect)*side.n2
   @param zint (input/output): Vertical coordinate of the intersection point

   The equations for the two segments (the one we are considering and one of the sides of the element) are:
   P1(tt1) = (x1,y1) + tt1*(x2-x1, y2-y1)  where 0 <= tt1 <= 1
   P2(tt2) = (x3,y3) + tt2*(x4-x3, y4-y3)  where 0 <= tt2 <= 1

   If the two segments are parallel, they do not intersect. This is checked by the determinant of the two vectors
   (Delta) being zero.

   The intersection2 point happens when:
   tt1(x2-x1) - tt2(x4-x3) = x3-x1
   tt1(y2-y1) - tt2(y4-y3) = y3-y1
   
   Assumption: the 4 points (including those of the side segment) are different.
   The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.	
   Output: false: no intersection or colinear; true: exactly 1 intersection.	
   If true, (xIntersect,yIntersect) is the intersection.				
   Modified to only check tt2, assuming 0<=tt1<=1 is already assured. 
**/


__device__ bool Element::findNextIntersection(Point2 x1, Point2 x2, Side *&side, Element *&nextElement,
                                              double &t, double& tIntersect, const Side *origSide) {
   assert((origSide == nullptr || ics == 2) && "origSide must be nullptr if ics == 1 and vice versa (1)");
// REMOVE!!!!
//    assert( containsPoint(x1, origSide) && "Starting point is not inside the element");                                    
//    assert(!containsPoint(x2, origSide) && "Ending point is not outside the element");

//    double closest_intersect_t = 1e10;
//    bool found_intersection = false;

//    for (int i = 0; i < i34; i++) {

//       Side *side2 = sides[i];

//       // (x3,y3)-(x4,y4) is the element side line
//       Point2 x3 = side2->n1->p;
//       Point2 x4 = side2->n2->p;

//       if constexpr(ics==2) {
         
//          x3 = toPoint2(origSide->toLocalFrame(toPoint3(x3, side2->n1->gz))); // XXXX?????
//          x4 = toPoint2(origSide->toLocalFrame(toPoint3(x4, side2->n2->gz)));  // XXXX?????

//       }

//       // delta is 0 if the two lines are parallel. In this case, they do not intersect.
//       const double delta = cross(x2 - x1, x3 - x4);

//       if (abs(delta) > 0) {

//          const double delta1 = cross(x3 - x1, x3 - x4);
//          const double delta2 = cross(x2 - x1, x3 - x1); 

//          // Solves the two linear equations for tt1 and tt2
//          const double tt1 = delta1/delta; // delta is not zero, so we can divide
//          const double tt2 = delta2/delta;

//          if (tt1 >= -small && tt2 >= -small && tt2 <= 1.0 + small && tt1 <= 1.0 + small) {
//             closest_intersect_t = min(closest_intersect_t, tt1);
//             side = side2;
//             nextElement = side->otherElem(this);

//             t =          clamp(tt2, small, 1.0-small);
//             tIntersect = clamp(tt1, small, 1.0-small);
//             found_intersection = true;

//          }

      
//       }
//       else {  // parallel segments

//          // continue; 
//          assert(false && "Parallel segments should not happen under the new code");
//          // In this case, we check if they are colinear. if they are, we compute the intersection point as:
//          // tt2 += length(x1 - x2) * sign(dot(x1 - x2, x3 - x4)) in {0, 1}

//          // Check if the segments are colinear. It is enough to check if x1 is on the line of the side segment.
//          const double dist = dot(x1 - x3, side2->n);
//          if (abs(dist) < small*length(x3 - x4)) {
//             const double2 unitarySide = normalize(x4 - x3);
//             const double len = dot(x2 - x1, unitarySide);
//             t = t + len;
//             t = max(0.0, min(t, 1.0));
// // printf("BEEEEEBB t = %.17e, dist = %.17e, len = %.17e, x1 = (%e,%e), x2 = (%e,%e), x3 = (%e,%e), x4 = (%e,%e)\n", t, dist, len, x1.x, x1.y, x2.x, x2.y, x3.x, x3.y, x4.x, x4.y);


//             // tIntersect = max(small, min(tt1, 1.0-small));
//             // t = tt2;
//             found_intersection = true;
//          }
//       }


//    }



//    // if (!found_intersection) {if (idx == 15876) printf("11111: No intersection found\n"); return false;}
//    assert(found_intersection && "No intersection found");

//    // found_intersection should always be true, since we checked all the sides.  
//    return found_intersection;
}

__inline__ __device__ bool Element::containsPoint(Vector2 xt, const Side *origSide) const { 
   // Check if the point (xt,yt) is inside the element.
   // if ics==2, origSide is used to convert the point to the starting point's frame.
   assert((ics==2 || origSide==nullptr) && "origSide must be nullptr if ics == 1 and vice versa (2)");

   double shpCoords[4];
   bool res = computeShapeCoords(xt, shpCoords, fiNo, origSide);

   return res;
}

__inline__ __device__ bool Element::containsPoint(Vector2 &xt, bool fixCoords, const Side *origSide) const { 
   // Check if the point (xt,yt) is inside the element.
   // If fixCoords is true, and *the point is inside the element*, the point is nudged to be inside the element.
   // if ics==2, origSide is used to convert the point to the starting point's frame.
   assert((ics==2 || origSide==nullptr) && "origSide must be nullptr if ics == 1 and vice versa (3)");

   double shpCoords[4];
   Vector2 tmp_xt = xt;
   const EFixIf fixIfCoords = fixCoords ? fiOnlyInside : fiNo;
   bool res = computeShapeCoords(tmp_xt, shpCoords, fixIfCoords, origSide);

   if (fixCoords && res) xt = tmp_xt;

   return res;
}


__device__ void Element::trackGoesOutOfWater(const Side *nextSide, const Side *origSide, Element *curElem, double tt, int jlev,
                                             Point3 &xt, double zint, Vector3 &u_int) {

   // When this happens, we stop our algorithm that usually assumes that we move to the next element.
   // Instead, we nudge the intersection point to the barycenter of the current element, and nextElement becomes
   // curElement.

   // Then it finds the tangential velocity at the element side (the side where the intersection was found),
   // And uses this velocity instead to continue the tracking.

   // computes the intersection point in local cartesian coordinates
   Point3 intersection = origSide->toLocalFrame(nextSide->p3_t(tt));  // XXXX I don't think this is necessary. The intersection is already computed in findNextIntersection.

   // Nudges the point to the barycenter of the element. 
   // This ensures that the point is inside the element and hence, we will find an intersection if the
   // next segment intersects the element sides. 
   // TODO: curElem->nudgeToBarycenter(intersection);  REMEMBER!!! RECOVER??? XXX KKKKK


   // Computes the tangential velocity at the side. Note that we are not using the x,y,z location
   // of the intersection point. We are approximating it by the velocity at the side center at jlev.
   // NOTE: With the new code it will be easy to compute the tangential velocity at the intersection point.
   double vtan = dot(nextSide->t(), nextSide->u2[jlev]);

   // If the side is an open boundary, the tangential velocity is set to zero.
   if (nextSide->isBnd()) vtan = 0.0;

   // The new velocity components are computed from the tangential velocity and the tangent vector.
   const double us = vtan*nextSide->t().x;
   const double vs = vtan*nextSide->t().y;
   const double ws = nextSide->w2(jlev); 

   u_int = {us, vs, ws};

   // The end of this tracking segment is now the intersection point.
   xt = {intersection.x, intersection.y, zint};
// printf("trackGoesOutOfWater: intersection = (%.17f, %.17f, %.17f), zint = %.17f, u_int = (%.17f, %.17f, %.17f)\n", intersection.x, intersection.y, zint, u_int.x, u_int.y, u_int.z);
   
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

__device__ __forceinline__ void Element::computePredicates(const Point2& xt, const Side* origSide,
                                Bool4& isLeft, Bool4& isObtuse, Bool4& isBisector) const
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

   // Side I goes from corners[I] to corners[nxq(I+1)] in CCW perimeter order.
   // (Note: with corners built from getNode in CCW order, no ssign factor is needed
   //  here — the orientation is baked into the corners array.)
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
   isLeft.set(0,  cross(s0, xt - C[ 0]) > -small);
   isLeft.set(1,  cross(s1, xt - C[ 1]) > -small);
   isLeft.set(2,  cross(s2, xt - C[ 2]) > -small);
   isLeft.set(3, (cross(s3, xt - C[i3]) > -small) || isTri());  

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
   isBisector.set(0, dot(xt - C[ 1], s0 + s1) < 0.0);
   isBisector.set(1, dot(xt - C[ 2], s1 + s2) < 0.0);
   isBisector.set(2, dot(xt - C[i3], s2 + s3_tri) < 0.0);
   isBisector.set(3, dot(xt - C[ 0], s3 + s0) < 0.0);   

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

      res += ar1 + ar2 + ar3;
   }

   return res;
}


__inline__ __device__ void Element::nudgeToBarycenter(Point2 &p, const Side *origSide) const {
// Nudges the local-frame point, p, to the barycenter of the element. 

   assert((origSide == nullptr || ics == 2) && "origSide must be nullptr if ics == 1 and vice versa (4)");

   Point3 c = origSide->toLocalFrame(toPoint3(center, centerZ())); // In local frame coordinates if ics == 2

   p = (1.0 - btrack_nudge)*p + btrack_nudge*toPoint2(c);
}

__inline__ __device__ bool Element::isDry2(void) const { return (nodes[0]->isDry2()) ||
                                                                (nodes[1]->isDry2()) ||
                                                                (nodes[2]->isDry2()) ||
                                                                (isQuad() && nodes[3]->isDry2());
}

__inline__ __device__ double Element::detJ_eta(void) const {
   return (nodes[1]->p.x - nodes[0]->p.x)*(nodes[2]->p.y - nodes[3]->p.y) 
        - (nodes[2]->p.x - nodes[3]->p.x)*(nodes[1]->p.y - nodes[0]->p.y);
}

__inline__ __device__ double Element::detJ_xi(void) const {
   return (nodes[2]->p.x - nodes[1]->p.x)*(nodes[3]->p.y - nodes[0]->p.y) 
        - (nodes[3]->p.x - nodes[0]->p.x)*(nodes[2]->p.y - nodes[1]->p.y);
}

__inline__ __device__ bool Element::isEdgeBoundary(Side *side) const { 
   // In boundary sides, the element 1 is nullptr
   return side->elems[1] == nullptr; 
}

__device__ bool Element::inverseBilinearQuad(Point2 p, double &xi, double &eta, double shp[4]) const {
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
   // cross(nb, ·) projects onto the direction PERPENDICULAR to W.
   // Since nb || W, cross(nb, vb) = 0 exactly, killing the xi*eta term.
   // What remains is a pure affine relationship between xi and eta.
   const double Lp = cross(nb, dp4);     // cross(nb, 4*H)
   const double Ls = cross(nb, vxi);     // cross(nb, 4*U)
   const double Lt = cross(nb, veta);    // cross(nb, 4*V)

   // ---- Quadratic equation (Q): Qp = Qs*xi + Qt*eta + Qw*xi*eta ----
   // dot(nb, ·) projects onto the W direction, where the nonlinearity lives.
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
   xi  = fmin(1.0, fmax(xi,  -1.0));
   eta = fmin(1.0, fmax(eta, -1.0));

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

   assert((origSide == nullptr || ics == 2) && "origSide must be nullptr if ics == 1 and vice versa (5)");

   if (isTri()) { // Triangular element

      // Horizontal local node coordinates. 
      Point2 x0 = nodes[0]->p;
      Point2 x1 = nodes[1]->p;
      Point2 x2 = nodes[2]->p;

      // If ics == 2, the node coordinates are transformed to local frame.
      if constexpr(ics==2) {
         const Point3 p0 = toPoint3(x0, nodes[0]->gz);
         const Point3 p1 = toPoint3(x1, nodes[1]->gz);
         const Point3 p2 = toPoint3(x2, nodes[2]->gz);


         // Convert the node coordinates to the local frame centered at the element center
         x0 = toPoint2(origSide->toLocalFrame(p0));
         x1 = toPoint2(origSide->toLocalFrame(p1));
         x2 = toPoint2(origSide->toLocalFrame(p2));
      }

      shpCoords[0] = signedArea(xt,x1,x2)/area;
      shpCoords[1] = signedArea(x0,xt,x2)/area; 
      shpCoords[2] = 1.0 - shpCoords[0] - shpCoords[1];
      shpCoords[3] = 0.0;


      double tmpmin = min3(shpCoords[0], shpCoords[1], shpCoords[2]);
      const bool res = (tmpmin >= -small);  
      
      // If the point is not inside the element (with a small tolerance), we do not need to fix the negative shape coordinates.
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

      // Check if the point is inside the element
      return res;      
   }
   else {
      assert(isQuad() && "computeShapeCoords");

   //    Point2 p = xt;

   //    // Reproject p the original side frame for ics=2
   //    if constexpr(ics == 2) p = origSide->frame.G2L(toPoint3(xt, 0.0));

   //    bool inside = false;



   //    if(ifl==0) then
   //      do i=1,4
   //        in1=nxq(1,i,i34(ie))
   //        in2=nxq(2,i,i34(ie))
   //        swild2(i)=signa2(xel(in1,ie),xel(in2,ie),x,yel(in1,ie),yel(in2,ie),y)
   //      enddo !i
   //      tmp=minval(swild2(1:4))/area(ie)
   //      if(tmp>-small2) then
   //        inside=1
   //        call ibilinear(itag,ie,area(ie),xel(1,ie),xel(2,ie),xel(3,ie),xel(4,ie), &
   //   &yel(1,ie),yel(2,ie),yel(3,ie),yel(4,ie),x,y,xi,eta,shapef,icaseno)
   //      endif !inside quad
   //    else !ifl=1 - already inside
   //      call ibilinear(itag,ie,area(ie),xel(1,ie),xel(2,ie),xel(3,ie),xel(4,ie), &
   //   &yel(1,ie),yel(2,ie),yel(3,ie),yel(4,ie),x,y,xi,eta,shapef,icaseno)
   //      !Update pt
   //      x=dot_product(xel(1:4,ie),shapef(1:4))
   //      y=dot_product(yel(1:4,ie),shapef(1:4))
   //    endif !ifl


   //    return quad_shape(1,6,nnel,xn2,yn2,inside2,arco)

   }

   assert(false && "computeShapeCoords: should never reach this point");
   return false;  // Never reached

}


__device__ void Element::shapeCoordsToXY(double shpCoords[4], Vector2 &xt) const {
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

__device__ double Element::intPhiPhiQuad(int i1, int i2) const {
   // return \int \phi_i1*\phi_i2 dA, where phi_i1 and phi_i2 are the shape functions of the nodes i1 and i2.
   // The integral is computed analytically, using the shape functions in xi and eta coordinates and the Jacobian.
   // The formula for the integral is:
   // (area/16)(1+\xi_i\xi_j/3)(1+\eta_i\eta_j/3) + (detJ_xi/12)(1+\eta_i\eta_j/3)(\xi_i+\xi_j) + (detJ_eta/12)(1+\xi_i\xi_j/3)(\eta_i+\eta_j)\]

   printf("intPhiPhi: TEST\n");

   const Node *n0 = nodes[0];
   const Node *n1 = nodes[1];
   const Node *n2 = nodes[2];
   const Node *n3 = nodes[3];

   const Vector2 p0 = n0->p;
   const Vector2 p1 = n1->p;
   const Vector2 p2 = n2->p;
   const Vector2 p3 = n3->p;
assert(false && "intPhiPhiQuad:CHECK?");

   double xi_term  = 1.0 + (ixi_n[i1]*ixi_n[i2])/3.0;
   double eta_term = 1.0 + (iet_n[i1]*iet_n[i2])/3.0;
   double eta_sum  = iet_n[i1] + iet_n[i2];
   double xi_sum   = ixi_n[i1] + ixi_n[i2];

   return (1.0/16.0)*area*xi_term*eta_term + (detJ_xi()/12.0)*xi_term*eta_sum + (detJ_eta()/12.0)*eta_term*xi_sum;
}
   

__device__ double Element::intGradPhiGradPhiQuad(int i1, int i2) const {
   // return \int (grad phi_i1).(grad phi_i2) dA, where phi_i1 and phi_i2 are the shape functions of the nodes i1 and i2.
   // The integral is computed by using the 2x2 Gauss-Legendre quadrature points in xi and eta coordinates.
   // Again, it uses the determinant of the Jacobian to compute the integral in x,y coordinates.
   // The formula is also based on the identity:
   //   dPhi/dx = (1/det J) * [ (dPhi/dxi)*(dy/deta) - (dPhi/deta)*(dy/dxi)]
   //   dPhi/dy = (1/det J) * [-(dPhi/dxi)*(dx/deta) + (dPhi/deta)*(dx/dxi)]

   printf("intGradPhiGradPhi: TEST\n");   
assert(false && "intGradPhiGradPhi:Do ICS==2");
   assert(isQuad() && i1 >= 0 && i2 >= 0 && i1 < 4 && i2 < 4 && "intGradPhiGradPhi: not quad");

   const double x0 = nodes[0]->p.x;
   const double y0 = nodes[0]->p.y;
   const double x1 = nodes[1]->p.x;
   const double y1 = nodes[1]->p.y;
   const double x2 = nodes[2]->p.x;
   const double y2 = nodes[2]->p.y;
   const double x3 = nodes[3]->p.x;
   const double y3 = nodes[3]->p.y;

   double res = 0.0;
   
   // Quadrature points in xi and eta coordinates
   #pragma unroll
   for (int i = 0; i < 4; i++) {
      // For each quadrature point...

      const Point2 p = GaussPoint(i);
      const double eta = p.x;
      const double xi  = p.y;

      // Jacobian entries at the Gauss point p
      const double dx_dxi  = 0.25 * ( (x1 + x2 - x0 - x3) + xi *(x0 + x2 - x1 - x3) );
      const double dx_deta = 0.25 * ( (x2 + x3 - x0 - x1) + eta*(x0 + x2 - x1 - x3) );
      const double dy_dxi  = 0.25 * ( (y1 + y2 - y0 - y3) + xi *(y0 + y2 - y1 - y3) );
      const double dy_deta = 0.25 * ( (y2 + y3 - y0 - y1) + eta*(y0 + y2 - y1 - y3) );      

      const double detJ = jacobianDet(p);
      assert(detJ > 0.0 && "Negative or zero determinant of the Jacobian in intGradPhiGradPhi");

      // Derivatives of the shape functions i1, i2 at the Gauss point p
      const double dPhi_dxi_1  = dPhi_dxi (i1, p);
      const double dPhi_deta_1 = dPhi_deta(i1, p);
      const double dPhi_dxi_2  = dPhi_dxi (i2, p);
      const double dPhi_deta_2 = dPhi_deta(i2, p);

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

   assert((origSide == nullptr || ics == 2) && "origSide must be nullptr if ics == 1 and vice versa (6)");

   // First step is to compute the shape coordinates of (xt,yt)
// printf("interpolateVelocities 1: xt = (%f, %f, %f)\n", xt.x, xt.y, xt.z);
   computeShapeCoords(toVector2(xt), shpCoord, fiOnlyInside, origSide);

if ((min4(shpCoord[0], shpCoord[1], shpCoord[2], shpCoord[3]) < -small)) printf("idx = %d, shpCoord = (%.16f, %.16f, %.16f, %.16f)\n", idx, shpCoord[0], shpCoord[1], shpCoord[2], shpCoord[3]);
   assert((min4(shpCoord[0], shpCoord[1], shpCoord[2], shpCoord[3]) >= -1e-4) && "Invalid shape coordinates");

   // Second step is to find the level that contains zt
   // Uses a binary search to find the level
   // Note that z coordinates are interpolated from their values at the nodes of the element.
   int L = kb;   // I use L(eft) and R(ight) as usual, but here L is the bottom level and R is the top level
   int R = nvrt - 1;

   // Get the z coordinates at the nodes of the element at the levels L and R
   double zLAtNodes[4], zRAtNodes[4], zMAtNodes[4];
   getNodalValues(L, &Node::compZCoord, zLAtNodes);
   getNodalValues(R, &Node::compZCoord, zRAtNodes);

   double zL = interpolate(shpCoord, zLAtNodes);
   double zR = interpolate(shpCoord, zRAtNodes);
   
   // Interpolates vertically the velocities
   double wght1, wght2;
   if (xt.z<=zL) {
      wght1 = 1.0;
      wght2 = 0.0;
      xt.z = zL;
   }
   else if (xt.z>=zR) {
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

if (R< 0 || R>49) printf("KKKK R = %d, L = %d, wght1 = %.17f, wght2 = %.17f\n", R, L, wght1, wght2);
   // Interpolate horizontally and vertically the velocities at the two levels L, R
   // Remember that interpolate takes frames into account for Vector2 fields and ics==2
   toVector2(u_int) = wght1*interpolate(L, shpCoord, &Node::u2) + wght2*interpolate(R, shpCoord, &Node::u2);
   u_int.z          = wght1*interpolate(L, shpCoord, &Node::w2) + wght2*interpolate(R, shpCoord, &Node::w2);

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
__inline__ __device__ Point3  Side::p3_t(double t) const { return {x(t), y(t), z(t)}; }
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
 * @return: True if the backtracking was successful, false otherwise
 */
__inline__ __device__ bool Side::btrack(
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


   // Initial velocities
   u_int = {u2[jlev].x, u2[jlev].y, w2(jlev)};

   // Find a non dry element for starting the backtracking. At this point there should be at least one.
   Element* curElem = (!elems[0]->isDry) ? elems[0] : elems[1];  // elems[0] is always non-null 

   assert(curElem != nullptr && !curElem->isDry && "Starting element is dry or boundary side");

// showPoints(x0, x0, jlev, this, 0xFF0000u, 'x', 1, 0, 0);


   const int MAX_ITERATIONS = 200; 
   Element* nextElem = nullptr;
   int iter = 0;
   Point3 xt = {-1.0e25, -1.0e25, -1.0e25};
   while(true) {

      // Calculate timestep
      double dtb = min(dtbk, time_rm);

      // Next point in the backtracking curve
      xt = x0 - dtb*u_int;

// showPoints(x0, xt, jlev, this, 0xFF00FFu, 'o', 1, 15, iter);

      Side *lastCrossSide = nullptr;
      const bool found = celestialWalk(curElem, x0, xt, this, lastCrossSide);

      // This means that the walk reached a boundary or a dry element.
      // This case is addressed at the end of the loop.
      if (!found) {

         // FIND CROSSING POINT xtp ON BOUNDARY SEGMENT
         // COMPUTE DISTANCE dtb = distance(x0, xtp) / length(u_int).
         // time_rm -= dtb;
         // MOVE time_rm along the boundary.
         time_rm = 0.0;
         u_int = {0.0, 0.0, 0.0};
         return true;
      }

      time_rm -= dtb; // if celestialWalk terminates early (if it intersects a boundary or a dry element),
                      // then time_rm needs to be increased by the unused time.

      u_int = nextElem->interpolateVelocities(xt, this);

      if (time_rm == 0.0) break; // from the definition of dtb, time_rm will reach exactly 0.0                      

      // If the velocity is too small, we stop the backtracking
      const double vmag = length(toVector2(u_int));
      if (vmag <= velmin_btrack) break;

      // The current point becomes the next segment's starting point.
      x0 = xt;
 
      curElem = nextElem;

      iter++;

      assert(iter < MAX_ITERATIONS && "Max iterations reached in btrack");
   }

   assert(curElem != nullptr && "Final element is nullptr");
   assert(!curElem->isDry && "Final element is dry");

// const Point3 p0 = toLocalFrame(toPoint3(curElem->nodes[0]->p, curElem->nodes[0]->gz));
// const Point3 p1 = toLocalFrame(toPoint3(curElem->nodes[1]->p, curElem->nodes[1]->gz));
// const Point3 p2 = toLocalFrame(toPoint3(curElem->nodes[2]->p, curElem->nodes[2]->gz));
// Vector3 xt_Before = xt;
// curElem->nudgeToBarycenter(toPoint2(xt), this);
// Vector3 xt_After = xt;
// showPoints(xt_Before, xt_After, jlev, this, 0x00006Fu, '+', 2, 10, iter);
// if (!curElem->containsPoint(toPoint2(xt), this)) printf("side->idx = %d, (%.16f,%.16f),(%.16f,%.16f),(%.16f,%.16f)  xt = (%.17f, %.17f, %.17f)\n", idx, p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, xt.x, xt.y, xt.z);



   // Interpolate the velocities at the final point
   u_int = curElem->interpolateVelocities(xt, this);

   // if (idx == 18867) printf("iter = %d, side->idx = %d, curElem->idx = %d, x0 = (%.17f, %.17f, %.17f), u_int = (%.17f, %.17f, %.17f)\n", iter, idx, curElem->idx, x0.x, x0.y, x0.z, u_int.x, u_int.y, u_int.z);         


}



// /**
//    For a point (x1,y1) inside the element and a line segment (x1,y1)-(x2,y2), find the intersection point 
//    (xIntersect,yIntersect) and the side (intersectingSide) that it intersects.

//    Inputs:
//    @param x1, y1: point inside the element
//    @param x2, y2: other end of the line segment

//    Outputs:
//    @param intersectingSide:  side that intersects the line segment
//    @param nextElement: element that contains the intersection point
//    @param xIntersect, yIntersect: intersection point  

//    The equations for the two segments (the one we are considering and one of the sides of the element) are:
//    P1(tt1) = (x1,y1) + tt1*(x2-x1, y2-y1)  where 0 <= tt1 <= 1
//    P2(tt2) = (x3,y3) + tt2*(x4-x3, y4-y3)  where 0 <= tt2 <= 1

//    If the two segments are parallel, they do not intersect. This is checked by the determinant of the two vectors
//    (Delta) being zero.

//    The intersection point happens when:
//    tt1(x2-x1) - tt2(x4-x3) = x3-x1
//    tt1(y2-y1) - tt2(y4-y3) = y3-y1
   
//    Assumption: the 4 points (including those of the side segment) are different.
//    The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.	
//    Output: false: no intersection or colinear; true: exactly 1 intersection.	
//    If true, (xIntersect,yIntersect) is the intersection.				
//    Modified to only check tt2, assuming 0<=tt1<=1 is already assured. 
// **/
// __inline__ __device__ bool Element::findPathIntersection(double x1, double y1, double x2, double y2, 
//                                                          Side& intersectingSide, Element& nextElement,
//                                                          double& xIntersect, double& yIntersect) {
//    for (int i = 0; i < i34; i++) {
//       const double x3 = nodes[i]->x;
//       const double y3 = nodes[i]->y;
//       const double x4 = nodes[nxq(i+1)]->x;
//       const double y4 = nodes[nxq(i+1)]->y;

//       const Side& side = sides[i];

//       xIntersect=-1.e25;
//       yIntersect=-1.e25;

//       const double delta  = (x2-x1)*(y3-y4)-(y2-y1)*(x3-x4);
//       const double delta1 = (x3-x1)*(y3-y4)-(y3-y1)*(x3-x4);
//       const double delta2 = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1); 

//       if (delta != 0.0) {
//          // Solves the two linear equations for tt1 and tt2
//          const double tt1 = delta1/delta; // delta is not zero, so we can divide
//          const double tt2 = delta2/delta;

//          // It is enough to check tt2, since the range [0,1] is the same for both tt1 and tt2.
//          // Similarly, we use tt2 to compute the intersection point.
//          if (tt2 >= -small && tt2 <= 1.0 + small) {
//             xIntersect = x3 + (x4 - x3)*tt2;
//             yIntersect = y3 + (y4 - y3)*tt2;

//             intersectingSide = side;
//             nextElement = nghbElems[i];

//             return true;
//          }
      
//       }
//       else {
//          // The two segments are parallel
//          return false;
//       }

//    }

//    return false;

// }



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

   return -1;
}


__device__ void Mesh::populateNodes(int *i34, int *isbnd, int *sharedIdxNodeBoundary,
                                    double *node_coords, double *node_coords_ll, bool use_ll,
                                    double *wind_rotate_angle, 
                                    int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                                    int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                                    int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                                    int *localNodePart,     int *localElemPart,     int *localSidePart,
                                    double *drag, double *roughness, double *rmanning, int *kb, double *localSigma) {

   NODE(pnode);  
   const int i = idx_node;
   Node &node(*pnode);

   node.idx = i;

   const int idxNodeStart = nodeNeighNodesStart[i];
   const int idxNodeEnd   = nodeNeighNodesStart[i+1];
   node.numNodes = idxNodeEnd - idxNodeStart;
   node.nodes = nodeNeighNodes + idxNodeStart;
   node.coeffMatrix = coeffMatrixMem + idxNodeStart;



   const int idxElemStart = nodeNeighElemsStart[i];
   const int idxElemEnd   = nodeNeighElemsStart[i+1];      
   node.numElements = idxElemEnd - idxElemStart;
   node.elements = nodeNeighElems + idxElemStart;


   for (int idx = idxNodeStart; idx < idxNodeEnd; idx++) {
      const int idxNodeGlobal = nodeNeighNodesIndices[idx];
      const int partNode = localNodePart[idxNodeGlobal]; 
      const int idxNode  = localNodeIndices[idxNodeGlobal];

      Node *node_ptr = nodes[partNode] + idxNode;
      
      nodeNeighNodes[idx] = node_ptr;
   }


   for (int idx = idxElemStart; idx < idxElemEnd; idx++) {
      const int idxElemGlobal = nodeNeighElemsIndices[idx];
      const int partElem = localElemPart[idxElemGlobal];
      const int idxElem  = localElemIndices[idxElemGlobal];

      Element *elem_ptr = elems[partElem] + idxElem;

      nodeNeighElems[idx] = elem_ptr;
   }      

   // Compute the cosine and sine of the mesh rotation angle. We use float here because memory efficiency is more important than precision in this case.
   assert(wind_rotate_angle[i] >= -pi && wind_rotate_angle[i] <= pi && "wind_rotate_angle out of range");
   sincosf((float)wind_rotate_angle[i], &node.sRot, &node.cRot);

   // Fills the node coordinates. If ics == 1, the coordinates are in the local frame. If ics == 2, the coordinates are in the global frame.
   if constexpr(ics == 1) {
      node.p.x   = node_coords[i * 3];
      node.p.y   = node_coords[i * 3 + 1];
      node.gz = 0.0;

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
   if      constexpr(nchi ==  0) node.drag     = drag[i];
   else if constexpr(nchi ==  1) node.rough    = roughness[i];
   else if constexpr(nchi == -1) node.rmanning = rmanning[i];
   else assert(false && "populateNodes: selected nchi not implemented.");
   
   if constexpr(ivcor == 1) 
   {  
      node.kb = kb[i];
      for (int k = 0; k < nvrt; k++) node.localSigma[k] = localSigma[k + i*nvrt];
   }

   for (int k = 0; k < nvrt; k++) node.rho[k] = rho0;
   for (int k = 0; k < nvrt; k++) node.z[k] = 0.0;
   for (int k = 0; k < nvrt; k++) node.dfv[k] = dfv0;   // XXXX 

   
   // Links the node to the boundary information
   if (isbnd[i] >= 0) {
     
      // Get the index of the node boundary
      const int idxNodeBoundary = atomicAdd(sharedIdxNodeBoundary, 1) + 1;

      node.nodeBnd = &(localNodeBoundaryMem[idxNodeBoundary]);
      node.nodeBnd->node = &node;

      node.nodeBnd->sideBnd1 = nullptr;
      node.nodeBnd->sideBnd2 = nullptr;

      if (isbnd[i]>=0) {

         node.nodeBnd->segment = &(localSegmentMem[isbnd[i]]);
         
         node.nodeBnd->segment->ifltype = 0 ; //(isbnd[i]==0)?3:0; // XXXXX KKKKK
         node.nodeBnd->segment->iettype = 0; // XXXXXXXXX

         // Initialize the length and area of the segment to nan to indicate that they are not yet computed
         node.nodeBnd->segment->clen  = nan("1");
         node.nodeBnd->segment->carea = nan("1");

      }
   } else {

      node.nodeBnd = nullptr; // For internal nodes

      if (isbnd[i] == landSegment) {
         node.nodeBnd = &landBoundaryNode;
         node.nodeBnd->segment = &landBoundarySegment;
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

   side.idx = i;  
   
   // Parameters used to compute a linear approximation of the Coriolis coefficient around the reference latitude sfea0
   const double fc   = 2.0*omega_e*sin(sfea0);
   const double beta = 2.0*omega_e*cos(sfea0);

   assert(sfea0 >= -pi/2.0 && sfea0 <= pi/2.0 && "GRID: sfea0 out of range");

   // finding the pointer to the nodes that form the side
   const int idx_n1 = isidenode[globalSideIndices[i]*2];
   const int idx_n2 = isidenode[globalSideIndices[i]*2 + 1];

   const int locidx_n1 = localNodeIndices[idx_n1];
   const int locidx_n2 = localNodeIndices[idx_n2];
   const int part_n1 =   localNodePart[idx_n1];
   const int part_n2 =   localNodePart[idx_n2];

   side.n1 = (Node *)&(mesh->nodes[part_n1][locidx_n1]);
   side.n2 = (Node *)&(mesh->nodes[part_n2][locidx_n2]);


   // Counts the number of sides connected to each node
   side.n1->numSides++;
   side.n2->numSides++;

   // Side boundary type
   NodeBoundary *n1bnd = side.n1->nodeBnd;
   NodeBoundary *n2bnd = side.n2->nodeBnd;

   if (n1bnd != nullptr && n2bnd != nullptr) 
   // If both nodes are on the boundary, then the side may be on the boundary (Notice that there are cases where the side is not on the boundary, but both nodes are on the boundary)
   {
      bool isLand = false;
      int isbs = findSideBoundaryType(side, iond, nond, ilnd, nlnd, isLand);
      
      if (isbs > -1) {
         // assert(isLand || isbs>=0);
         const int idxSideBoundary = atomicAdd(sharedIdxSideBoundary, 1) + 1;

         side.sideBnd = &(localSideBoundaryMem[idxSideBoundary]);

         if   (isLand) side.sideBnd->segment = &landBoundarySegment;
         else          side.sideBnd->segment = &(localSegmentMem[isbs]);

         assert (isLand == side.sideBnd->segment->isLand);

         // Links each node boundary to the side boundary using the first available side boundary pointer
         if (n1bnd->sideBnd1 == nullptr) n1bnd->sideBnd1 = side.sideBnd;
         else                            n1bnd->sideBnd2 = side.sideBnd;

         if (n2bnd->sideBnd1 == nullptr) n2bnd->sideBnd1 = side.sideBnd;
         else                            n2bnd->sideBnd2 = side.sideBnd;
         
         side.sideBnd->side = &side;

      }
      else {
         side.sideBnd = nullptr;
      }
   }
   else 
   {
      side.sideBnd = nullptr;
   }


   // finding the pointer to the elements that share the side
   const int idx_e1 = isdel[globalSideIndices[i]*2];
   const int idx_e2 = isdel[globalSideIndices[i]*2 + 1];

   if (idx_e1 > -1) {
      const int locidx_e1 = localElemIndices[idx_e1];
      const int part_e1 =  localElemPart[idx_e1];
      side.elems[0] = elems[part_e1] + locidx_e1;

   }
   else {
      side.elems[0] = nullptr;
   }

   if (idx_e2 > -1) {
      const int locidx_e2 = localElemIndices[idx_e2];
      const int part_e2 =  localElemPart[idx_e2];
      side.elems[1] = elems[part_e2] + locidx_e2;
   }
   else {
      side.elems[1] = nullptr;
   }
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

   // Normal, unitary vector to the side
   if constexpr(ics == 1) {
      side.n = { -(side.n1->p.y - side.n2->p.y), (side.n1->p.x - side.n2->p.x)};
      side.n = normalize(side.n);
   }
   else {
      // Unit side vector in 3D cartesian coordinates
      Vector3 nSide = { side.n1->p.x - side.n2->p.x, side.n1->p.y - side.n2->p.y, side.n1->gz - side.n2->gz };
      nSide = normalize(nSide);

      // It cross product with the unit local z-axis vector (nUp) is the
      // normal vector we are looking for.
      Vector3 n3D = normalize(cross(nSide, side.frame.nUp()));  

      // The 2D components can be obtained by projecting the 3D normal vector to the plane defined by the unit vectors nZonal and nMeridional
      const Vector2 n2D = {dot(n3D, side.frame.nZonal()), dot(n3D, side.frame.nMeridional())};
      side.n = normalize(n2D);      

   }

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
   const Vector2 vLen = side.n2->p - side.n1->p;
   side.len = sqrt(dot(vLen, vLen));


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
         // a_1 = [(x2 - x)(y1 - y) - (x1 - x)(y2 - y)]/(2*area)
         // da1/dx = (y1 - y2)/(2*area)

         // Get the horizontal coordinates of the nodes id2 and id3 in the local frame
         const Point2 p2 = elem.getLocalHorizontalNodePosition(id2);
         const Point2 p3 = elem.getLocalHorizontalNodePosition(id3);

         const double dPhi_j_dx = (p2.y - p3.y)/(2*area); 
         const double dPhi_j_dy = (p3.x - p2.x)/(2*area); 

         elem.gradPhi[j] = {dPhi_j_dx, dPhi_j_dy};

      }
      else  { // Quad; evaluate at centroid

         // Get the horizontal coordinates of the nodes id2 and id4 in the local frame
         const Point2 p2 = elem.getLocalHorizontalNodePosition(id2);
         const Point2 p4 = elem.getLocalHorizontalNodePosition(id4);

         const double dPhi_j_dx = (p2.y - p4.y)/(2*area); 
         const double dPhi_j_dy = (p4.x - p2.x)/(2*area); 

         elem.gradPhi[j] = {dPhi_j_dx, dPhi_j_dy};
      } // if Tri

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

      const int idx2 = iseg*mnond*nbfr + i + jfr*mnond; //XXXXX CHECK
      nodeBnd->efa[jfr] = efa[idx2];
      nodeBnd->emo[jfr] = emo[idx2];
      nodeBnd->umo[jfr] = {umo[idx2], vmo[idx2]};
      nodeBnd->ufa[jfr] = {ufa[idx2], vfa[idx2]};

   }
        

}

__device__ void Mesh::populateBoundarySideLinks(void) {

   SIDE_BND(sideBnd);
   const int idx = idx_side_bnd;    

   Side *side = sideBnd->side;

   // We assume that the leftSide is the side s that has s->n2 == side->n1
   // and the rightSide is the side s that has s->n1 == side->n2
   for (int i = 0; i < mesh->nbs; i++) {

      Side *s = &mesh->localSideData[i];

      if (s == side) continue;  // Do not link the side to itself
   
      if (s->n1 == side->n2) sideBnd->leftSide = s;
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

      Element *elem = node1->elements[j];



      // Compute the average gravity
      const double grav3 = elem->mean(&Node::grav2);

      // Look at the documentation for the description of the terms in the following integrals

      // Integral I_1 
      // \bar{\breve{H}} 
      const double hhatb = elem->mean(&Side::hhat);        

      assert((hhatb >= 0.0 || ihhat != 1) && "Invalid hhat");

      const int id1 = elem->nodeLocalIdx(node1); 
      const double area = elem->area;

      if (elem->isTri()) {

         // To clarify the notation in this section: 
         //   i stands for the node under consideration (node1)
         //   j (node2) and k (node3) are the other nodes of the triangle.
         const int id2 = elem->nxq(1 + id1);
         const int id3 = elem->nxq(2 + id1);
         Node *node2 = elem->nodes[id2];
         Node *node3 = elem->nodes[id3];

         Point2 x1, x2, x3;
         elem->getLocalHorizontalNodePositionsTri(x1, x2, x3, id1);


         // integral of phi_i*phi_j
         const double ipp_ii = (1.0/6.0 )*area;
         const double ipp_ij = (1.0/12.0)*area;
         const double ipp_ik = (1.0/12.0)*area;

         // integral of [(grad phi_i) * (grad phi_j)]
         const double igpgp_ii = ((1.0/4.0)/area)*dot(x3-x2, x3-x2);
         const double igpgp_ij = ((1.0/4.0)/area)*dot(x3-x2, x1-x3);
         const double igpgp_ik = - igpgp_ii - igpgp_ij; // This comes from the fact that if a,b,c are vectors forming a triangle (in order), then a+b+c = 0, hence: -a*c = a*(a+b+c-c) = a*a + a*b.

         // Computes the different components of integral I_1 of each node with respect to node i.
         const double I1_ii = ipp_ii + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ii;
         const double I1_ij = ipp_ij + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ij;
         const double I1_ik = ipp_ik + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ik;

         node1->addToSparseMatrix(node1, I1_ii);
         node1->addToSparseMatrix(node2, I1_ij); 
         node1->addToSparseMatrix(node3, I1_ik);

      }
      else { // Quad
         // XXX IMPLEMENT THIS
         assert(false && "Quad not implemented in setupCoefficientMatrix");

         // // To clarify the notation in this section: 
         // //   i stands for the node under consideration (node1)
         // //   j (node2), k (node3), l (node4) are the other nodes of the quad.
         // const int id2 = elem->nxq(1 + id1);
         // const int id3 = elem->nxq(2 + id1);
         // const int id4 = elem->nxq(3 + id1);
         // Node *node2 = elem->nodes[id2];
         // Node *node3 = elem->nodes[id3];
         // Node *node4 = elem->nodes[id4];         

         // // integrals of phi_i*phi_j, phi_i*phi_k, phi_i*phi_l
         // const double ipp_ii = elem->intPhiPhi(id1, id1);
         // const double ipp_ij = elem->intPhiPhi(id1, id2);
         // const double ipp_ik = elem->intPhiPhi(id1, id3);
         // const double ipp_il = elem->intPhiPhi(id1, id4);

         // // integrals of (grad phi_i) * (grad phi_j), (grad phi_i) * (grad phi_k), (grad phi_i) * (grad phi_l)
         // const double igpgp_ii = elem->intGradPhiGradPhi(id1, id1);
         // const double igpgp_ij = elem->intGradPhiGradPhi(id1, id2);
         // const double igpgp_ik = elem->intGradPhiGradPhi(id1, id3);
         // const double igpgp_il = elem->intGradPhiGradPhi(id1, id4);

         // // Computes the different components of integral I_1 of each node with respect to node i.
         // const double I1_ii = ipp_ii + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ii;
         // const double I1_ij = ipp_ij + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ij;
         // const double I1_ik = ipp_ik + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ik;
         // const double I1_il = ipp_il + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_il;

         // node1->addToSparseMatrix(node1, I1_ii);
         // node1->addToSparseMatrix(node2, I1_ij); 
         // node1->addToSparseMatrix(node3, I1_ik);
         // node1->addToSparseMatrix(node4, I1_il);
  
      } 

 
      // Integral I_4
      // I_4 = \int_\Omega \left[ \phi_i\eta^n + \theta\Delta t \nabla\phi_i\cdot\pmb{E} + (1-\theta)\Delta t \nabla\phi_i \cdot\pmb{U}^n \right]d\Omega
      // Where
      // U is pmb{U}^{n+1} = \pmb{E} - g\theta\breve{H}\Delta t\nabla\eta^{n+1}

      // Since \nabla\phi_i is constant over the element, for the integral of U we can compute the mean value
      // of U over the element.
      Vector2 U = {0.0, 0.0};
      for (int m=0; m < elem->i34; m++) {      
         const Side *side = elem->sides[m];

         if constexpr(ics == 1) U += side->U; 
         else                   U += side->frame.projHorVec(side->U, elem->frame);  // In global coordinates, we must project the side velocity to the element frame

      }
      U = U/elem->i34;


      const Vector2 gradPhi_1 = elem->gradPhi[id1];


      // Integral of \nabla\phi_i \cdot \pmb{U}^n
      const double I_4_1 = dot(gradPhi_1, U);


      // Integral of \nabla\phi_i \cdot \pmb{E}
      const double I_4_2 = dot(gradPhi_1, elem->ghat1);

      node1->qel += dt*area*((1 - thetai)*I_4_1 + thetai*I_4_2);      

      // Additional terms 
      double vnorm = 0.0;
      if constexpr(imm == 2)  { // Pre-compute vnorm on the element
         Vector3 ubed = {0.0, 0.0, 0.0};
         // update_bdef(time, elem->center, elem->depth, ubed); XXX
         assert(false && "computeVerticalVelocity: not implemented for imm==2");

         Vector3 gradP = {0.0, 0.0, 1.0};
         for (int m=0; m < elem->i34; m++) {
           const Vector2 gradP_2D = elem->nodes[m]->depth*elem->gradPhi[m];
           gradP += toVector3(gradP_2D, 0.0);
         }   
         vnorm = dot(ubed, gradP)/sqrt(dot(gradP, gradP));
      } // imm==2

      if (elem->isTri()) {
         for (int l=0; l<3; l++) {

            Node *node = elem->nodes[l];

            double fac = (node == node1) ? 2.0 : 1.0;
   
            if      constexpr(imm == 2) node1->qel += (1.0/12.0)*elem->area*fac*(node->eta2 + dt*vnorm);
            else if constexpr(imm == 1) node1->qel += (1.0/12.0)*elem->area*fac*(node->eta2 + (node->bdef2 - node->bdef1)); 
            else                        node1->qel += (1.0/12.0)*elem->area*fac*(node->eta2); 
           
         }


         if (!elem->isDry) {
           const Vector2 gradEta = elem->grad2(&Node::eta2);

           node1->qel -= elem->area*grav3*dt*dt*thetai*(1.0 - thetai)*hhatb*dot(elem->gradPhi[id1], gradEta);

         }

      } // isTri
      else { // isQuad
         assert(false && "Quad not implemented in setupCoefficientMatrix");
         // for (int l=0; l<4; l++) {
         //   const Node *node = elem->nodes[l];

         //   if (imm==2) tmp2 = dt*vnorm; // b_t*dt
         //   else        tmp2 = node->bdef2 - node->bdef1;

         //   node1->qel += (node->eta2 + tmp2)*quad_int(1,ie,id,l) - 
         //                  node->eta2*grav3*thetai*(1.0 - thetai)*dt*dt*hhatb*quad_int(2,ie,id,l)
         // }

      } // isQuad


     node1->qel += (dt/elem->i34)*elem->vsource;
      

   } // for all elements connected to the node
   

   // Boundary integrals I_{2,3,5,6}; they all vanish at land boundaries
   // I_2 and I_6 are not needed if essential b.c. are enforced by elminating rows and columns
   NodeBoundary *nodeBnd = node1->nodeBnd;
   if (nodeBnd != nullptr && nodeBnd != &(mesh->landBoundaryNode))  { 
      // If it is an open boundary node...

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
   __global__ void setupCGSolvers(double* data_d, double* b_d, 
                                  int* indices, int* indptr) {

      NODE(node);
      
      const int idx = node->idx;
      
      // Computes the index in the nzz array      
      const int nzzIdx = (node->nodes - mesh->nodeNeighNodes) + idx;


      indptr[idx] = nzzIdx;
      indices[nzzIdx] = idx;

      data_d[nzzIdx]    = node->coeffDiag;

      b_d[idx] = node->qel;

      for (int k = 0; k < node->numNodes; ++k) {
         const Node* nb = node->nodes[k];

         const int j = nb->idx;

         data_d[nzzIdx + k + 1]    = node->coeffMatrix[k];

         indices[nzzIdx + k + 1] = j;
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
      //  Initialize (xt,yt,zt), nnel and vel.
      //  For ics=2, the coordinates are all in the side's local frame
      //  Caution! nnel must be initialized inside this loop as it is updated inside.

         const int iadvf = side->iadv();  // Advection flag   
         assert(nadv == 1 && "Advection flag is only implemented for nadv = 1");

         // If on land or open boundary...
         const bool ifl_bnd =  (side->sideBnd != nullptr);
         
         const double vmag = length(side->u2[k]);

         if (vmag <= velmin_btrack) {
            // No activity 

            side->btu[k] = side->u2[k];  // I don't think this is necessary because it has already been set to this value XXXX
// addDebugPlot(side->center, side->center, 0x00000F, '+');     

            
         } else {  // do btrack

            // Compute number of sub-divisions based on local gradients of the velocity field

            // XXX IS THIS NECESSARY?
            // The gradient at a side is the average of the gradients at the elements sharing the side.
            double suma = 0.0;
            int icount = 0;
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
            } // for ii

            assert(icount > 0 && "Impossible 77");

            // '4' is somewhat arbitrary
            const int ndelt = clamp(int(suma/icount)*4, ndelt_min, ndelt_max); 
            const double dtbk = dt/ndelt; // Target btrack step; may be smaller sometimes
            Vector3 u_int = {0.0, 0.0, 0.0};

            bool res =side->btrack(k, dt, dtbk, u_int);

            side->btu[k] = toVector2(u_int);
               
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

      NODE_BND(nodeBnd);

      Node *node = nodeBnd->node;

      node->ETP = 0.0;

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
   //      if(rat<-small1.or.rat>1.d0+small1) then
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
      if (nramp_elev == 1) eic = nodeBnd->etaic;
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
         node->coeffDiag = 1.0;

         for (int j = 0; j < node->numNodes; j++) {
               node->coeffMatrix[j] = 0.0;
         }
         
         // Set solution and RHS to BC value
         node->mX  = nodeBnd->elbc;
         node->qel = nodeBnd->elbc;

      } else {
         // Not an essential BC node. Still need to modify RHS and columns that belong to neighbors with essential BCs

         for (int j = 0; j < node->numNodes; j++) {
               Node* nbNode = node->nodes[j];
               const NodeBoundary* nbNodeBnd = nbNode->nodeBnd;
               if (nbNodeBnd && nbNodeBnd->isEssentialBC) {
                  // Verify neighbor's BC value
                  assert(nbNodeBnd->elbc >= -9998.0 && "setEssentialBC: wrong b.c. at node");

                  // Subtract contribution to RHS
                  node->qel -= node->coeffMatrix[j]*nbNodeBnd->elbc;
                  node->coeffMatrix[j] = 0.0;  
               }
         }
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

      // Compute the drag coefficients chi2 and chi (\tilde\chi in notes) // XXX change names, chi->chitilde, chi2->chi
      // Velocity magnitude at the bottom level
      // XXXX Check if we have added -dt*(bcc(1,k,n1)+bcc(1,k,n2))/2.d0-dt*swild2(k,1) to u2.
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
      if (iveg != 0 && side->is3D() && veg_h > 0.0) {
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
               node->u2[k]     = toVector2(utmp);
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

      for (int k = node->kb; k < nvrt; k++) node->w2[k] = node->w2[k] / node->tmp.d;  // node->tmp.d contains the total area of the elements sharing the node

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

                  Vector2 u = s2u2 + s3u2 - s1u2; // u interpolated linearly at node j.   
                  u = clampVelocity(u);
                  atomicAdd(&(node->utmp[k].x), u.x);
                  atomicAdd(&(node->utmp[k].y), u.y);

               }

            }
         }
         else { // For quads
            assert(false && "computeNodalHorizontalVelocities1: not implemented for quads");

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
         Vector3 elem_top_n = -elem->n(l + 1, topArea3D);  // Get the normal vector and the area of the top face (as a 3D surface)

         assert(topArea3D > 0.0 && "computeVerticalVelocity: topArea3D <= 0");
         assert(elem_top_n.z < 0.0 && "computeVerticalVelocity: elem_top_n.z >= 0");

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
         const double topHorizontalFlux = dot(u2meanTop, toVector2(elem_top_n))*topArea3D; // Part of the flux related to the horizontal velocity

         // From the flux balance, (topVerticalFlux + topHorizontalFlux) + bottomFlux = 0, so:
         const double topVerticalFlux = -sideFlux - topHorizontalFlux - bottomFlux;  

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
                           Node **nodesPtr, Element **elemsPtr, Side **sidesPtr,
                           Node *nodes, Element *elems, Side *sides,
                           Node **partNodeNeighNodesMem, Element **partNodeNeighElemsMem,
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
                   nodesPtr, elemsPtr, sidesPtr, 
                   nodes, elems, sides,
                   partNodeNeighNodesMem, partNodeNeighElemsMem, 
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
   int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
   int *localNodePart,     int *localElemPart,     int *localSidePart,
   double *drag, double *roughness, double *rmanning, int *kb, double *localSigma) {

   mesh->populateNodes(i34, isbnd, sharedIdxNodeBoundary, 
                       node_coords, node_coords_ll, use_ll, wind_rotate_angle, 
                       nodeNeighNodesStart, nodeNeighNodesIndices,
                       nodeNeighElemsStart, nodeNeighElemsIndices,
                       localNodeIndices, localElemIndices, localSideIndices,
                       localNodePart, localElemPart, localSidePart,
                       drag, roughness, rmanning, kb, localSigma);
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
   mesh->populateBoundarySideLinks();

}





extern "C"  __global__ void destroyMeshCPU(void) {
   delete mesh;
}


extern "C" __global__ void test(double x, double y, double z) {
   mesh->test(x, y, z);
}


__global__ void getStructSizes(int* sizes)  {
   // Get the size of each struct. This is used to allocate memory for the structs in the host code.
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
      res[idx_node*size + idx_node] = node->coeffDiag;
      for (int j = 0; j < node->numNodes; j++) {
         Node* nbNode = node->nodes[j];
         if (node->coeffMatrix[j] != 0.0) res[nbNode->idx*size + idx_node] =  node->coeffMatrix[j];
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

#ifdef OLD_PCG_SOLVER
__global__ void computeMatrixProductAndDotProducts() {
   // Computes A*p, r'Mr and p'Ap

   // 1) shared buffer: first half for rMr, second for pAp
   __shared__ double sdata[2*BLOCK_SIZE];
   double* s_rMr = sdata;
   double* s_pAp = sdata + BLOCK_SIZE;

   const int tid = threadIdx.x;

   double l_pAp = 0.0;

   s_rMr[tid] = 0.0;
   s_pAp[tid] = 0.0;

   // 2) get your node pointer
   // NODE(node);   
   const int idx_node = threadIdx.x + blockIdx.x * blockDim.x; 
   Node *node = &(mesh->localNodeData[idx_node]); 
      if (idx_node < mesh->np) {

      const double diag = node->coeffDiag;
      const double r    = node->mR;
      const double p    = node->mP;  

      // apply preconditioner on the fly:
      // const double Mr   = r/diag;   // M is the preconditioner (inverse of the diagonal)

      // 4) compute Ap = A.p
      double Ap = node->coeffDiag * node->mP;
      for (int j = 0; j < node->numNodes; ++j) {
         const Node *nbNode = node->nodes[j];
         Ap += node->coeffMatrix[j] * nbNode->mP;

      }

      node->mAp = Ap;
      

      // 5) dot-products individual products (r'Mr, p'Ap) (Will be reduced later)
      // const double l_rMr = r*Mr;  // r'Mr
      l_pAp = p*Ap;  // p'Ap


   }

   // 6) write into shared mem
   // s_rMr[tid] = l_rMr;
   s_pAp[tid] = l_pAp;

   __syncthreads();

   // 7) block-wide tree reduction
   for (int stride = BLOCK_SIZE/2; stride > 0; stride >>= 1) {
      if (tid < stride)  {
         s_pAp[tid] += s_pAp[tid + stride];
      }
      __syncthreads();
   }

   // 8) thread 0 atomically accumulates into the globals
   // XXX: This is slower than it needs to be, but it is simple and works
   if (tid == 0) {
      atomicAdd(pAp, s_pAp[0]);
   }


}


__global__ void updateSolutionAndDotProducts(void) {
   extern __shared__ double s_rMr[BLOCK_SIZE];       

   const int tid = threadIdx.x;

   s_rMr[tid] = 0.0;

   // get node & bounds-check:
   NODE(node);

   
   // Compute alpha
   const double alpha = (*rMr_old)/(*pAp); 


   assert(*pAp != 0.0 && "updateSolutionAndDotProducts: pAp = 0.0");


   // load scalars
   const double diag = node->coeffDiag;
   const double p    = node->mP;
   const double Ap   = node->mAp;
   const double x    = node->mX;
   const double r    = node->mR;
   
   assert(diag != 0.0 && "updateSolutionAndDotProducts: diag = 0.0");

   // 1) update solution & residual
   const double x_new = x + alpha*p;
   const double r_new = r - alpha*Ap;
   node->mX = x_new;
   node->mR = r_new;

   // 2) apply Jacobi preconditioner on the fly
   const double Mr_new = r_new/diag;


   // 3) local contributions for the two dot-products
   const double l_rMr = r_new*Mr_new;  // r'M^-1r

   // 4) write into shared-mem
   s_rMr[tid] = l_rMr;

   __syncthreads();


   // 5) block-wide tree reduction (O(log BLOCK_SIZE))
   for (int stride = BLOCK_SIZE/2; stride > 0; stride >>= 1) {
      if (tid < stride) s_rMr[tid] += s_rMr[tid + stride];
      __syncthreads();
   }

   // 6) thread 0 atomically folds into the global accumulators
   if (tid == 0) atomicAdd(rMr_new, s_rMr[0]);


}


__global__ void updateSearchDirection() {

    NODE(node);

    // recompute z on-the-fly
    const double z = node->mR / node->coeffDiag;
    const double p = node->mP;

    const double beta = (*rMr_new)/(*rMr_old);


    // p =z + beta*p
    node->mP = z + beta*p;

}

__global__ void initializeScalars(double *scalars) {

   // These point to a cupy array in the host code
   rMr_old = scalars + 0;
   rMr_new = scalars + 1;
   pAp     = scalars + 2;
}

__global__ void initializePreconditionedResidual(void) {

    __shared__ double s_rMr[BLOCK_SIZE];  

    s_rMr[threadIdx.x] = 0.0;

    NODE_ALL(node, valid);
    const int tid = threadIdx.x;
    // For the first iteration:
    // 1. Initial solution x = 0 (already set)
    // 2. Initial residual r = b (stored in qel)
    // 3. Initial search direction p = M^-1 * r = M^-1 * b

// XXX we should use the guess solution here    

    double l_rMr = 0.0;
    if (valid) {
      // Set initial residual r = b (from qel)
      node->mR = node->qel;

      // Compute z = M^-1 * r (preconditioned residual)
      const double z = node->mR / node->coeffDiag;

      // Set initial search direction p = z
      node->mP = z;

      // Initialize other variables needed for CG
      node->mX = 0.0;  // Initial solution
      node->mAp = 0.0; // Will be computed in first iteration


      // Computes the first r'Mr (reduction)
      l_rMr = node->mR*node->mR/node->coeffDiag;  // r'M^-1r
    }
    // 4) write into shared-mem
    s_rMr[tid] = l_rMr;

    __syncthreads();


    // 5) block-wide tree reduction (O(log BLOCK_SIZE))
    for (int stride = BLOCK_SIZE/2; stride > 0; stride >>= 1) {
        if (tid < stride) s_rMr[tid] += s_rMr[tid + stride];
        __syncthreads();
    }

    // 6) thread 0 atomically folds into the global accumulators
    if (tid == 0) atomicAdd(rMr_new, s_rMr[0]);        

    *rMr_old = 0.0;
    *pAp = 0.0;
}
#endif

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