 
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

#ifdef DEBUG
  #define assertFmt(cond, fmt, ...)                                      \
    do {                                                                   \
      if (!(cond)) {                                                       \
        printf("ASSERT FAILED: " fmt "\n", ##__VA_ARGS__);                \
        assert(false && "ASSERT FAILED"); \
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
   __device__ double& operator[](int k) { return (double &)*((double *)nullptr); }
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
   __device__ __forceinline__ void getLocalHorizontalNodePositionsQuad(Point2 &x0, Point2 &x1, Point2 &x2, Point2 &x3, int idx = 0) const;


   // (eta,xi) coordinates of the Gauss-Legendre quadrature points (idxNode = 0,1,2,3)
   __device__ __forceinline__ Point2 GaussPoint(int idxNode) const {
       constexpr double g = 0.5773502691896257; // 1/sqrt(3)
       assert(idxNode >= 0 && idxNode < 4 && "GaussPoint: idxNode out of range");
       const double eta = (idxNode & 1) ? +g : -g;
       const double xi  = (idxNode & 2) ? +g : -g;
       return {eta, xi};
   }  



   // Shape function Phi at a given point p = (eta,xi), idxNode = node index 0,1,2,3
   __device__ __forceinline__ double Phi(int idxNode, Point2 p) const {
      return 0.25*(1.0 + p.x*eta_n(idxNode))*(1.0 + p.y*xi_n(idxNode));
   }

   // Derivative of shape function Phi_idx/dxi at a given point p = (eta,xi)
   __device__ __forceinline__ double dPhi_dxi(int idxNode, Point2 p) const {
      return 0.25*xi_n(idxNode)*(1.0 + p.x*eta_n(idxNode));
   }   

   // Derivative of shape function Phi_idx/deta at a given point p = (eta,xi) 
   __device__ __forceinline__ double dPhi_deta(int idxNode, Point2 p) const {
      return 0.25*(1.0 + p.y*xi_n(idxNode))*eta_n(idxNode);
   }      
   
   __inline__ __device__ double detJ_eta_quad(Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const;

   __inline__ __device__ double detJ_xi_quad(Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const;

   __inline__ __device__ double jacobianDetQuad(Point2 p, Point2 const &x0, Point2 const &x1, Point2 const &x2, Point2 const &x3) const 
         { return area/4.0 + detJ_eta_quad(x0, x1, x2, x3)*p.x + detJ_xi_quad(x0, x1, x2, x3)*p.y; }

   // Check if the element is dry using the more up-to-date elevation.
   __inline__ __device__ bool isDry2(void) const;

   __inline__ __device__ bool isEdgeBoundary(Side *side) const;

   __device__ bool inverseBilinearQuad(Point2 p, double &xi, double &eta, double shp[4]) const;

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
   double depth;

   Frame<ics> frame;

   LocalSigma<ivcor> localSigma;  

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
   double bdef1 = 0.0; // Bottom deformation at previous timestep
   double bdef2 = 0.0; // Bottom deformation at current timestep

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


