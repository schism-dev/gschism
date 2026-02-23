


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


__device__ __forceinline__ double clampVelocity(double x, double min, double max);

// VECTOR/POINT operations ---------------------------------------------------------------------------

using Vector2 = double2; 
using Vector3 = double3;
using Point2  = double2;
using Point3  = double3;

__device__ __forceinline__ double2 operator+(const double2& a, const double2& b) {
    return make_double2(a.x + b.x, a.y + b.y);
}
using Vector2 = double2; 
using Vector3 = double3;
__device__ __forceinline__ double2 operator+=(double2& a, const double2& b) {
    a.x += b.x;
    a.y += b.y;
    return a;
}

__device__ __forceinline__ double2 operator-(const double2& a, const double2& b) {
    return make_double2(a.x - b.x, a.y - b.y);
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

__device__ __forceinline__ double norm2(const double3& v) {
    return v.x*v.x + v.y*v.y + v.z*v.z;
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
}


__forceinline__ __device__ void computeLonLat(Point3 const &p, double &rlon, double &rlat) {
// Given a 3D global point (may not be on surface of earth), find its longitude and latitude in radians

   const double radius2 = dot(p, p);
   assert(radius2 != 0.0     && "computeLonLat: radius = 0");
   assert(p.z*p.z <= radius2 && "computeLonLat: p.z <= radius");

   rlon = atan2(p.y, p.x); // (-pi,pi]

   if (abs(rearth_pole - rearth_eq) < 1.0e-2) { // for backward compatibility
      rlat = asin(p.z/sqrt(radius2));
   } else {
      rlat = asin(p.z/rearth_pole);
   }
  
}

__inline__ __device__ void Error(const char* msg) {
   printf("Error: %s\n", msg);
   assert(false);
}

struct Tracers {
    double tracer[N_TRACERS];
};



struct Frame {
   double _t[9]; // t stands for tensor.

   // Access the tensor in Fortran order
   __device__ __forceinline__ double& t(int i, int j) { return _t[i + 3*j]; }

   __device__ __forceinline__ double t(int i, int j) const { return _t[i + 3*j]; }

   __device__ void update(double ar1, double ar2) { // XXX use lat, lon instead of ar1, ar2?

      double sin_ar1, cos_ar1, sin_ar2, cos_ar2;
      sincos(ar1, &sin_ar1, &cos_ar1); // These functions are faster than sin(ar1) and cos(ar1).
      sincos(ar2, &sin_ar2, &cos_ar2);

      t(0,0) = -sin_ar1;
      t(1,0) =  cos_ar1;
      t(2,0) =  0.0;
      t(0,1) = -cos_ar1*sin_ar2;
      t(1,1) = -sin_ar1*sin_ar2;

      t(2,1) = (rearth_pole/rearth_eq)*cos_ar2;

      const double ar4 = sqrt(t(0,1)*t(0,1) + t(1,1)*t(1,1) + t(2,1)*t(2,1));
      assert (ar4!=0.0 && "GRID: ar4=0");

      const double ar4_inv = 1.0/ar4;
      t(0,1) = t(0,1)*ar4_inv;
      t(1,1) = t(1,1)*ar4_inv;
      t(2,1) = t(2,1)*ar4_inv;

      // The third vector is the cross product of the first two
      crossProduct(t(0,0), t(1,0), t(2,0), 
                   t(0,1), t(1,1), t(2,1), 
                   t(0,2), t(1,2), t(2,2));
   }


   __forceinline__ __device__ Vector2 project_hvec(Vector2 const &vec0, Frame const &frameout) const {
   // Do plane rotation of a vector (i.e., assuming z-axes are same)
   // Inputs:
   //        vec0: vel. in this frame;
   //        frameout: Output frame;
   // Output: 
   //        vec1: vel. in new frame  

      auto const fo = frameout;
      Vector2 vec1;

      // Since frame tensors are orthogonal, the transpose is the inverse. That's why the indices
      // of t and fo.t are the same.
      vec1.x = vec0.x*(t(0,0)*fo.t(0,0) + t(1,0)*fo.t(1,0) + t(2,0)*fo.t(2,0)) +
               vec0.y*(t(0,1)*fo.t(0,0) + t(1,1)*fo.t(1,0) + t(2,1)*fo.t(2,0));

      vec1.y = vec0.x*(t(0,0)*fo.t(0,1) + t(1,0)*fo.t(1,1) + t(2,0)*fo.t(2,1)) +
               vec0.y*(t(0,1)*fo.t(0,1) + t(1,1)*fo.t(1,1) + t(2,1)*fo.t(2,1));

      return vec1;
   }

   __forceinline__ __device__ void G2L(Vector3 const &vi, Vector3 &vo) const {
      // Transform a vector from global to local frame
      vo.x = vi.x*t(0,0) + vi.y*t(1,0) + vi.z*t(2,0);
      vo.y = vi.x*t(0,1) + vi.y*t(1,1) + vi.z*t(2,1);
      vo.z = vi.x*t(0,2) + vi.y*t(1,2) + vi.z*t(2,2);
   }

   __forceinline__ __device__ void L2G(Vector3 const &vi, Vector3 &vo) const {
      // Transform a vector from local to global frame
      vo.x = vi.x*t(0,0) + vi.y*t(0,1) + vi.z*t(0,2);
      vo.y = vi.x*t(1,0) + vi.y*t(1,1) + vi.z*t(1,2);
      vo.z = vi.x*t(2,0) + vi.y*t(2,1) + vi.z*t(2,2);
   }   

   __forceinline__ __device__ void G2L(Point3 const &pi, Point3 const &origin0, Point3 &po) const {
   // Tranforms a point from global to local frame     
      G2L(pi - origin0, po);
   }

   __forceinline__ __device__ void L2G(Point3 const &pi, Point3 const &origin0, Point3 &po) const {
   // Transforms a point from local to global frame      
      L2G(pi, po);
      po = po + origin0;
   }

};


// Forward declarations
struct Side;
struct Element;
struct Node;
struct Mesh;
struct NodeBoundary;

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

   Frame frame;

   Vector2 gradPhi[4];        // Derivatives of shape functions
   Vector2 gradEta;           // Gradient of free surface elevation
   Vector2 gradPr;            // Gradient of pressure    
   Vector2 gradETP;           // Gradient of Earth Tidal Potential    XXXXX Really used?

   bool isDry;                // True if the element is dry

   Vector2 center;         // Element center coordinates
   double xlon, ylat;         // Element center lat/lon coordinates in _degrees_
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
   

   __inline__ __device__ double demeanedDensity(int k) const { return rho[k] - rho_mean[k]; }   

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

   // (eta,xi) coordinates of the Gauss-Legendre quadrature points (idx = 0,1,2,3)
   __device__ __forceinline__ Point2 GaussPoint(int idx) const {
       constexpr double g = 0.5773502691896257; // 1/sqrt(3)
       assert(idx >= 0 && idx < 4 && "GaussPoint: idx out of range");
       const double eta = (idx & 1) ? +g : -g;
       const double xi  = (idx & 2) ? +g : -g;
       return {eta, xi};
   }  

   // (eta,xi) coordinates of the Gauss-Legendre quadrature points (idx = 0,1,2,3)
   __device__ __forceinline__ double etaSign(int idx) const {
      // bit1: 0 -> -1, 1 -> +1
      assert(idx >= 0 && idx < 4 && "etaSign: idx out of range");
      return (idx & 2) ? 1.0 : -1.0;     // [-1,-1,+1,+1]
   }

   __device__ __forceinline__ double xiSign(int idx) const {
      // +1 if (bit0 XOR bit1), else -1
      assert(idx >= 0 && idx < 4 && "xiSign: idx out of range");
      return ((idx & 1) ^ (idx >> 1)) ? 1.0 : -1.0;   // [-1,+1,+1,-1]

   }

   // Shape function Phi_idx at a given point p = (eta,xi), idx = 0,1,2,3
   __device__ __forceinline__ double Phi(int idx, Point2 p) const {
      return 0.25*(1.0 + p.x*etaSign(idx))*(1.0 + p.y*xiSign(idx));
   }

   // Derivative of shape function Phi_idx/dxi at a given point p = (eta,xi) 
   __device__ __forceinline__ double dPhi_dxi(int idx, Point2 p) const {
      return 0.25*xiSign(idx)*(1.0 + p.x*etaSign(idx));
   }   

   // Derivative of shape function Phi_idx/deta at a given point p = (eta,xi) 
   __device__ __forceinline__ double dPhi_deta(int idx, Point2 p) const {
      return 0.25*(1.0 + p.y*xiSign(idx))*etaSign(idx);
   }      
   
   __inline__ __device__ double detJ_eta(void) const;

   __inline__ __device__ double detJ_xi(void) const;

   __inline__ __device__ double jacobianDet(Point2 p) const { return area/4.0 + detJ_eta()*p.x + detJ_xi()*p.y; }

   // Check if the element is dry using the more up-to-date elevation.
   __inline__ __device__ bool isDry2(void) const;

   __inline__ __device__ bool isEdgeBoundary(Side *side) const;

   __inline__ __device__ double computeSignedArea(void) const;

   __inline__ __device__ void nudgeToBarycenter(Point2 &p) const {
      // Nudges the point (x,y) to the barycenter of the element

      p = (1.0 - btrack_nudge)*p + btrack_nudge*center;
   }

   enum EFixIf { fiNo, fiOnlyInside, fiAnywhere };

   __device__ bool computeShapeCoords(Point2 &xt, double shpCoords[4], EFixIf fixNegShapeCoords = fiNo) const;
   
   __device__ void shapeCoordsToXY(double shpCoords[4], Vector2 &xt) const;

   __device__ double intPhiPhi(int i1, int i2) const;

   __device__ double intGradPhiGradPhi(int i1, int i2) const;

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


   __device__ void setFrame(void) {

      double ar1, ar2;
      assert(false && "setFrame: center should be a Vector3???");
     // computeLonLat(center, ar1, ar2); XXXX

      frame.update(ar1, ar2);        
    }

   __device__ bool findNextIntersection(Point2 x1, Point2 x2, Side *&side, Element *&nextElement, double &t, 
                                        double& tIntersect);

   __inline__ __device__ bool containsPoint(Vector2 xt) const { 
      // Check if the point (xt,yt) is inside the element
      double shpCoords[4];

      bool res = computeShapeCoords(xt, shpCoords, fiNo);

      return res;
   }

   __inline__ __device__ bool containsPoint(Vector2 &xt, bool fixCoords) const { 
      // Check if the point (xt,yt) is inside the element.
      // If fixCoords is true, and *the point is inside the element*, the point is nudged to be inside the element.
      double shpCoords[4];

      Vector2 tmp_xt = xt;
      const EFixIf fixIfCoords = fixCoords ? fiOnlyInside : fiNo;
      bool res = computeShapeCoords(tmp_xt, shpCoords, fixIfCoords);

      if (fixCoords && res) xt = tmp_xt;

      return res;
   }

   __device__ double interpolate(const double shapeCoords[4], const double zAtNodes[4]) const
   {
      // Interpolate the value at the point (xt,yt) using the shape coordinates and the values at the nodes
      // This function works for both triangular and quadrilateral elements, since shapeCoords[3] is always 0.0 for triangular elements
      return shapeCoords[0]*zAtNodes[0] + shapeCoords[1]*zAtNodes[1] + shapeCoords[2]*zAtNodes[2] + shapeCoords[3]*zAtNodes[3];
      
   }

   __device__ void interpolateVelocities(Point3 xt, Vector3& u_int) const;

   // ------------------------------------------------------------------------------------------------
   // These are functions that work on all nodes of the element.
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
   __device__ void getNodalValues(T Node::* field, T res[4]) const {
      // Get the values of the field at the nodes of the element
      res[0] = nodes[0]->*field;
      res[1] = nodes[1]->*field;
      res[2] = nodes[2]->*field;
      if (isQuad()) res[3] = nodes[3]->*field;
      else          res[3] = T();
   }

   template<typename T>
   __device__ void getNodalValues(int level, T (Node::* field)[nvrt], T res[4]) const {
      // Get the values of the field at the nodes of the element
      res[0] = (nodes[0]->*field)[level];
      res[1] = (nodes[1]->*field)[level];
      res[2] = (nodes[2]->*field)[level];
      if (isQuad()) res[3] = (nodes[3]->*field)[level];
      else          res[3] = T();
   }    

   template<typename T>
   __device__ void getNodalValues(int level, T (Node::*func)(int) const, T res[4]) const {
      // Get the values of the func at the nodes of the element
      res[0] = (nodes[0]->*func)(level);
      res[1] = (nodes[1]->*func)(level);
      res[2] = (nodes[2]->*func)(level);
      if (isQuad()) res[3] = (nodes[3]->*func)(level);
      else          res[3] = T();  // This is 0.0 or a zero vector

   }      

   template<typename T>
   __device__ void getSideValues(T Side::* field, T res[4]) const {
      // Get the values of the field at the nodes of the element
      res[0] = sides[0]->*field;
      res[1] = sides[1]->*field;
      res[2] = sides[2]->*field;
      if (isQuad()) res[3] = sides[3]->*field;
      else          res[3] = T();
   }

   template<typename T>
   __device__ void getSideValues(int level, T (Side::* field)[nvrt], T res[4]) const {
      // Get the values of the field at the nodes of the element
      res[0] = (sides[0]->*field)[level];
      res[1] = (sides[1]->*field)[level];
      res[2] = (sides[2]->*field)[level];
      if (isQuad()) res[3] = (sides[3]->*field)[level];
      else          res[3] = T();
   }    

   template<typename T>
   __device__ void getSideValues(int level, T (Side::*func)(int) const, T res[4]) const {
      // Get the values of the func at the nodes of the element
      res[0] = (sides[0]->*func)(level);
      res[1] = (sides[1]->*func)(level);
      res[2] = (sides[2]->*func)(level);
      if (isQuad()) res[3] = (sides[3]->*func)(level);
      else          res[3] = T();  // This is 0.0 or a zero vector

   }      

   template<typename T>
   __device__ T interpolate(int level, const double shapeCoords[4], T (Node::*func)(int) const) const {
      T values[4];
      getNodalValues(level, func, values);

      return shapeCoords[0]*values[0] + shapeCoords[1]*values[1] + shapeCoords[2]*values[2] + shapeCoords[3]*values[3];
   }

   template<typename T>
   __device__ T interpolate(const double shapeCoords[4], T Node::* field) const {
      T values[4];
      getNodalValues(field, values);
      
      return shapeCoords[0]*values[0] + shapeCoords[1]*values[1] + shapeCoords[2]*values[2] + shapeCoords[3]*values[3];
   }


   template<typename T>
   __device__ T interpolate(int level, const double shapeCoords[4], T (Node::* field)[nvrt]) const {
      T values[4];
      getNodalValues(level, field, values);

      return shapeCoords[0]*values[0] + shapeCoords[1]*values[1] + shapeCoords[2]*values[2] + shapeCoords[3]*values[3];
   }

   template<typename T>
   __device__ T mean(T Node::* field) const {
      T values[4];
      getNodalValues(field, values);

      assert(i34 == 4 || values[3] == T() && "mean: values[3] is not zero with i34 == 4");

      return (values[0] + values[1] + values[2] + values[3])/i34;
   }

   template<typename T>
   __device__ T mean(T Side::* field) const {
      T values[4];
      getSideValues(field, values);
      
      assert((i34 == 4) || (values[3] == T()) && "mean: values[3] is not zero with i34 == 4");

      return (values[0] + values[1] + values[2] + values[3])/i34;
   }   

   template<typename T>
   __device__ T mean(int level, T (Node::* field)[nvrt]) const {
      T values[4];
      getNodalValues(level, field, values);

      assert((i34 == 4) || (values[3] == T()) && "mean: values[3] is not zero with i34 == 4");

      return (values[0] + values[1] + values[2] + values[3])/i34;
   }

   template<typename T>
   __device__ T mean(int level, T (Side::* field)[nvrt]) const {
      T values[4];
      getSideValues(level, field, values);
      
      // assert((i34 == 4) || (values[3] == T()) && "mean: values[3] is not zero with i34 == 4");

      return (values[0] + values[1] + values[2] + values[3])/i34;
   }   

   // __device__ double mean(double (Node::*func)(int) const) const {
   //    double values[4];
   //    getNodalValues(level, func, values);
   //    return shapeCoords[0]*values[0] + shapeCoords[1]*values[1] + shapeCoords[2]*values[2] + shapeCoords[3]*values[3];
   // }
   // ------------------------------------------------------------------------------------------------

   __device__ bool globalToLocal(double xt, double yt, double &xlocal, double &ylocal) const {
      // Convert global coordinates to local coordinates
      printf("NOT IMPLEMENTED: globalToLocal\n");
      return true;
   }

   __inline__ __device__ bool localToGlobal(double xlocal, double ylocal, double &xt, double &yt) const {
      // Convert local coordinates to global coordinates
      printf("NOT IMPLEMENTED: localToGlobal\n");
      return true;
   }

   __device__ void trackGoesOutOfWater(Side *nextSide, Element *curElem, double tt, int jlev,
                                       Point3 &xt, double zint, Vector3 &u_int);
    
   __device__ bool trackThroughElementsAndLevels(Point3 x0, Point3 &xt, Vector3 &u_int, 
                                                 Element *&nextElem, int &jlev, double &dt, double time_rm,
                                                 const Side *origSide);

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

   Frame frame;

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

   __inline__ __device__ double x(double t) const;
   __inline__ __device__ double y(double t) const;
   __inline__ __device__ double lon(void) const;
   __inline__ __device__ double lat(void) const;
   __inline__ __device__ Point2 p_t(double t) const;

   __inline__ __device__ Point3 center3D(int k) const { assert(k < nvrt); return {center.x, center.y, zc[k]}; }
   __inline__ __device__ double Dz(int k)       const { assert(k < nvrt); return zc[k] - zc[k - 1]; }
   __inline__ __device__ double area(int k)     const { assert(k < nvrt - 1); return len*Dz(k + 1); }

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

   __inline__ __device__ Element* findElementTowards(Vector2 u) const {
      // From the elements that share this side, find the one that is in the direction opposite to the velocity vector
      // (unless it is a boundary side or a dry element, in which case we return the only wet element)

      assert(elems[0] != nullptr && "elems[0] is nullptr in findElementTowards");

      // Dot product of the negative velocity vector with the side normal vector
      const double dot_prod1 = -dot(u, n);

      // Vector from side center to element 0 centroid
      const Vector2 v2 = elems[0]->center - center;

      // Dot product of the vector with the side normal vector
      const double dot_prod2 = dot(v2, n);

      // If the dot product is positive, returns the element whose vector to the centroid is in the same direction as the target point
      // Otherwise, returns the other element (otherwise includes dot_prod1 == 0)

      const bool isDry0 = elems[0]==nullptr || elems[0]->isDry;
      const bool isDry1 = elems[1]==nullptr || elems[1]->isDry;

      assert(!isDry0 || !isDry1 && "Both elements are dry in findElementTowards");

      if (dot_prod1*dot_prod2 > 0.0) return !isDry0 ? elems[0] : elems[1];
      else                           return !isDry1 ? elems[1] : elems[0];

 
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

   Point2 p;
   Vector2 u2[nvrt];
   double  w2[nvrt];
   Vector2 utmp[nvrt];  // Temporary array to store velocities
   double z[nvrt];
   double lon, lat;
   double fun_lat[ntip + 1]; // XXXXX CHECK IF THIS MAKES SENSE. Also, the +1 is only to avoid erros when ntip ==1
   double depth;

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

   double bdef1; // Bottom deformation at previous timestep
   double bdef2; // Bottom deformation at current timestep

   double qel; // RHS
   double coeffDiag; // Diagonal coefficient of the matrix
   double *coeffMatrix; // Matrix coefficients for neighbors. [0..numNodes-1]

   double mP; // Search direction for conjugate gradient method
   double mR; // Residual for conjugate gradient method
   double mAp; // Product of matrix and search direction for conjugate gradient method
   double mX; // Solution for conjugate gradient method

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
      const Vector2 u_air = {u.x*cRot - u.y*sRot, u.x*sRot + u.y*cRot};

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


   __device__ Point3 p3D(int k) const { return {p.x, p.y, z[k]}; }

   __device__ double sigmaTransitionDepth(void) const;

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

   __device__ double compZCoord(int k) const;

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

   // // Get's the pointer to the node using the array of global pointers to nodes
   // __device__ const Node *getGlobalNode(int jn) {
   //    const int idx = jn / max_np; // Get the index of the partition that contains the node
   //    return nodes[idx] + idx;
   // }

   // // Get's the pointer to the element using the array of global pointers to elements
   // __device__ const Element *getGlobalElem(int je) {
   //    const int idx = je / max_ne; // Get the index of the partition that contains the element
   //    return elems[idx] + idx;
   // }

   // // Get's the pointer to the side using the array of global pointers to sides
   // __device__ const Side *getGlobalSide(int js) {
   //    const int idx = js / max_ns; // Get the index of the partition that contains the side
   //    return sides[idx] + idx;
   // }

   // // Get's the global index of local node idx
   // __device__ int getGlobalNodeIndex(int idx) {
   //    return idx + partNum * max_np;
   // }

   // // Get's the global index of local element idx
   // __device__ int getGlobalElemIndex(int idx) {
   //    return idx + partNum * max_ne;
   // }

   // // Get's the global index of local side idx
   // __device__ int getGlobalSideIndex(int idx) {
   //    return idx + partNum * max_ns;
   // }  

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
                                          double *vmo, double *vfa);                    

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
         auto isInside = localElemData[i].containsPoint({x, y});

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

   __global__ void initializeScalars(double *scalars);
 
   __global__ void initializePreconditionedResidual(void);

   __global__ void computeMatrixProductAndDotProducts(void);

   __global__ void updateSolutionAndDotProducts(void);

   __global__ void updateSearchDirection(void);
   
   __global__ void computeHorGradientAtSide(void);  


   __global__ void setupCGSolvers(double* data_d, double* b_d, 
                                  int* indices, int* indptr);



}


