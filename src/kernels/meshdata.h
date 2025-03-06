
// XXX: Pure constants that are to be set before compiling the code
constexpr int N_TRACERS = 2;

constexpr int ivcor = 2; // Vertical coordinate system
constexpr int ics = 1;   // 1: global, 2: local
constexpr int ihh = 1;   // Option to limit \hat{H} to enhance stability for large friction in shallow area

constexpr int nvrt = 54; // Number of vertical levels
constexpr int nzlev = 5; // Number of z levels (in contrast to sigma levels)
constexpr int nsiglev = nvrt - nzlev + 1; // Number of sigma levels (in contrast to z levels)

constexpr double small1 = 1.0e-10;       // Small non-negative number
constexpr double small2 = small1*100.0;  // Slightly larger than small1
constexpr double small = 0.0;            // Small positive number or 0

constexpr int internalNode = -1;
constexpr int landNode = -2;

constexpr double pi = 3.14159265358979323846;

constexpr double rearth_eq   = 6378206.4; //[m] 
constexpr double rearth_pole = 6378206.4; //[m]

constexpr double thetai = 0.6;

constexpr double eps = 1.019e-2; // Nudging factor XXX WHY THIS VALUE?

constexpr double dt = -1e10;  // Time step  XXX

// VECTOR/POINT operations ---------------------------------------------------------------------------

using Vector2 = double2;
using Vector3 = double3;
using Point2  = double2;
using Point3  = double3;

__device__ __host__ __forceinline__ double2 operator+(const double2& a, const double2& b) {
    return make_double2(a.x + b.x, a.y + b.y);
}

__device__ __host__ __forceinline__ double2 operator+=(double2& a, const double2& b) {
    a.x += b.x;
    a.y += b.y;
    return a;
}

__device__ __host__ __forceinline__ double2 operator-(const double2& a, const double2& b) {
    return make_double2(a.x - b.x, a.y - b.y);
}

__device__ __host__ __forceinline__ double2 operator*(const double2& v, double s) {
    return make_double2(v.x * s, v.y * s);
}

__device__ __host__ __forceinline__ double2 operator*(double s, const double2& v) {
    return v * s;  // reuse the previous operator
}

__device__ __host__ __forceinline__ double2 operator/(const double2& v, double s) {
    return make_double2(v.x / s, v.y / s);
}

__device__ __host__ __forceinline__ double dot(const double2& v, const double2& w) {
    return v.x*w.x + v.y*w.y;
}

__device__ __host__ __forceinline__ double cross(const double2& v, const double2& w) {
    return v.x*w.y - v.y*w.x;
}

__device__ __host__ __forceinline__ double3 operator+(const double3& a, const double3& b) {
    return make_double3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__device__ __host__ __forceinline__ double3 operator+=(double3& a, const double3& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

__device__ __host__ __forceinline__ double3 operator-(const double3& a, const double3& b) {
    return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__device__ __host__ __forceinline__ double3 operator*(const double3& v, double s) {
    return make_double3(v.x * s, v.y * s, v.z * s);
}

__device__ __host__ __forceinline__ double3 operator*(double s, const double3& v) {
    return v * s;
}

__device__ __host__ __forceinline__ double3 operator/(const double3& v, double s) {
    return make_double3(v.x / s, v.y / s, v.z / s);
}

__device__ __host__ __forceinline__ double dot(const double3& v, const double3& w) {
    return v.x*w.x + v.y*w.y + v.z*w.z;
}

__device__ __host__ __forceinline__ Vector2 toVector2(const double3& v) {
    return make_double2(v.x, v.y);
}

__device__ __host__ __forceinline__ Point2 toPoint2(const double3& p) {
    return make_double2(p.x, p.y);
}

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


struct Element {

   // Geometry information
   int idx;               // Index of the element
   int i34;               // Number of sides for the element 
   Side *sides[4];        // Sides of the element
   Node *nodes[4];        // Nodes of the element
   Element *nghbElems[4]; // Neighbors of the element

   
   double area;               // Area of the element
   double ssign[4];           // Sign associated with each side of an element
   double radius;             // Element equivalent radii

   double wind_angle;         // Wind rotation angle

   Frame frame;

   Vector2 gradPhi[4];       // Derivatives of shape functions

   int idry;                  // True if the element is dry
   int kbe;                   // Bottom index of the element

   double xc, yc, zc;         // Element center coordinates
   double xlon_el, ylat_el;   // Element center lat/lon coordinates in _degrees_
   double depth;              // Depth at element (min of all nodes)
   double z[nvrt];            // z-coord  (local frame - vertical up)

   // Physical information
   Tracers tracers[nvrt];     // Tracer values at each vertical level  

   Vector2 ghat1;             // Element averaged ghat (\breve{G} in notes)

   __inline__ __device__ bool isTri (void) const { return i34 == 3; }
   __inline__ __device__ bool isQuad(void) const { return i34 == 4; }

   __inline__ __device__ int nxq(int index) const { assert(index + 16*i34 >= 0); return (index + 16*i34) % i34; }
   __inline__ __device__ int nodeLocalIdx(Node *node) const { 
      int res = (nodes[0] == node) ? 0 : ((nodes[1] == node) ? 1 : ((nodes[2] == node) ? 2 : 3)); 
      assert(res >= 0 && res < i34 && nodes[res] == node && "nodeLocalIdx: node not found in element");
      return res;
   }

   __inline__ __device__ bool isEdgeBoundary(Side *side) const;

   __inline__ __device__ double computeSignedArea(void) const;

   __inline__ __device__ void nudgeToBarycenter(double &x, double &y) const {
      // Nudges the point (x,y) to the barycenter of the element
      x = (1-eps)*x + eps*xc;
      y = (1-eps)*y + eps*yc;
   }

   __device__ bool computeShapeCoords(Vector2 xt, double shpCoords[4], bool fixNegShapeCoords = false) const;

   __device__ void setFrame(void) {

      double ar1, ar2;
      computeLonLat({xc, yc, zc}, ar1, ar2);

      frame.update(ar1, ar2);        
    }

   __device__ bool findNextIntersection(Point2 x2, Side *&side, Element *&nextElement, double &t, 
                                        double& tIntersect);
   __device__ bool findNextIntersection(Point2 x1, Point2 x2, Side *&side, Element *&nextElement, double &t, 
                                        double& tIntersect);

   __inline__ __device__ bool containsPoint(Vector2 xt) const { 
      // Check if the point (xt,yt) is inside the element
      double shpCoords[4];
      bool res = computeShapeCoords(xt, shpCoords, false);
      return res;
   }

   __device__ double interpolate(const double shapeCoords[4], const double zAtNodes[4])
   {
      // Interpolate the value at the point (xt,yt) using the shape coordinates and the values at the nodes
      // This function works for both triangular and quadrilateral elements, since shapeCoords[3] is always 0.0 for triangular elements
      return shapeCoords[0]*zAtNodes[0] + shapeCoords[1]*zAtNodes[1] + shapeCoords[2]*zAtNodes[2] + shapeCoords[3]*zAtNodes[3];
      
   }

   __device__ void interpolateVelocities(Point3 xt, Vector3& u_int);

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


   __device__ double interpolate(int level, const double shapeCoords[4], double (Node::*func)(int) const) const {
      double values[4];
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
                                       double &dt, Point3 &xt, double zint, Vector3 &u_int);
    
   __device__ bool trackThroughElementsAndLevels(Point3 x0, Point3 &xt, Vector3 &u_int, 
                                                 Element *&nextElem, int &jlev, double &dt);

};


struct BoundarySegment {
   int ifltype; // Type of boundary condition
   
};

struct SideBoundary {
   Side *side; // Pointer to the side that this boundary segment belongs to
   BoundarySegment *segment; // Pointer to the boundary segment
  
};


struct Side {

   // Geometry information
   int idx;               // Index of the side
   Node *n1, *n2;
   Element *elems[2];     // Elements that share this side XXX (why not e1,e2?)

   double length;         // Length of the side. Used to be called distj; but I think this name is less clear now.
   double xc, yc;         // Location of side center. These variables used to have a j in their name: xcj, ycj, zcj
   double zc[nvrt];       // Side center z-coord at each vertical level
   double depth;          // Depth at side center

   Frame frame;

   SideBoundary *sideBnd; // Pointer to additional boundary information. For non-boundary sides, 
                          // this is nullptr

   
   int kb;                // Bottom index of the side
   bool idry;             // True if the side is dry
 
   Vector2 n;             // Side unit normal vector (for ics=2: these are done inside local lon/lat frame)

   // Physical information
   Vector2 u2[nvrt];            // Side unit normal vector (for ics=2: these are done inside local lon/lat frame)
   Vector2 U;                   // Side vertically averaged velocity vector

   double flux_adv_hface[nvrt]; // Flux of horizontal advection at each vertical level

   __inline__ __device__ double x(double t) const;
   __inline__ __device__ double y(double t) const;
   __inline__ __device__ Point2 point(double t) const { return {x(t), y(t)}; }
   __inline__ __device__ Point2 center(void) const { return {xc, yc}; }

   __forceinline__ __device__ Vector2 t(void) const {
      return {-n.y, n.x};
   }

   __inline__ __device__ bool isBnd(void) const {
      return sideBnd != nullptr;
   }

   // Get valid starting element (for backtracking)
   __inline__ __device__ Element* getStartElement() {
      return elems[0] ? elems[0] : elems[1];
   }

   __inline__ __device__ Element* findElementTowards(Vector2 u) const {
      // From the elements that share this side, find the one that is in the direction of the target point (xt,yt)

      // If one of the elements is null, return the other one
      if      (elems[0] == nullptr) return elems[1];
      else if (elems[1] == nullptr) return elems[0];


      // Dot product of the negative velocity vector with the side normal vector
      const double dot_prod1 = -dot(u, n);
      
      // Vector from side center to element 0 centroid
      const Vector2 v2 = {elems[0]->xc - xc, elems[0]->yc - yc};

      // Dot product of the vector with the side normal vector
      const double dot_prod2 = dot(v2, n);

      // If the dot product is positive, returns the element whose vector to the centroid is in the same direction as the target point
      // Otherwise, returns the other element (otherwise includes dot_prod1 == 0)

      if (dot_prod1*dot_prod2 > 0.0) return elems[0];
      else                           return elems[1];

 
   }

   __forceinline__ __device__ Element* otherElement(Element* elem) {
      // Given an element, returns the other element that shares this side
      return (elems[0] == elem) ? elems[1] : elems[0];
   }

   __forceinline__ __device__ Node* otherNode(Node* node) {
      // Given a node, returns the other node that shares this side
      return (n1 == node) ? n2 : n1;
   }
 
   __inline__ __device__ bool btrack(int jlev, double dt, double dtbk, double* sclr,
                                     Vector3 &u_int) const;




};


// This struct is used to store the boundary information for a node.
struct NodeBoundary {   
   BoundarySegment *segment;
   
   SideBoundary *sideBnd1; // Pointer to the boundary side object (can be nullptr)
   SideBoundary *sideBnd2; // Pointer to the boundary side object (can be nullptr)

   __forceinline__ __device__ SideBoundary *sideBnd(int index) { 
      assert(index == 0 || index == 1); 
      return (index == 0) ? sideBnd1 : sideBnd2; 
   }
};

struct Node {

   // Geometry information
   int idx;        // Index of the node
   Node **nodes;   // pointer to the nodes that are neighbors of this node (there is no guarantee that they are connected)
   int numNodes;   // number of nodes neighbors of this node

   Element **elements; // pointer to the elements that share this node
   int numElements;   // number of elements that share this node

   double x, y, z;
   Vector3 u2[nvrt];
   double znl[nvrt];
   double lon, lat;
   double depth;

   int kb;       // Bottom index of the node

   int idry;     // True if the node is dry

   NodeBoundary *nodeBnd; // Pointer to additional boundary information. For non-boundary nodes, 

   // Physical information
   Tracers tr_nd0[nvrt]; // Initial tracer conc. at nodes

   double eta1; // Elevation at nodes at previous timestep
   double eta2; // Elevation at nodes at current timestep

   double grav2; // Gravity 
   double hhat;  // \breve{H} in notes

   double qel; // RHS

   __forceinline__ __device__ Point2 point2(void) const { return {x, y}; }

   __host__ __device__ double sigmaTransitionDepth(void) const;

   __device__ void clearMatrixRow(void) {
      // Clear the matrix row for this node
      printf("NOT IMPLEMENTED: clearMatrixRow\n");
   }

   __device__ void addToSparseMatrix(Node *node, double value) {
      // Add the value to the sparse matrix for this node
      printf("NOT IMPLEMENTED: addToSparseMatrix\n");
   }

#ifdef BACKTRACKING
   __device__ void updateZCoordData() {
   // XXX: To be called before compZCoord(), preferibly right after Depth and/or eta2 are updated

      // Check if node is wet
      idry = (depth + eta2 <= h0);  // XXX Check if we are doing this elsewhere

      if (idry) return;

      if constexpr(ivcor == 2) {  // SZ coordinates
         
         // Compute bottom level index
         if (depth <= mesh->h_s) {
               kb = mesh->nzlevels;
         } else {
               kb = -1;
               for (int i = 0; i < mesh->kz-1; i++) {
                  if (-depth >= mesh->ztot[i] && -depth < mesh->ztot[i+1]) {
                     kb = i;
                     break;
                  }
               }
         }
      } 
      else if (ivcor == 1) {  // Localized sigma
         // Do nothing
      }

   }
#endif   

   __device__ double compZCoord(int k) const;



};




class Mesh {
   public:  

   int partNum; // Partition number for this mesh.

   // Pointers to each GPU's local memory (size is the number of partitions).
   const Node**    nodes;  
   const Element** elems;  
   const Side**    sides;  
   
   // Local memory
   Node*    localNodeData; 
   Element* localElemData;
   Side*    localSideData; 

   NodeBoundary *localNodeBoundaryMem; // Array of objects with the boundary information of each node
   SideBoundary *localSideBoundaryMem; // Array of objects with the boundary information of each side
   BoundarySegment *localSegmentMem;   // Array of objects with the boundary segments
   
   const int npgb;  // Global number of horizontal nodes
   const int negb;  // Global number of elements
   const int nsgb;  // Global number of sides

   const int np;   // Local number of horizontal nodes 
   const int ne;   // Local number of elements
   const int ns;   // Local number of sides

   const int max_np; // Maximum number of horizontal nodes per partition (shared by all partitions)
   const int max_ne; // Maximum number of elements per partition (shared by all partitions)
   const int max_ns; // Maximum number of sides per partition (shared by all partitions)

  
   // Additional topology information
   Node    **nodeNeighNodes; // Array of pointers to the nodes neighbors of each node
   Element **nodeNeighElems; // Array of pointers to the elements neighbors of each node



   // Vertical coordinate system
   // Sigma levels
   double h_s; 
   double h_c; 
   int    kz;
   double theta_b;
   double theta_f, sinh_theta_f;
   double *sigma_lcl;
   double *sigma;
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
                   int np, int ne, int ns, int nvrt, int kz,
                   int max_np, int max_ne, int max_ns,
                   const Node **nodesPtr, const Element **elemsPtr, const Side **sidesPtr,
                   Node *nodesMem, Element *elemsMem, Side *sidesMem,
                   Node    **localNodeNeighNodesMem,
                   Element **localNodeNeighElemsMem,
                   NodeBoundary *localNodeBoundaryMem,
                   SideBoundary *localSideBoundaryMem,
                   BoundarySegment *localSegmentMem,
                   double h_s, double theta_b, double theta_f, double *sigma_lcl, double *sigma, double *ztot,
                   double *cs, double *dcs
                   ) : partNum(partNum), npgb(npgb), negb(negb), nsgb(nsgb),
                       np(np), ne(ne), ns(ns),
                       max_np(max_np), max_ne(max_ne), max_ns(max_ns),
                       nodes(nodesPtr), elems(elemsPtr), sides(sidesPtr),
                       localNodeData(nodesMem), localElemData(elemsMem), localSideData(sidesMem),
                       nodeNeighNodes(localNodeNeighNodesMem), nodeNeighElems(localNodeNeighElemsMem),
                       localNodeBoundaryMem(localNodeBoundaryMem), localSideBoundaryMem(localSideBoundaryMem),
                       localSegmentMem(localSegmentMem),
                       h_s(h_s), theta_b(theta_b), theta_f(theta_f), sigma_lcl(sigma_lcl), sigma(sigma), ztot(ztot),
                       cs(cs), dcs(dcs)
   {
      printf("Mesh constructor for partition %d\n", partNum);
      
      sinh_theta_f = sinh(theta_f);

   }

   __device__ void populateMesh(int *i34, int *elnode, int *elside,
                   double *node_coords, double *wind_rotate_angle, 
                   int *isdel, int *isidenode, int *isbnd,
                   int *nond, int *nlnd, int *iond, int *ilnd,
                   int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                   int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                   int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                   int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                   int *localNodePart,     int *localElemPart,     int *localSidePart);


   __device__ ~Mesh(void)
   {
   }

   __device__ void test(double x, double y, double z)
   {
      printf("Testing mesh for partition %d, x=%f, y=%f\n", partNum, x, y);
      for (int i = 0; i < ne; i++) {
         auto isInside = localElemData[i].containsPoint({x, y});
         // printf("i = %d, isInside = %d, x = %f, y = %f\n", i, isInside, localElemData[i].xc, localElemData[i].yc);
         if (isInside) {
            printf("Point %f %f %f is inside element %d\n", x, y, z, i);
            Element *elem = &localElemData[i];
            
            Vector3 u_int;
            // auto res = elem->trackThroughElementsAndLevels(x,y,z, 
            //                                                x+0.01,y+0.01,z,
            //                                                uuint, vvint, wwint, elem);
            // Choose the first side
            Side *side = elem->sides[0];
            int jlev = 5;  
            // RECORDAR:
            // - con y=0 no encuentra el triangulo
            // - jlev = 5, z = 0.1...
            side->btrack(5, 3000.0, 3.1, nullptr, u_int);

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
   __global__ void getStructSizes(int* sizes)  {
      // Get the size of each struct. This is used to allocate memory for the structs in the host code.
       sizes[0] = (int)sizeof(Node);
       sizes[1] = (int)sizeof(Element);
       sizes[2] = (int)sizeof(Side);
       sizes[3] = (int)sizeof(NodeBoundary);
       sizes[4] = (int)sizeof(SideBoundary);
       sizes[5] = (int)sizeof(BoundarySegment);
   }
}