
// XXX: Pure constants that are to be set before compiling the code
constexpr int N_TRACERS = 2;

constexpr int ivcor = 2; // Vertical coordinate system
constexpr int ics = 1; // 1: global, 2: local

constexpr int nvrt = 54; // Number of vertical levels
constexpr int nzlev = 5; // Number of z levels (in contrast to sigma levels)
constexpr int nsiglev = nvrt - nzlev + 1; // Number of sigma levels (in contrast to z levels)


constexpr double small1 = 1.0e-10;       // Small non-negative number
constexpr double small2 = small1*100.0;  // Slightly larger than small1
constexpr double small = 0.0;            // Small positive number or 0

constexpr int internalNode = -1;
constexpr int landNode = -2;

constexpr double pi = 3.14159265358979323846;

enum trackExitCodes {
   tec_success = 0,
   tec_dry_element = 1,
   tec_zero_length_track = 2,
   tec_iteration_limit = 3,
   tec_no_intersection = 4,
   tec_edge_boundary = 5,
   tec_no_element = 6,
   tec_failed_to_compute_area_coords = 7,
   tec_failed_to_find_vertical_level = 8,
   tec_failed_to_interpolate_velocities = 9
};


__inline__ __device__ double signedArea(double x1, double x2, double x3, double y1, double y2, double y3) {
  // Compute signed area formed by pts 1,2,3 (positive counter-clockwise)
  return 0.5*((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3));
}

__inline__ __device__ void Error(const char* msg) {
   printf("Error: %s\n", msg);
   assert(false);
}

struct Tracers {
    double tracer[N_TRACERS];
};


struct Frame {
}; 


// Forward declarations
struct Side;
struct Element;
struct Node;
struct Mesh;


struct Element {
   static Mesh* mesh;

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

   int idry;                  // True if the element is dry
   int kbe;                   // Bottom index of the element

   double xc, yc, zc;         // Element center coordinates
   double xlon_el, ylat_el;   // Element center lat/lon coordinates in _degrees_
   double depth;              // Depth at element (min of all nodes)
   double z[nvrt];            // z-coord  (local frame - vertical up)

   // Physical information
   Tracers tracers[nvrt];  // Tracer values at each vertical level  

    static __device__ void setMesh(Mesh* _mesh) { mesh = _mesh; }

   __inline__ __device__ bool isTri (void) const { return i34 == 3; }
   __inline__ __device__ bool isQuad(void) const { return i34 == 4; }

   __inline__ __device__ int nqx(int index) const { assert(index + 16*i34 >= 0); return (index + 16*i34) % i34; }

   __inline__ __device__ bool isEdgeBoundary(Side *side) const;

   __inline__ __device__ double computeSignedArea(void) const;

   __device__ bool computeShapeCoords(double &xt, double &yt, double shpCoords[4], bool fixNegShapeCoords = false) const;

   __device__ bool findNextIntersection(double x2, double y2, Side *&side, Element *&nextElement, double &t, double& tIntersect);
   __device__ bool findNextIntersection(double x1, double y1, double x2, double y2, 
                                        Side *&side, Element *&nextElement, double &t, double& tIntersect);

   __inline__ __device__ bool containsPoint(double xt, double yt) const { 
      // Check if the point (xt,yt) is inside the element
      double shpCoords[4];
      bool res = computeShapeCoords(xt, yt, shpCoords, false);
      return res;
   }

   __device__ double interpolate(const double shapeCoords[4], const double zAtNodes[4])
   {
      // Interpolate the value at the point (xt,yt) using the shape coordinates and the values at the nodes
      // This function works for both triangular and quadrilateral elements, since shapeCoords[3] is always 0.0 for triangular elements
      return shapeCoords[0]*zAtNodes[0] + shapeCoords[1]*zAtNodes[1] + shapeCoords[2]*zAtNodes[2] + shapeCoords[3]*zAtNodes[3];
      
   }
   __device__ bool interpolateVelocities(double xt, double yt, double zt, 
                                         double& uuint, double& vvint, double& wwint);


   // template<typename T>
   // __device__ T[4] getNodalValues(T Node::* field) const {
   //    // Get the values of the field at the nodes of the element
   //    T values[4];
   //    values[0] = nodes[0]->*field;
   //    values[1] = nodes[1]->*field;
   //    values[2] = nodes[2]->*field;
   //    if (isQuad()) values[3] = nodes[3]->*field;
   //    else          values[3] = 0.0;

   //    return values;
   //  }

   template<typename T>
   __device__ void getNodalValues(int level, T (Node::* field)[nvrt], T *res) const {
      // Get the values of the field at the nodes of the element
      res[0] = (nodes[0]->*field)[level];
      res[1] = (nodes[1]->*field)[level];
      res[2] = (nodes[2]->*field)[level];
      if (isQuad()) res[3] = (nodes[3]->*field)[level];
      else          res[3] = 0.0;

   }    

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


    
   __device__ trackExitCodes trackThroughElementsAndLevels(double x0, double y0, double z0,
                                                           double xt, double yt, double zt, double dtb,
                                                           double time,
                                                           double uuint, double vvint, double wwint, 
                                                           int& jlev, Side *side, Element *nextElem);
   
   


#ifdef BACKTRACKING
   __inline__ __device__ bool btrack(double xt, double yt, double zt,
                                     double uuint, double vvint, double wwint,
                                     int jlev,
                                     double dt, double dtbk,  // Tracking timesteps
                                     double time_rm,        // Remaining time
                                     double* sclr             // Array for backtracked values
   ) {
      const Element* curElem = this;

      // Initialize outputs
      for(int i = 0; i < 4 + N_TRACERS; i++) {
         sclr[i] = 0.0;
      }
        
      // Start Euler tracking
      double x0 = xt;
      double y0 = yt; 
      double z0 = zt;

      while(true) {
         // Calculate timestep
         const double dtb = min(dtbk, time_rm);

         // Euler step
         xt = x0 - dtb * uuint;
         yt = y0 - dtb * vvint;
         zt = z0 - dtb * wwint;

         Element* nextElem = nullptr;
         // Find containing element and interpolate velocities
         if (!curElem->trackThroughElementsAndLevels(x0, y0, z0, 
                                                     xt, yt, zt, dtb,
                                                     uuint, vvint, wwint, 
                                                     jlev, nextElem)) {
            // Hit boundary or dry element
            time_rm = 0.0;
            break; // XXX This is not needed on GPU???
         }

         // Check if velocity too small to continue
         const double vmag_2 = uuint*uuint + vvint*vvint;
         if (vmag_2 <= VELMIN_BTRACK_2) {   // XXX I think this is not needed on GPU
            time_rm = 0.0;
            break;
         }

         // Update remaining time and position
         time_rm -= dtb;
         if (time_rm <= 1e-6*dt) break;

         x0 = xt;
         y0 = yt;
         z0 = zt;
      }

      // Calculate final interpolated values for sclr array
      // (Simplified for now - need to add proper interpolation)
      interpolateTracers(sclr, xt, yt, zt, jlev);
      
      
      return false; // Did not exit domain
   
      

   }
#endif

};


struct Side {
   static Mesh* mesh;

   // Geometry information
   int idx;               // Index of the side
   Node *n1, *n2;
   Element *elems[2];     // Elements that share this side XXX (why not e1,e2?)

   double length;         // Length of the side. Used to be called distj; but I think this name is less clear now.
   double xc, yc;         // Location of side center. These variables used to have a j in their name: xcj, ycj, zcj
   double zc[nvrt];       // Side center z-coord at each vertical level
   double depth;          // Depth at side center

   int bdysd;             // Boundary side flag
   int kb;                // Bottom index of the side
   bool idry;             // True if the side is dry
 
   double nx, ny;         // Side unit normal vector (for ics=2: these are done inside local lon/lat frame)


   // Physical information
   double u2[nvrt], v2[nvrt]; // Horizontal velocity

   double flux_adv_hface[nvrt]; // Flux of horizontal advection at each vertical level

   static __device__ void setMesh(Mesh* _mesh) { mesh = _mesh; }

   // Get valid starting element (for backtracking)
   __inline__ __device__ Element* getStartElement() {
      return elems[0] ? elems[0] : elems[1];
   }

   __inline__ __device__ Element* findElementTowards(double xt, double yt) {
      // From the elements that share this side, find the one that is in the direction of the target point (xt,yt)

      // If one of the elements is null, return the other one
      if      (elems[0] == nullptr) return elems[1];
      else if (elems[1] == nullptr) return elems[0];

      // Vector from side center to target point
      const double vx = xt - xc;
      const double vy = yt - yc;

      // Dot product of the vector with the side normal vector
      const double dot_prod1 = vx*nx + vy*ny;
      
      // Vector from side center to element 0 centroid
      const double vx2 = elems[0]->xc - xc;
      const double vy2 = elems[0]->yc - yc;

      // Dot product of the vector with the side normal vector
      const double dot_prod2 = vx2*nx + vy2*ny;

      // If the dot product is positive, returns the element whose vector to the centroid is in the same direction as the target point
      // Otherwise, returns the other element (otherwise includes dot_prod1 == 0)

      if (dot_prod1*dot_prod2 > 0.0) return elems[0];
      else                           return elems[1];

  
  
  
   }

   __inline__ __device__ Element* otherElement(Element* elem) {
      // Given an element, returns the other element that shares this side
      return elems[0] == elem ? elems[1] : elems[0];
   }

#ifdef BACKTRACKING
   __inline__ __device__ bool btrack(
         int jlev,
         double dt,            // Tracking timestep
         double* sclr         // Array for backtracked values
   ) {
      // Start from side center
      double x0 = xc, y0 = yc, z0 = zc;

      
      Element* curElem = findElementTowards(xt, yt);
      
      // Initial velocities
      double u0 = u2[jlev];
      double v0 = v2[jlev];
      double w0 = 0.5*(n1->w2[jlev] + n2->w2[jlev]);  // vertical velocity at the side center

      Element* nextElem = curElem;
      while(true) {
         bool found = curElem->trackThroughElementsAndLevels(x0, y0, z0, xt, yt, zt, dt, u0, v0, w0, jlev, this, nextElem);
         if (found) break;
 

      }
   }
#endif








};


struct Node {
   static Mesh* mesh;

   // Geometry information
   int idx;        // Index of the node
   Node **nodes;   // pointer to the nodes that are neighbors of this node (there is no guarantee that they are connected)
   int numNodes;   // number of nodes neighbors of this node

   Element **elements; // pointer to the elements that share this node
   int numElements;   // number of elements that share this node

   double x, y, z;
   double u[nvrt], v[nvrt], w[nvrt];
   double znl[nvrt];
   double lon, lat;
   double depth;

   int kb;       // Bottom index of the node

   int idry;     // True if the node is dry

   // Physical information
   Tracers tr_nd0[nvrt]; // Initial tracer conc. at nodes

   double eta1; // Elevation at nodes at previous timestep
   double eta2; // Elevation at nodes at current timestep

   __host__ __device__ double sigmaTransitionDepth(void) const;



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
                   double h_s, double theta_b, double theta_f, double *sigma_lcl, double *sigma, double *ztot,
                   double *cs, double *dcs
                   ) : partNum(partNum), npgb(npgb), negb(negb), nsgb(nsgb),
                       np(np), ne(ne), ns(ns),
                       max_np(max_np), max_ne(max_ne), max_ns(max_ns),
                       nodes(nodesPtr), elems(elemsPtr), sides(sidesPtr),
                       localNodeData(nodesMem), localElemData(elemsMem), localSideData(sidesMem),
                       nodeNeighNodes(localNodeNeighNodesMem), nodeNeighElems(localNodeNeighElemsMem),
                       h_s(h_s), theta_b(theta_b), theta_f(theta_f), sigma_lcl(sigma_lcl), sigma(sigma), ztot(ztot),
                       cs(cs), dcs(dcs)
   {
      printf("Mesh constructor for partition %d\n", partNum);
      
      sinh_theta_f = sinh(theta_f);

   }

   __device__ void populateMesh(int *i34, int *elnode, int *elside,
                   double *node_coords, double *wind_rotate_angle, 
                   int *isdel, int *isidenode, int *isbnd,
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
      printf("Testing mesh for partition %d\n", partNum);
      for (int i = 0; i < ne; i++) {
         auto isInside = localElemData[i].containsPoint(x, y);
         if (isInside) {
            printf("Point %f %f %f is inside element %d\n", x, y, z, i);
            Element *elem = &localElemData[i];
            
            int jlev = 0;
            Side *side = nullptr;
            
            auto res = elem->trackThroughElementsAndLevels(x,y,z, x+0.01,y+0.01,z,
                                                           0.1, 1.0, 1.0,0.0,0.0,jlev, side, elem);
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
__device__ __inline__ int nqx(int offset, int index, int size) {
   // In c++, the % operator can return negative numbers, for this reason we add 16*size 
   // to the index to make it positive
   assert(index + offset + 16*size >= 0);
   return (index + offset + 16*size) % size;  
}

__device__ __inline__ int nqx(int index, int size) {
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
   }
}