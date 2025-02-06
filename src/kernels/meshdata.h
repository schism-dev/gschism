
// XXX: Pure constants that are to be set before compiling the code
constexpr int N_TRACERS = 2;

constexpr int ivcor = 1; // Vertical coordinate system
constexpr int ics = 1; // 1: global, 2: local

constexpr int nvrt = 10; // Number of vertical levels
constexpr int nzlev = 5; // Number of z levels (in contrast to sigma levels)
constexpr int nsiglev = nvrt - nzlev + 1; // Number of sigma levels (in contrast to z levels)


constexpr double small1 = 1.0e-10;       // Small non-negative number
constexpr double small2 = small1*100.0;  // Slightly larger than small1
constexpr double small = 0.0;            // Small positive number or 0


__inline__ __device__ double signedArea(double x1, double x2, double x3, double y1, double y2, double y3) {
  // Compute signed area formed by pts 1,2,3 (positive counter-clockwise)
  return 0.5*((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3));
}

__inline__ __device__ void Error(const char* msg) {
   printf("Error: %s\n", msg);
   exit(1);
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
   int i34;                 // Number of sides for the element 
   Side *sides[4];      // Sides of the element
   Node *nodes[4];      // Nodes of the element
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

   __inline__ __device__ double computeSignedArea(void) const {
      // Compute the area of the element

      double res = signedArea(nodes[0]->x, nodes[1]->x, nodes[2]->x, nodes[0]->y, nodes[1]->y, nodes[2]->y);

      if (isQuad()) {
         double ar1 = signedArea(nodes[0]->x, nodes[2]->x, nodes[3]->x, nodes[0]->y, nodes[2]->y, nodes[3]->y);
         double ar2 = signedArea(nodes[0]->x, nodes[1]->x, nodes[3]->x, nodes[0]->y, nodes[1]->y, nodes[3]->y);
         double ar3 = signedArea(nodes[1]->x, nodes[2]->x, nodes[3]->x, nodes[1]->y, nodes[2]->y, nodes[3]->y);
         if (ar1 < 0.0 || ar2 < 0.0 || ar3 < 0.0) {
            Error("Negative area found in quad element");
         }
         res += ar1;
      }

      return res;
   }

   __device__ bool findNextIntersection(double t, double x2, double y2, 
                                              Side *side, Element *nextElement,
                                              double& tIntersect);

   __inline__ __device__ bool containsPoint(double xt, double yt) const {
      // Check if the point (xt,yt) is inside the element
      printf("NOT IMPLEMENTED: containsPoint\n");
      return true;
   }


   __inline__ __device__ bool globalToLocal(double xt, double yt, double &xlocal, double &ylocal) const {
      // Convert global coordinates to local coordinates
      printf("NOT IMPLEMENTED: globalToLocal\n");
      return true;
   }

   __inline__ __device__ bool localToGlobal(double xlocal, double ylocal, double &xt, double &yt) const {
      // Convert local coordinates to global coordinates
      printf("NOT IMPLEMENTED: localToGlobal\n");
      return true;
   }

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
   Node *n1, *n2;
   Element *elems[2];     // Elements that share this side XXX (why not e1,e2?)

   double length;         // Length of the side. Used to be called distj; but I think this name is less clear now.
   double xc, yc;         // Location of side center. These variables used to have a j in their name: xcj, ycj, zcj
   double zc[nvrt];       // Side center z-coord at each vertical level
   double depth;          // Depth at side center

   int bdysd;             // Boundary side flag
   int kb;               // Bottom index of the side
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




   /*
      Tracks a particle through elements and vertical levels in the mesh.

      @param mesh         Reference to mesh object containing grid information
      @param x0,y0,z0     Starting position coordinates
      @param xt,yt,zt     Target position coordinates
      @param dtb          Time step for tracking
      @param uuint,vvint,wwint Initial velocity components
      @param jlev         Current/resulting vertical level
      @param nextElem     Starting/resulting element index

      @return Status code of type trackExitCodes: 

      Horizontal tracking:
      - We start at a side center (x0,y0)
      - We start on the element that contains the side center (nextElem). This is one of the elements that share this side.   
      - We need to find the element that contains the target point (xt,yt)

      At every step we have:
         - A current element (curElem)
         - A current position (xCur,yCur)
         - A target position (xt,yt) (always the same)
         - The side that we are currently on (At the start it is the side that contains the starting point, otherwise it is the side that contains the intersection point)

      Now, there are two cases: (a) the element at the other side of the current side contains the target point, or (b) it does not.
      In which case we know for sure that the target point is in a different element. And, for this reason we
      know for sure that the segment (xCur,yCur)-(xt,yt) intersects the element sides.
      
      - We continue this process until we hit a boundary or the element that contains the target point.
      - Once the target element is found, we compute the shape coordinates for the final position.

      Vertical tracking:
      - We then find the vertical level containing the target point.
      - We interpolate the velocities at the final position.
   */ 
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

#ifdef BACKTRACKING
   /**
    * Tracks a particle through elements and vertical levels in the mesh.
    * 
    * @param x0,y0,z0     Starting position coordinates
    * @param xt,yt,zt     Target position coordinates
    * @param dtb          Time step for tracking
    * @param uuint,vvint,wwint Initial velocity components
    * @param jlev         Current/resulting vertical level
    * @param nextElem     Starting/resulting element index
    * 
    * @return Status code of type trackExitCodes
    * 
    * This function uses the same ideas as the original Fortran code (function QuickSearch), but some
    * of the details are different. The idea is that the code requires less branches and is better fit
    * for GPU execution.
    * 
    * The main ideas here are:
    * - The starting point is guaranteed to be in the starting element.
    * - since it is inside the element, there are only two cases: (a) segment (x0,y0)-(xt,yt) intersects the element sides, or (b) it does not and (xt,yt) is in the same element.
    * - If it is not inside the element, we know that an interection must exist. We compute the intersection
    *   of the segment with the *lines* (infinite length) that form the element sides.
    * - The intersection is computed as the parametric value tt1 of the intersection point on the segment (x0,y0)-(xt,yt).
    * - tt1 can be positive or negative. We only consider the positive values.
    * - The next intersection will occur at the smallest positive tt1.
    * - This intersection will be the next (x0,y0) position.
    * - Proceeding this way it is guaranteed that we are always moving forward, so there is no chance
    *   of getting stuck in a loop.
    * - Additionally, this function never finds intersections with the side segment where the last intersection
    *   was found.
    */
   __device__ trackExitCodes trackThroughElementsAndLevels(double x0, double y0, double z0,
                                                           double xt, double yt, double zt, double dtb,
                                                           double uuint, double vvint, double wwint, 
                                                           int& jlev, Side *side, Element *nextElem)
   {

      // 1. Initial checks
      // Verify starting element is valid (wet)
      if (nextElem->isDry) return tec_dry_element;  // Error code for invalid starting conditions

      // Checks if track has a zero length.
      double trm = time; //time remaining


      double xcg = x0, ycg=y0;
      double pathl2 = (xt - xcg)*(xt - xcg) + (yt - ycg)*(yt - ycg);
      
      if (pathl2==0.0 || trm==0.0)  return tec_zero_length_track; 
      

      const int MAX_ITERATIONS = 1000;
      double xCur = x0, yCur = y0;        // Cur means current
      double xInt = -1e10, yInt = -1e10;  // Int means intersection
      
      double t = 0.5; // Initial t value for intersection search is in the middle of the side
      Side *nextSide = side;
      double tIntersect = 1e10;

      // Horizontal tracking loop
      for (int iterations = 0; iterations < MAX_ITERATIONS; iterations++) {

         bool found_intersection = curElem->findNextIntersection(t, xt, yt, nextSide, nextElement, tIntersect);

         // If we are already in the target element, we are done.
         // Remember that, if we are here, the point (xt,yt) is inside curElem.
         if (!found_intersection) break;

         // Additional error conditions.
         if (nextElement->isEdgeBoundary(nextSide)) return tec_edge_boundary;
         if (nextElement->isDry) return tec_dry_element;
         
         // Update position and element for next iteration
         side = nextSide;
         t = tIntersect;
      }


      if (iterations >= MAX_ITERATIONS) return tec_iteration_limit;  
   
      
      // Interpolate velocities at the final position. This is a conceptually simple step,
      // but the function is not simple, because it has to search for the vertical level containing zt.
      if (!currentElem->interpolateVelocities(shapeCoords, jlev, 
                                              uuint, vvint, wwint)) {
         return tec_failed_to_interpolate_velocities;
      }

      side = nullptr;

      return tec_success;
   }
#endif
};


struct Node {
   static Mesh* mesh;

   // Geometry information
   Node **nodes;   // pointer to the nodes connected to this node
   int numNodes;   // number of nodes connected to this node

   Element **elements; // pointer to the elements connected to this node
   int numElements;   // number of elements shared by this node

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


   __device__ __inline__ double Node::sigmaTransitionDepth(void) const;



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

   __device__ double compZCoord(int k) {
      // Computes the z coordinate at level k.
      // Depending on the vertical coordinate system, the z coordinate is computed differently.
      // When ivcor=1, the z coordinate is computed using the localized sigma levels.
      // When ivcor=2, the z coordinate is computed using the zeta/sigma levels. Z levels are
      // used 
      
      // NOTE: Assumes that updateZCoordData() has been called before this function.
        
      if (idry) {
         printf("compZCoord: dry location: %f %f\n", depth, eta2);
         return 0.0;
      }

      if constexpr(ivcor == 2) {
         // Special levels
         const double hmod2 = sigmaTransitionDepth(); // Depth transition between z and sigma levels

         if (k == mesh->kz)    return -hmod2;
         else if (k == nvrt-1) return eta2;
         else if (k >= kb) {
            const int ksig = k - mesh->nzlevels + 1;

            if (hmod2 <= mesh->h_c) {
               return mesh->sigma[ksig] * (hmod2 + eta2) + eta2;
            } 
            else if (eta2 <= - mesh->h_c - (depth - mesh->h_c) * mesh->theta_f/mesh->s_con1) {
               printf("ZCOOR: Please choose larger h_c: %f %f\n", eta2, mesh->h_c);
               return 0.0;
            } 
            else {
               return eta2*(1.0 + mesh->sigma[ksig]) + 
                      mesh->h_c * mesh->sigma[ksig] + 
                      (hmod2 - mesh->h_c) * mesh->cs[ksig];
            }
         }
      }
      else if (mesh->ivcor == 1) {
         if (k == kb) return -depth;
         else if (k == nvrt-1) return eta2;
         else if (k >= kb) {
               return (eta2 + depth) * node->sigma_lcl[k] + eta2;
         }
      }
      
      printf("ZCOOR: invalid level or coordinate system\n");
      return 0.0;
   }
#endif

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
   Node    **connectedNodes; // Array of pointers to the nodes connected to each node
   Element **sharedElements; // Array of pointers to the elements shared by each node

   // Vertical coordinate system
   // Sigma levels
   double h_s; 
   double h_c; 





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
                   int np, int ne, int ns, int nvrt,
                   int max_np, int max_ne, int max_ns,
                   const Node **nodesPtr, const Element **elemsPtr, const Side **sidesPtr,
                   Node *nodesMem, Element *elemsMem, Side *sidesMem,
                   Node    **localNodeNeighNodesMem,
                   Element **localNodeNeighElemsMem
                   ) : partNum(partNum), npgb(npgb), negb(negb), nsgb(nsgb),
                       np(np), ne(ne), ns(ns),
                       max_np(max_np), max_ne(max_ne), max_ns(max_ns),
                       nodes(nodesPtr), elems(elemsPtr), sides(sidesPtr),
                       localNodeData(nodesMem), localElemData(elemsMem), localSideData(sidesMem),
                       connectedNodes(localNodeNeighNodesMem), sharedElements(localNodeNeighElemsMem)
   {
      printf("Mesh constructor for partition %d\n", partNum);
   }

   __device__ void populateMesh(int *i34, int *elnode, int *elside,
                   double *node_coords, double *wind_rotate_angle, 
                   int *isdel, int *isidenode, int *isbnd,
                   int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                   int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                   int *localNodePart,     int *localElemPart,     int *localSidePart) {

      printf("Populating mesh for partition %d\n", partNum);

      printf("Processing %d elements\n", ne);
      for (int i = 0; i < ne; i++)
      {
         const int i34 = i34[i];

         Element &elem = localElemData[i];

         elem.i34 = i34;
         elem.wind_angle = wind_rotate_angle[i];

         elem.xc = 0.0;
         elem.yc = 0.0;
         elem.depth = 0.0;
         for (int j = 0; j < i34; j++) {
            elem.xc    += elem.nodes[j]->x/i34;
            elem.yc    += elem.nodes[j]->y/i34;
            elem.depth += elem.nodes[j]->z/i34;
         }

         // if constexpr(ics==2) elem.zc += elem.znd[j]/i34;
         // if(dp(elnode(j,ie))<dpe(ie)) dpe(ie)=dp(elnode(j,ie))
 
      

   
         elem.area = elem.computeSignedArea(); 
         elem.radius = sqrt(elem.area)/M_PI; 
         // Determines which way the normal vector to each side is oriented.
         // It only ensures that it is consistent for all elements. Otherwise, it is arbitrary.
         for (int j = 0; j < i34; j++)  
         {
            if (sides[j]->elems[0] == &elem) {
               elem.ssign[j] = 1;
            }
            else {
               elem.ssign[j] = -1;
            }

         }

      }


      printf("Processing %d nodes\n", np);      
      for (int i = 0; i < np; i++)
      {
         // localNodeData[i].nodes = connectedNodes[i];
         Node &node = localNodeData[i];

         node.x = node_coords[i * 3];
         node.y = node_coords[i * 3 + 1];
         node.z = node_coords[i * 3 + 2];
      }

      printf("Processing %d sides\n", ns);
      for (int i = 0; i < ns; i++)
      {
         // localNodeData[i].nodes = connectedNodes[i];
         Side &side = localSideData[i];

         // finding the pointer to the nodes that form the side
         const int idx_n1 = isidenode[globalSideIndices[i]*2];
         const int idx_n2 = isidenode[globalSideIndices[i]*2 + 1];

         const int locidx_n1 = localNodeIndices[idx_n1];
         const int locidx_n2 = localNodeIndices[idx_n2];
         const int part_n1 =   localNodePart[idx_n1];
         const int part_n2 =   localNodePart[idx_n2];

         side.n1 = (Node *)&(nodes[part_n1][locidx_n1]);
         side.n2 = (Node *)&(nodes[part_n2][locidx_n2]);

         if (i<10) printf("Side %d: %p %p\n", i, side.n1, side.n2);

         // finding the pointer to the elements that share the side
         const int idx_e1 = isdel[globalSideIndices[i]*2];
         const int idx_e2 = isdel[globalSideIndices[i]*2 + 1];

         const int locidx_e1 = localElemIndices[idx_e1];
         const int locidx_e2 = localElemIndices[idx_e2];
         const int part_e1 =  localElemPart[idx_e1];
         const int part_e2 =  localElemPart[idx_e2];

         // if (i<10) printf("Side %d: %d %d %d %d\n", i, idx_e1, idx_e2, part_e1, part_e2);
         side.elems[0] = (Element *)&(elems[part_e1][locidx_e1]);
         side.elems[1] = (Element *)&(elems[part_e2][locidx_e2]);

         side.xc = 0.5*(side.n1->x + side.n2->x);
         side.yc = 0.5*(side.n1->y + side.n2->y);
         for (int k = 0; k < nvrt; k++) side.zc[k] = 0.5*(side.n1->znl[k] + side.n2->znl[k]);

         side.length = sqrt((side.n2->x - side.n1->x)*(side.n2->x - side.n1->x) + 
                            (side.n2->y - side.n1->y)*(side.n2->y - side.n1->y));
      }
   }



   __device__ ~Mesh()
   {
   }
};

__device__ void Error(const char* msg) {
   // TODO: Implement error handling
}


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