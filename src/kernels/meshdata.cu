// C++ POD struct that stores variable values at each node ot a FEM mesh (e.g. temperature, pressure, velocity, etc.). It is for
// an array of structs. This if for CUDA, and we are working on a extremely optimized model.



constexpr int N_TRACERS = 2;
constexpr int mvtx = 4;  // Maximum number of vertices per element (4 for quads)

// XXX These are to be set before compiling the code
constexpr int nvrt = 10; // Number of vertical levels


// Helper functions for computing the local node index at a given offset from another.
inline int nxQ(int offset, int index, int size) {
    return (index + offset) % size;
}

// Helper function for computing the local node index at a given offset from another.
inline int nxQ(int index, int size) {
    return index % size;
}


struct Tracers {
    double tracer[N_TRACERS];
};


// Forward declarations
struct Side;
struct Element;
struct Node;


struct Element {
   // Geometry information
   int i34;                 // Number of sides for the element 
   Side *sides[mvtx];      // Sides of the element
   Node *nodes[mvtx];      // Nodes of the element
   Element *nghbElems[mvtx]; // Neighbors of the element

   
   double area;               // Area of the element
   double ssign[mvtx];        // Sign associated with each side of an element
   double radiel;             // Element equivalent radii

   int idry;                 // True if the element is dry
   int kbe;                  // Bottom index of the element

   double xc, yc, zc;         // Element center coordinates
   double xlon_el, ylat_el;   // Element center lat/lon coordinates in _degrees_
   double depth;              // Depth at element (min of all nodes)
   double z[nvrt];            // z-coord  (local frame - vertical up)

   // Physical information
   Tracers tracers[nvrt];  // Tracer values at each vertical level  


   __inline__ __device__ bool isTri (void) const { return i34 == 3; }
   __inline__ __device__ bool isQuad(void) const { return i34 == 4; }

   __inline__ __device__ int nxQ(int index) const { return index % i34; }

   __inline__ __device__ bool btrack(
         int jlev,
         double dt,            // Tracking timestep
         double* sclr         // Array for backtracked values
   ) {
      const Element* currentElem = this;

      // Initialize outputs
      for(int i = 0; i < 4 + NUM_TRACERS; i++) {
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
         if (!currentElem->trackThroughElementsAndLevels(x0, y0, z0, 
                                                        xt, yt, zt, dtb,
                                                        uuint, vvint, wwint, 
                                                        jlev, nextElem)) {
            // Hit boundary or dry element
            time_rm = 0.0;
            break; // XXX This is not needed on GPU
         }

         // Check if velocity too small to continue
         const double vmag_2 = uuint*uuint + vvint*vvint;
         if(vmag_2 <= VELMIN_BTRACK_2) {   // XXX I think this is not needed on GPU
            time_rm = 0.0;
            break;
         }

         // Update remaining time and position
         time_rm -= dtb;
         if(time_rm <= 1e-6 * dt) break;

         x0 = xt;
         y0 = yt;
         z0 = zt;
      }

      // Calculate final interpolated values for sclr array
      // (Simplified for now - need to add proper interpolation)
      interpolateTracers(sclr, xt, yt, zt, jlev);
      
      
      return false; // Did not exit domain
   
      

   }


   /**
      For a point (x1,y1) inside the element and a line segment (x1,y1)-(x2,y2), find the intersection point 
      (xIntersect,yIntersect) and the side (intersectingSide) that it intersects.

      Inputs:
      @param x1, y1: point inside the element
      @param x2, y2: other end of the line segment

      Outputs:
      @param intersectingSide:  side that intersects the line segment
      @param nextElement: element that contains the intersection point
      @param xIntersect, yIntersect: intersection point  

      The equations for the two segments (the one we are considering and one of the sides of the element) are:
      P1(tt1) = (x1,y1) + tt1*(x2-x1, y2-y1)  where 0 ≤ tt1 ≤ 1
      P2(tt2) = (x3,y3) + tt2*(x4-x3, y4-y3)  where 0 ≤ tt2 ≤ 1

      If the two segments are parallel, they do not intersect. This is checked by the determinant of the two vectors
      (Delta) being zero.

      The intersection point happens when:
      tt1(x2-x1) - tt2(x4-x3) = x3-x1
      tt1(y2-y1) - tt2(y4-y3) = y3-y1
      
      Assumption: the 4 points (including those of the side segment) are different.
      The eqs. for the 2 lines are: X=X1+(X2-X1)*tt1 and X=X3+(X4-X3)*tt2.	
      Output: false: no intersection or colinear; true: exactly 1 intersection.	
      If true, (xIntersect,yIntersect) is the intersection.				
      Modified to only check tt2, assuming 0<=tt1<=1 is already assured. 
   **/
   __inline__ __device__ bool findPathIntersection(double x1, double y1, double x2, double y2, 
                                                   Side& intersectingSide, Element& nextElement,
                                                   double& xIntersect, double& yIntersect) {
      for (int i = 0; i < i34; i++) {
         const double x3 = nodes[i]->x;
         const double y3 = nodes[i]->y;
         const double x4 = nodes[nxQ(i+1)]->x;
         const double y4 = nodes[nxQ(i+1)]->y;

         const Side& side = sides[i];

         xIntersect=-1.e25;
         yIntersect=-1.e25;

         const double delta  = (x2-x1)*(y3-y4)-(y2-y1)*(x3-x4);
         const double delta1 = (x3-x1)*(y3-y4)-(y3-y1)*(x3-x4);
         const double delta2 = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1); 

         if (delta != 0.0) {
            // Solves the two linear equations for tt1 and tt2
            const double tt1 = delta1/delta; // delta is not zero, so we can divide
            const double tt2 = delta2/delta;

            // It is enough to check tt2, since the range [0,1] is the same for both tt1 and tt2.
            // Similarly, we use tt2 to compute the intersection point.
            if (tt2 >= -small && tt2 <= 1.0 + small) {
               xIntersect = x3 + (x4 - x3)*tt2;
               yIntersect = y3 + (y4 - y3)*tt2;

               intersectingSide = side;
               nextElement = nghbElems[i];

               return true;
            }
         
         }
         else {
            // The two segments are parallel
            return false;
         }

      }

      return false;

   }


   /**  Interpolate velocities at the final position

      @param shapeCoords: shape coordinates for the final position in the element.
      @param uuint, vvint, wwint: velocities interpolated at the final position (xt,yt,zt)

      This is done in four steps:

      1. Find the shape coordinates for the final position in the element.
      2. Interpolate the z level of the top and bottom levels at xt, yt. 
      3. Perform a *binary search* to find the vertical level containing zt. It will have to interpolate the 
         z levels at xt, yt at each step.   
      4. Interpolate the velocities at the final position (xt,yt,zt). This is done by interpolating the velocities
         horizontally and vertically.
   **/
   __inline__ __device__ bool interpolateVelocities(double *shapeCoordsOfPoint,  
                                                    double &uuint, double &vvint, double &wwint) {


   }



};


struct Side {
   // Geometry information
   Node *n1, *n2;
   Element *elems[2];     // Elements that share this side XXX (why not e1,e2?)

   double length;         // Length of the side. Used to be called distj; but I think this name is less clear now.
   double xc, yc;         // Location of side center. These variables used to have a j in their name: xcj, ycj, zcj
   double zc[nvrt];       // Side center z-coord at each vertical level
   double depth;          // Depth at side center

   int bdysd;             // Boundary side flag
   int kbs;               // Bottom index of the side
   bool idry;             // True if the side is dry
 
   double nx, ny;         // Side unit normal vector (for ics=2: these are done inside local lon/lat frame)

   double z[nvrt];        // Depth at side center.

   // Physical information
   double su2[nvrt], sv2[nvrt]; // Horizontal velocity

   double flux_adv_hface[nvrt]; // Flux of horizontal advection at each vertical level

   // Get valid starting element (for backtracking)
   __inline__ __device__ Element* getStartElement() {
      return elems[0] ? elems[0] : elems[1];
   }

   __inline__ __device__ bool btrack(
         int jlev,
         double dt,            // Tracking timestep
         double* sclr         // Array for backtracked values
   ) {
      // Set starting position to element center
      double xt = 0.0;
      double yt = 0.0;
      double zt = 0.0;
      
      // Call generic backtracking function
      return btrack(dt, xc, yc, zc, xt, yt, zt, jlev, sclr);
   }





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

   __device__ trackExitCodes trackThroughElementsAndLevels(Mesh& mesh,
                                                           double x0, double y0, double z0,
                                                           double xt, double yt, double zt, double dtb,
                                                           double uuint, double vvint, double wwint, 
                                                           int& jlev, Side& side, Element *nextElem)
   {

      // We start with the element in the direction of the target point.
      // Also, nudge the starting point towards the element centroid to avoid ambiguity
      Element *curElem = findElementTowards(side, x0, y0, xt, yt);

      // 1. Initial checks
      // Verify starting element is valid (wet)
      if (curElem->isDry) return tec_dry_element;  // Error code for invalid starting conditions

      // Checks if track has a zero length.
      double trm = time; //time remaining

      // Nudging weights
      constexpr double nudgeW1 = 0.01;  
      constexpr double nudgeW2 = 1.0 - nudgeW1;  


      double xctr3, yctr3;
      if (ics == 1) {
         double xctr3 = curElem->xc;
         double yctr3 = curElem->yc;
      } else {
         // TODO: Handle lat/lon projection case
         // project_pt('g2l', xctr, yctr, zctr, gcor0, frame0, xctr3, yctr3, tmp)
      }


      
      // When close to an element side, nudge starting point towards element centroid to avoid ambiguity (nudgeW1 + nudgeW2 = 1.0)
      // const double x0 = nudgeW2*x0 + nudgeW1*xctr3;
      // const double y0 = nudgeW2*y0 + nudgeW1*yctr3;

      double xcg = x0, ycg=y0;
      double pathl2 = (xt - xcg)*(xt - xcg) + (yt - ycg)*(yt - ycg);
      
      if (pathl2==0.0 || trm==0.0)  return tec_zero_length_track; 
      

      const int MAX_ITERATIONS = 1000;
      double xCur = x0, yCur = y0;        // Cur means current
      double xInt = -1e10, yInt = -1e10;  // Int means intersection
      
      // Horizontal tracking loop
      for (int iterations = 0; iterations < MAX_ITERATIONS; iterations++) {

         // If we are already in the target element, we are done.
          if (!curElem->containsPoint(xt, yt, zt)) break;


         // Find the intersection of the segment (xCur,yCur)-(xt,yt) with the element sides.
         int intersectingEdge;
         double xIntersect, yIntersect;
 

         // Find intersection with current element
         int intersectingEdge;
         double xIntersect, yIntersect;
         if (!curElem->findPathIntersection(xCur, yCur, xt, yt, 
                                            side, nextElem,
                                            xInt, yInt)) {
            // findPathIntersection is written in such a way that it must find an intersection.
            // If it doesn't find an intersection, then we have a problem.
            return tec_no_intersection;  
         }

         // Additional error conditions.
         if (curElem->isEdgeBoundary(intersectingEdge)) return tec_edge_boundary;
         if (nextElement->isDry) return tec_dry_element;
         

         // Update position and element for next iteration
         xCur = xInt;
         yCur = yInt;
         curElem = nextElem;
      }

      if (iterations >= MAX_ITERATIONS) return tec_iteration_limit;  
   

      // Compute shape coordinates for the final position
      double shapeCoords[4];  // For both tri (3) and quad (4) cases
      if (!currentElem->computeShapeCoords(xt, yt, shapeCoords)) {
         return tec_failed_to_compute_area_coords;  
      }

      
      // Interpolate velocities at the final position. This is a conceptually simple step,
      // but the function is not simple, because it has to search for the vertical level containing zt.
      if (!currentElem->interpolateVelocities(shapeCoords, jlev, 
                                             uuint, vvint, wwint)) {
         return tec_failed_to_interpolate_velocities;
      }

      return 0;
   }
   
};


struct Node {
   // Geometry information
   double x, y, z;
   double lon, lat;
   double depth;

   int idry;                 // True if the node is dry

   // Physical information
   // double w; // Vertical velocity

   Tracers tr_nd0[nvrt]; // Initial tracer conc. at nodes

   double eta1; // Elevation at nodes at previous timestep
   double eta2; // Elevation at nodes at current timestep
};



class Mesh {
   public:  
   const Node**    nodes;  // Array of pointers to NodeValues structs
   const Element** elems;  // Array of pointers to ElementValues structs
   const Side**    sides;  // Array of pointers to SideValues structs
   
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
   

   // Get's the pointer to the node using the array of global pointers to nodes
   __device__ const Node *node(int jn) {
      const int idx = jn / max_np; // Get the index of the partition that contains the node
      return nodes[idx] + idx;
   }

   // Get's the pointer to the element using the array of global pointers to elements
   __device__ const Element *elem(int je) {
      const int idx = je / max_ne; // Get the index of the partition that contains the element
      return elems[idx] + idx;
   }

   // Get's the pointer to the side using the array of global pointers to sides
   __device__ const Side *side(int js) {
      const int idx = js / max_ns; // Get the index of the partition that contains the side
      return sides[idx] + idx;
   }


  
   __device__ Mesh(int partNum,
       int npgb, int negb, int nsgb, 
       int np, int ne, int ns, int nvrt, 
       int max_np, int max_ne, int max_ns,
       const Node **nodesPtr, const Element **elemsPtr, const Side **sidesPtr,
       Node *_nodes, Element *_elems, Side *_sides, 
       int *i34, int *elnode, int *elside) : 
          npgb(npgb), negb(negb), nsgb(nsgb),
          np(np), ne(ne), ns(ns), 
          max_np(max_np), max_ne(max_ne), max_ns(max_ns),
          nodes(nodesPtr), elems(elemsPtr), sides(sidesPtr),
          localNodeData(_nodes), localElemData(_elems), localSideData(_sides)
          {

          for (int i = 0; i < ne; i++) {
             localElemData[i].i34 = i34[i];
             for (int j = 0; j < localElemData[i].i34; j++) {
                localElemData[i].sides[j] = &(localSideData[elside[i*mvtx + j]]);
                localElemData[i].nodes[j] = &(localNodeData[elnode[i*mvtx + j]]);
             }
             printf("i34[%d] = %d\n", i, i34[i]);
          }



   }

   __device__ ~Mesh() {

   }
};

// Each device has a global mesh object
__device__ Mesh* mesh = nullptr;



extern "C" {
   __global__ void getStructSizes(int* sizes)  {
      // Get the size of each struct. This is used to allocate memory for the structs in the host code.
       sizes[0] = (int)sizeof(Node);
       sizes[1] = (int)sizeof(Element);
       sizes[2] = (int)sizeof(Side);
   }
}





extern "C" {
__global__ void backtrackKernel(
    const Mesh* mesh,
    const double dt,
    BacktrackResults* results  // Structure to store results
) {
    // Get side index
    int sIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (sIdx >= mesh->ns) return;  // Check bounds for local sides
    
    // Get local side data
    const Side* side = mesh->localSideData + sIdx;
    if (side->isDry()) return;  // Skip dry sides
    
    // Get starting position (side center)
    double x0 = side->xc;
    double y0 = side->yc;
    double z0 = side->zc;
    
    // Get starting element
    Element* startElement = side->getStartElement();
    
    // Do backtracking
    double xt, yt, zt;  // Final positions
    bool exitedDomain = startElement->backtrack(
        x0, y0, z0,
        dt,
        xt, yt, zt
    );
    
    // Store results
    results[sid].finalPos = {xt, yt, zt};
    results[sid].exitedDomain = exitedDomain;
}


__device__ bool btrack(Mesh* mesh) {

     
};


extern "C" { 
   __global__ void computeFluxAdvHface(void) {
      // Get global thread ID
      int j = blockIdx.x * blockDim.x + threadIdx.x;
      
      if (j >= mesh->ns) return;  // Guard against out-of-bounds access

      // Get the side object
      Side *side = &(mesh->localSideData[j]);

      // Skip dry sides
      if (side->idry) return;
      
      // Check if it's land 
      const bool isLand = (side->elems[1] == 0 && side->bdysd <= 0);  // XXX What does isbs (bdysd) mean?
      
      // Get the bottom index for this side
      const int kbs = side->kbs;
      
      // Process each vertical level
      for (int k = kbs + 1; k < nvrt; k++) {
          if (isLand) {
              side->flux_adv_hface[k] = 0.0;

          } else {
              // Calculate normal component of the velocity at current and previous levels
              const double vnor1 = side->su2[k  ]*side->nx + side->sv2[k  ]*side->ny;
              const double vnor2 = side->su2[k-1]*side->nx + side->sv2[k-1]*side->ny;
              
              // Calculate flux (normal * area in local x-direction)
              const double height = side->zc[k] - side->zc[k-1]; // Vertical distance between the two side centers
              const double vmean = 0.5*(vnor1 + vnor2);
              const double length = side->length;
              
              side->flux_adv_hface[k] = height*length*vmean;
          }
      }
   }
}



__global__ void checkMassConservation(Tracers* sumtr) {

   int i = blockIdx.x * blockDim.x + threadIdx.x;

   if (i >= mesh->ne) return;  // Guard against out-of-bounds access

   Element *elem = &(mesh->localElemData[i]);

   // Skip dry elements
   if (elem->idry) return;

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

#ifdef DONTCOMPILE
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
//           endif !nproc>1
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



// factory functions
extern "C" 
__global__ void createMesh(int partNum, 
                           int npgb, int negb, int nsgb, int nvrt,
                           int np, int ne, int ns, 
                           int max_np, int max_ne, int max_ns,
                           const Node **nodesPtr, const Element **elemsPtr, const Side **sidesPtr,
                           Node    *nodes, Element *elems, Side *sides,
                           int *i34, int *elnode, int *elside) {

   if (mesh != nullptr) 
   {
      printf("ERROR: Mesh is already initialized\n");
      return;
   };

   mesh = new Mesh(partNum, npgb, negb, nsgb, nvrt, np, ne, ns, max_np, max_ne, max_ns,
                  //  start_p, start_e, start_s, end_p, end_e, end_s,
                   nodesPtr, elemsPtr, sidesPtr, nodes, elems, sides, i34, elnode, elside);
   
}

extern "C" 
__global__ void destroyMeshCPU(void) {
   delete mesh;
}
