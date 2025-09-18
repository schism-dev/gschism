#include "constants.h"
#include "meshdata.h"
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

__device__ void lock(int *mutex) {
    while (atomicCAS(mutex, 0, 1) != 0); // spin until lock acquired
}

__device__ void unlock(int *mutex) {
    atomicExch(mutex, 0);
}

__device__ void addDebugPlot(Vector3 x1, Vector3 x2, unsigned int color = 0xFF0000) {
   return;
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
      (*debugPlotIdx)++;
   }

   __threadfence();               // ensure all above writes are visible
   // Unlock the mutex
   unlock(&debugPlotMutex);

}

__device__ void addDebugPlot(Vector2 x1, Vector2 x2, unsigned int color = 0xFF0000) {
   addDebugPlot(toPoint3(x1, 0.0), toPoint3(x2, 0.0), color);
}

//===============================================================================
//     
//===============================================================================
__device__ __forceinline__ void update_bdef(double time, double x0, double y0, double &dep, Vector3 &vel) {
   // For imm=2, user needs to update bottom depth and velocity

   dep = min(1.0, 7.0 - (x0 + time));
   vel = { -1.0, 0.0, 0.0 };

}


__device__ __forceinline__ Vector3 Element::n(int k, double &area3D) const {
   // Computes the 3D normal vector to the element at level k
   // Returns the normal vector and the area of the element

      if (isTri()) {
         // For triangles...

         // Gets the normal vector to the element (with a length equal to twice the area)
         if constexpr(ics==1) {
            const int kb0 = max(k, nodes[0]->kb);
            const int kb1 = max(k, nodes[1]->kb);
            const int kb2 = max(k, nodes[2]->kb);

            const Vector3 n20 = nodes[2]->p3D(k) - nodes[0]->p3D(k);
            const Vector3 n10 = nodes[1]->p3D(k) - nodes[0]->p3D(k);
            
            // Computes the cross product of vectors (n2 - n0) and (n1 - n0). 
            const double nx = (n10.y)*(n20.z) - (n20.y)*(n10.z);
            const double ny = (n20.x)*(n10.z) - (n10.x)*(n20.z);
            const double nz = 2.0*area; // = (n10.x)*(n20.y) - (n20.x)*(n10.y);

            const double len = sqrt(nx*nx + ny*ny + nz*nz);
            area3D = 0.5*len;

            assert(area3D > 0.0 && "Element normal vector: area3D is negative");

            return {nx/len, ny/len, nz/len};

           
         }
         else {  // lat/lon
            assert(false && "Element normal vector: lat/lon not implemented for triangles");

         } // if ics
      }
      else {
         // For quads...

         assert(false && "Element normal vector: quads not yet implemented");


      }
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
__device__ bool Element::findNextIntersection(Point2 x2, Side *&side, Element *&nextElement,
                                              double &t, double& tIntersect) {

   double closest_intersect_t = 1e10;
   bool found_intersection = false;

   // (x1,y1)-(x2,y2) is the track line segment. Every step we move (x1,y1) to the next
   // intersection point, at position t on "side". The final position is (x2,y2)=(xt,yt).
   const Point2 x1 = side->point(t);
   

   Side *nextSide = nullptr;
   for (int i = 0; i < i34; i++) {

      // Element's i-th side 
      Side *side2 = sides[i];

      // Skip the side that contains the initial point
      if (side2 == side) continue;

      // (x3,y3)-(x4,y4) is the element side line
      const Point2 x3 = side2->n1->p;
      const Point2 x4 = side2->n2->p;

      // delta is 0 if the two lines are parallel. In this case, they do not intersect.
      const double delta = cross(x2 - x1, x3 - x4);

      if (delta != 0.0) {

         const double delta1 = cross(x3 - x1, x3 - x4);
         const double delta2 = cross(x2 - x1, x3 - x1); 

         // Solves the two linear equations for tt1 and tt2
         const double tt1 = delta1/delta; // delta is not zero, so we can divide
         const double tt2 = delta2/delta;

         if (tt1 > 0.0 && tt2 >= -small && tt2 <= 1.0 + small && tt1 <= 1.0) {
            closest_intersect_t = min(closest_intersect_t, tt1);
            nextSide = side2;

            nextElement = nextSide->otherElem(this);
            tIntersect = max(small, min(tt1, 1.0-small));
            t = tt2;

            found_intersection = true;
         }
      
      }

   }

   // found_intersection should always be true, since we checked all the sides.  
   side = nextSide;

   return found_intersection;
}




__device__ bool Element::findNextIntersection(Point2 x1, Point2 x2, Side *&side, Element *&nextElement,
                                              double &t, double& tIntersect) {

   double closest_intersect_t = 1e10;
   bool found_intersection = false;

   for (int i = 0; i < i34; i++) {

      Side *side2 = sides[i];

      // (x3,y3)-(x4,y4) is the element side line
      const Point2 x3 = side2->n1->p;
      const Point2 x4 = side2->n2->p;

      // delta is 0 if the two lines are parallel. In this case, they do not intersect.
      const double delta = cross(x2 - x1, x3 - x4);

      if (abs(delta) > 0) {

         const double delta1 = cross(x3 - x1, x3 - x4);
         const double delta2 = cross(x2 - x1, x3 - x1); 

         // Solves the two linear equations for tt1 and tt2
         const double tt1 = delta1/delta; // delta is not zero, so we can divide
         const double tt2 = delta2/delta;

         if (tt1 >= -small && tt2 >= -small && tt2 <= 1.0 + small && tt1 <= 1.0 + small) {
            closest_intersect_t = min(closest_intersect_t, tt1);
            side = side2;
            nextElement = side->otherElem(this);
            t = tt2;
            tIntersect = max(small, min(tt1, 1.0-small));
            found_intersection = true;
         }
      
      }
      else {  // parallel segments
         assert(false && "Parallel segments");
         // // In this case, we theck if they are colinear. if they are, we compute the intersection point as:
         // // tt2 += length(x1 - x2) * sign(dot(x1 - x2, x3 - x4)) in {0, 1}

         // // Check if the segments are colinear. It is enough to check if x1 is on the line of the side segment.

         // const double dist = dot(x1 - x3, side2->n);
         // if (abs(dist) < small*length(x3 - x4)) {
         //    const double2 unitarySide = normalize(x4 - x3);
         //    const double len = dot(x2 - x1, unitarySide);
         //    t = t + len;
         //    t = max(0.0, min(t, 1.0));
         //    printf("BEEEEEBB t = %.17e, dist = %.17e, len = %.17e, x1 = (%e,%e), x2 = (%e,%e), x3 = (%e,%e), x4 = (%e,%e)\n", t, dist, len, x1.x, x1.y, x2.x, x2.y, x3.x, x3.y, x4.x, x4.y);


         //    // tIntersect = max(small, min(tt1, 1.0-small));
         //    // t = tt2;
         //    found_intersection = true;
         // }
      }


   }

      
   assert(found_intersection && "No intersection found");

   // found_intersection should always be true, since we checked all the sides.  
   return found_intersection;
}


__device__ void Element::trackGoesOutOfWater(Side *nextSide, Element *curElem, double tt, int jlev,
                                             double &dt, Point3 &xt, double zint, Vector3 &u_int) {

   printf("Tracking segment goes out of water at element %i, nextSide = %i, xt = %.17g, yt = %.17g, ux = %.17g, uy = %.17g\n", curElem->idx, nextSide->idx, xt.x, xt.y, u_int.x, u_int.y);

   // When this happens, we stop our algorithm that usually assumes that we move to the next element.
   // Instead, we nudge the intersection point to the barycenter of the element, and nextElement becomes
   // curElement.

   // Then it finds the tangential velocity at the element side (the side where the intersection was found),
   // And uses this velocity instead to continue the tracking.

   // computes the intersection point in cartesian coordinates
   Point2 intersection = nextSide->p_t(tt);

   // Nudges the point to the barycenter of the element. 
   // This ensures that the point is inside the element and hence, we will find an intersection if the
   // next segment intersects the element sides.  
   curElem->nudgeToBarycenter(intersection);

   // Computes the tangential velocity at the side. Note that we are not using the x,y,z location
   // of the intersection point. We are approximating it by the velocity at the side center at jlev.
   // NOTE: With the new code it will be easy to compute the tangential velocity at the intersection point.
   double vtan = dot(nextSide->t(), nextSide->u2[jlev]);

   // If the side is an open boundary, the tangential velocity is zero.
   if (nextSide->isBnd()) vtan = 0.0;

   // The new velocity components are computed from the tangential velocity and the tangent vector.
   const double us = vtan*nextSide->t().x;
   const double vs = vtan*nextSide->t().y;
   const double ws = 0.5*(nextSide->n1->w2[jlev] + nextSide->n2->w2[jlev]);

   // Compute the new "tracked" time based on the remaining segment length and the original velocity.
   const double len2 = norm2(intersection - toPoint2(xt));
   const double vel2 = dot(u_int, u_int);        
   if (vel2 > 0.0) dt = dt - sqrt(len2/vel2); 
   else            dt = 0.0;

   // The end of this tracking segment is now the intersection point.
   xt = {intersection.x, intersection.y, zint};

   u_int = {us, vs, ws};
   
}

/**
 * Tracks a particle through elements and vertical levels in the mesh.
 * 
 * @param x0,y0,z0     Starting position coordinates
 * @param xt,yt,zt     Target position coordinates
 * @param dtb          Time step for tracking
 * @param uuint,vvint,wwint Initial velocity components
 * @param nextElem     Starting/resulting element index
 * 
 * @return Status code of type trackExitCodes
 * 
 * This function uses the same ideas as the original Fortran code (function QuickSearch), but some
 * of the details are different. The idea is that the code requires less branches and is better fit
 * for GPU execution.
 * 
 * The main ideas here are:
 * - The starting point is expected to be in the starting element.
 * - The end point is guaranteed to be inside an element.
 * - since the starting point is inside the element, there are only two cases: 
 *     (a) segment (x0,y0)-(xt,yt) intersects the element sides, or 
 *     (b) it does not and (xt,yt) is already in the starting element.
 * - If (xt,yt) is not inside the element, we know that an interection must exist. We compute the intersection
 *   of the segment with the *lines* (infinite length) that form the element sides.
 * - The intersection is computed as the parametric value tt1 of the intersection point on the segment (x0,y0)-(xt,yt).
 * - tt1 can be positive or negative. We only consider the positive values.
 * - tt2 is the parametric value of the intersection point on the element side segment.
 * - The next intersection will occur at the smallest positive tt1.
 * - This intersection will be the next (x0,y0) position.
 * - Proceeding this way it is guaranteed that we are always moving forward, so there is no chance
 *   of getting stuck in a loop.
 * - Additionally, this function never finds intersections with the side segment where the last intersection
 *   was found.
 */
__device__ bool Element::trackThroughElementsAndLevels(Point3 x0, Point3 &xt, Vector3 &u_int, 
                                                       Element *&nextElem, int &jlev, double &dt,
                                                       const Side *origSide) {



   // First thing is to check if the target point is inside the element. If it is, there is nothing
   // else to do.
   // If it is, nudges the point to be inside the element.
   const Element *curElem = nextElem;
   if (containsPoint(toPoint2(xt), true)) {

      curElem->interpolateVelocities(xt, u_int);

addDebugPlot(x0, xt);
      return true;
   }

   Side *side = nullptr;

   // 1. Initial checks
   // Verify starting element is valid (wet)
   assert(!curElem->isDry && "Starting element is dry");

   // Checks if track has a zero length.
   const double len2 = (xt.x - x0.x)*(xt.x - x0.x) + (xt.y - x0.y)*(xt.y - x0.y);   
   // assert(len2 > 0.0 && "Zero length track"); 
   
   const int MAX_ITERATIONS = 1000;
   
   double t;
   Side *nextSide = nullptr;
   double tIntersect = 1e10;
   double zint = x0.z;

   // Horizontal tracking loop
   // Finds the first intersection with the element sides. It is guaranteed to exist because
   // (x0,y0) is inside the element but the target point is not.

   bool found_intersection = findNextIntersection(toPoint2(x0), toPoint2(xt), nextSide, nextElem, t, tIntersect);

   assert(found_intersection && "No intersection found at %d");

   zint = (1 - tIntersect)*zint + tIntersect*xt.z;

   // if the tracking segment goes out of water (either dry or edge boundary)
   if (isEdgeBoundary(nextSide) || (nextElem == nullptr || nextElem->isDry)) {

addDebugPlot(x0, xt, 0x00A0FF);
      trackGoesOutOfWater(nextSide, this, t, jlev, dt, xt, zint, u_int);
addDebugPlot(x0, xt, 0x00A0FF);      

      nextElem = this;
      // Returns "early", meaning that we stopped before we were able to reach the target point.
      // In this case, the caller (btrack) needs to take into account the remaining time.
      return true;    
   }


   // Keeps searching for the next intersection until we find xt, yt inside nextElem.
   int iterations = 0;
   while (iterations < MAX_ITERATIONS) {
      iterations++;
      
      Element *curElem = nextElem;

      bool found_intersection = curElem->findNextIntersection(toPoint2(xt), nextSide, nextElem, t, tIntersect);
addDebugPlot(x0, xt);
      zint = (1 - tIntersect)*zint + tIntersect*xt.z;

      // If we are already in the target element, we are done.
      // Remember that, if we are here, the point (xt,yt) is inside curElem.
      if (!found_intersection) break;



      // if the tracking segment goes out of water (either dry or edge boundary).
      if (curElem->isEdgeBoundary(nextSide) || (nextElem == nullptr || nextElem->isDry)) {
addDebugPlot(x0, xt, 0x00A0FF);
         trackGoesOutOfWater(nextSide, curElem, t, jlev, dt, xt, zint, u_int);
addDebugPlot(x0, xt, 0x00A0FF);   
         nextElem = curElem;

         // Returns "early", meaning that we stopped before we were able to reach the target point.
         // In this case, the caller (btrack) needs to take into account the remaining time.
         return true;
      }
      
      // Update position and element for next iteration
      side = nextSide;

   }

   assert(iterations < MAX_ITERATIONS && "Iteration limit reached");  

 
   // Interpolate velocities at the final position. This is a conceptually simple step,
   // but the function is not simple, because it has to search for the vertical level containing zt.
addDebugPlot(x0, xt);   
   nextElem->interpolateVelocities(xt, u_int); 

   return true;
}


__inline__ __device__ double Element::computeSignedArea(void) const {
   // Compute the area of the element

   double res = signedArea(nodes[0]->p, nodes[1]->p, nodes[2]->p);

   // In the quad case we need to add the area of the other triangle.
   if (isQuad()) {
      double ar1 = signedArea(nodes[0]->p, nodes[2]->p, nodes[3]->p);
      double ar2 = signedArea(nodes[0]->p, nodes[1]->p, nodes[3]->p);
      double ar3 = signedArea(nodes[1]->p, nodes[2]->p, nodes[3]->p);
      assert (ar1 >= 0.0 && ar2 >= 0.0 && ar3 >= 0.0 && "Negative area found in quad element");

      res += ar1;
   }

   return res;
}


__inline__ __device__ bool Element::isDry2(void) const { return (nodes[0]->isDry2()) ||
                                                                (nodes[1]->isDry2()) ||
                                                                (nodes[2]->isDry2()) ||
                                                                (i34==4 && nodes[3]->isDry2());
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



/**
   Compute the shape coordinates of the point (xt,yt) in the element.
   @param xt (input): x-coordinate of the point
   @param shpCoords (output): shape coordinates of the point
   @param fixNegshpCoord (input): if true, fix negative shape coordinates
   @return true if the point is inside the element, false otherwise
**/
__device__ bool Element::computeShapeCoords(Point2 &xt, double shpCoords[4],
                                            bool fixNegshpCoord) const {

   // Compute the shape coordinates of the point xt in the element.
   if (i34==3) {
      // Compute area coordinates for a given pt w.r.t. to a triangular element
      // If ifl=1, will fix 0 or negative area coord. (assuming it's not too negative)
      // and in this case, the pt will be nudged into the element
      // It's not reliable to use the area coord for ics=2 when the pt is far away
      // from the local frame, so use additional check.

      Point2 x0, x1, x2;

      if constexpr(ics==1) {
         x0 = nodes[0]->p;
         x1 = nodes[1]->p;
         x2 = nodes[2]->p;
      }
      else // lat/lon
      {
         Error("Not implemented yet");
#ifdef ICS2
         globalToLocal(nodes[0]->x, nodes[0]->y, nodes[0]->z, x0, y0, z0);
         globalToLocal(nodes[1]->x, nodes[1]->y, nodes[1]->z, x1, y1, z1);
         globalToLocal(nodes[2]->x, nodes[2]->y, nodes[2]->z, x2, y2, z2);
#endif
      }

      shpCoords[0] = signedArea(xt,x1,x2)/area;
      shpCoords[1] = signedArea(x0,xt,x2)/area; 
      shpCoords[2] = 1.0 - shpCoords[0] - shpCoords[1];
      shpCoords[3] = 0.0;
      double tmpmin  = min3(shpCoords[0], shpCoords[1], shpCoords[2]);

      // If at least one of the shape coordinates is negative, and we have the flag fixNegshpCoord set,
      // then we fix the negative shape coordinates.
      if (fixNegshpCoord && tmpmin<=0.0) {

         // NOTE: I have changed this algorithm. Now it scales and offsets all the shape coordinates by the same factor,
         // keeping the sum of the shape coordinates equal to 1.
         // I have tested and compared the new algorithm with the old one, and they give
         // very similar results in terms of how close the "fixed" point is to the original point in cartesian coordinates.
         // The new algorithm should be faster since it is branchless.
         const double a = - tmpmin + 1.0e-4;  // XXX it used to be 1.0e-4 in the comments of the original code
         const double b = 1.0/(shpCoords[0] + shpCoords[1] + shpCoords[2] + 3.0*a);
         
         shpCoords[0] = (shpCoords[0] + a)*b;
         shpCoords[1] = (shpCoords[1] + a)*b;
         shpCoords[2] = (shpCoords[2] + a)*b;          

         if (shpCoords[0] < 0.0 || shpCoords[1] < 0.0 || shpCoords[2] < 0.0 ||
             shpCoords[0] > 1.0 || shpCoords[1] > 1.0 || shpCoords[2] > 1.0) {
            Error("Negative shape coordinates in computeshpCoord, when fixNegshpCoord is true");
         }

         // Update the point coordinates for the new shape coordinates
         double xcoords[3] = {x0.x, x1.x, x2.x};
         double ycoords[3] = {x0.y, x1.y, x2.y};
         xt = {dotprod3(xcoords, shpCoords), dotprod3(ycoords, shpCoords)};

      }

      // Check if the point is inside the element
      return min(min(shpCoords[0], shpCoords[1]), shpCoords[2]) > -small;      
   }
   else {
      Error("Not implemented yet");
#ifdef QUADS
      // Reproject pt in eframe of nel for ics=2
      if constexpr(ics==1) {
         xn2 = xt; 
         yn2 = yt;
      }
      else {
         localToGlobal(xt, yt, gcor0, xcg2, ycg2, zcg2)
         globalToLocal(xcg2, ycg2, zcg2, xctr, yctr, zctr, xn2, yn2, zn2)
      }

      return quad_shape(1,6,nnel,xn2,yn2,inside2,arco)
#endif
   }

   return false;  // Never reached

}


__device__ void Element::shapeCoordsToXY(double shpCoords[4], Vector2 &xt) const {
// Convert shape coordinates to cartesian coordinates

   if (i34 == 3) {
      // Triangular element
      double xcoords[3] = {nodes[0]->p.x, nodes[1]->p.x, nodes[2]->p.x};
      double ycoords[3] = {nodes[0]->p.y, nodes[1]->p.y, nodes[2]->p.y};
   }
   else {
      assert(false && "shapeCoordsToXY: quadrilateral element not implemented");
   }

}

__device__ double Element::intPhiPhi(int i1, int i2) const {
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

   double xi_term  = 1.0 + (ixi_n[i1]*ixi_n[i2])/3.0;
   double eta_term = 1.0 + (iet_n[i1]*iet_n[i2])/3.0;
   double eta_sum  = iet_n[i1] + iet_n[i2];
   double xi_sum   = ixi_n[i1] + ixi_n[i2];

   return (1.0/16.0)*area*xi_term*eta_term + (detJ_xi()/12.0)*xi_term*eta_sum + (detJ_eta()/12.0)*eta_term*xi_sum;
}
   

__device__ double Element::intGradPhiGradPhi(int i1, int i2) const {
   // return \int (grad phi_i1).(grad phi_i2) dA, where phi_i1 and phi_i2 are the shape functions of the nodes i1 and i2.
   // The integral is computed by using the 2x2 Gauss-Legendre quadrature points in xi and eta coordinates.
   // Again, it uses the determinant of the Jacobian to compute the integral in x,y coordinates.
   // The formula is also based on the identity:
   //   dPhi/dx = (1/det J) * [ (dPhi/dxi)*(dy/deta) - (dPhi/deta)*(dy/dxi)]
   //   dPhi/dy = (1/det J) * [-(dPhi/dxi)*(dx/deta) + (dPhi/deta)*(dx/dxi)]

   printf("intGradPhiGradPhi: TEST\n");   

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
__device__ void Element::interpolateVelocities(Point3 xt, Vector3& u_int) const {
   // Interpolate the velocities at the point xt. 
   // This function works the same for triangular and quadrilateral elements.
   double shpCoord[4];

   // First step is to compute the shape coordinates of (xt,yt)
   computeShapeCoords(toVector2(xt), shpCoord, false);
   assert((min4(shpCoord[0], shpCoord[1], shpCoord[2], shpCoord[3]) >= -small) && "Invalid shape coordinates");

   // Second step is to find the level that contains zt
   // Uses a binary search to find the level
   // Note that z coordinates are interpolated from their values at the nodes of the element.
   int L = kb;   // I use L(eft) and R(ight) as usual, but here L is the bottom level and R is the top level
   int R = nvrt - 1;

   // Get the z coordinates at the nodes of the element at the levels L and R
   double zLAtNodes[4], zRAtNodes[4];
   getNodalValues(L, &Node::compZCoord, zLAtNodes);
   getNodalValues(R, &Node::compZCoord, zRAtNodes);

   double zL = interpolate(shpCoord, zLAtNodes);
   double zR = interpolate(shpCoord, zRAtNodes);
   
   // Interpolates vertically the velocities
   double wght1, wght2;
   if (xt.z<=zL) {
      wght1 = 1.0;
      wght2 = 0.0;
   }
   else if (xt.z>=zR) {
      wght1 = 0.0;
      wght2 = 1.0;
   }
   else {
      // Perform a binary search to find the levels between L and R that contain xt.z
      while (R - L > 1) {
         int M = (L + R) / 2;

         double zMAtNodes[4];
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
      assert(L == R-1);
      // At this point, L and R are the two levels that contain zt, and zL and zR are the z coordinates at the nodes of the element at these levels

      wght2 = (xt.z - zL)/(zR - zL);
      wght1 = 1.0 - wght2;
   }

   // Interpolate horizontally and vertically the velocities at the two levels L, R
   u_int = wght1*interpolate(L, shpCoord, &Node::uvw2) + wght2*interpolate(R, shpCoord, &Node::uvw2);

   // Neither of the interpolated velocities should be NaN
   assert(!isnan(u_int.x) && !isnan(u_int.y) && !isnan(u_int.z) && "NaN velocities in interpolateVelocities");

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



__inline__ __device__ double  Side::x(double t)   const { assert(t >= 0.0 && t <= 1.0); return (1-t)*n1->p.x + t*n2->p.x; }
__inline__ __device__ double  Side::y(double t)   const { assert(t >= 0.0 && t <= 1.0); return (1-t)*n1->p.y + t*n2->p.y; }
__inline__ __device__ double  Side::lon(void) const { return 0.5*(n1->lon + n2->lon); }
__inline__ __device__ double  Side::lat(void) const { return 0.5*(n1->lat + n2->lat); }
__inline__ __device__ Point2  Side::p_t(double t) const { return {x(t), y(t)}; }
__inline__ __device__ Vector2 Side::grad2(double Node::* field, bool (*condition)(const Element*)) const {
   Vector2 grad = {0.0, 0.0};
   int count = 0;

   if constexpr(ics == 2) assert(false && "grad2: ics = 2 not implemented");

   // First element
   if (elems[0] != nullptr && condition(elems[0])) {
      grad += elems[0]->grad2(field);
      count++;
   }

   // Second element
   if (elems[1] != nullptr && condition(elems[1])) {
      grad += elems[1]->grad2(field);
      count++;
   }

   if (count > 0) {
      grad = grad/count;
   }
   
   return grad;
}

__inline__ __device__ double  Side::grav3(void)   const { return 0.5*(n1->grav2  + n2->grav2); }
__inline__ __device__ int     Side::iadv(void)    const { return min(n1->iadv,     n2->iadv); }   
__inline__ __device__ double  Side::eta2(void)    const { return 0.5*(n1->eta2   + n2->eta2); }
__inline__ __device__ double  Side::eta1(void)    const { return 0.5*(n1->eta1   + n2->eta1); }
__inline__ __device__ double  Side::depth(void)   const { return 0.5*(n1->depth  + n2->depth); } 
__inline__ __device__ double  Side::dfv(int k)    const { return 0.5*(n1->dfv[k] + n2->dfv[k]); }

__inline__ __device__ double  Side::veg_h(void)   const { return 0.5*(n1->veg_h  + n2->veg_h); }
__inline__ __device__ double  Side::veg_nv(void)  const { return 0.5*(n1->veg_nv + n2->veg_nv); }
__inline__ __device__ double  Side::veg_alpha_vert_mean(void) const { return 0.5*(n1->veg_alpha_vert_mean + n2->veg_alpha_vert_mean); }
__inline__ __device__ double  Side::veg_alpha3D(int k) const { return 0.5*(n1->veg_alpha3D[k] + n2->veg_alpha3D[k]); }


__inline__ __device__ double  Side::Cd(void)     const { if (!isDry) return 0.5*(n1->Cd + n2->Cd); else return 0.0; }

__inline__ __device__ Vector2 Side::tau(void)    const { return 0.5*(n1->tau + n2->tau); }

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
      int jlev,
      double dt,            // Tracking time
      double dtbk,          // Tracking timestep
      Vector3 &u_int        // Velocities at the initial point
) const {

   // Start from side center
   Point3 x0 = {center.x, center.y, zc[jlev]};

 
   double time_rm = dt; // Remaining time


   // // Initialize outputs XXX
   // for(int i = 0; i < 4 + N_TRACERS; i++) {
   //    sclr[i] = 0.0;
   // }


   // Initial velocities
   u_int = {u2[jlev].x, u2[jlev].y, 0.5*(n1->w2[jlev] + n2->w2[jlev])};

   // Find the element that is in the direction opposite to the velocity
   Element* curElem = findElementTowards(toVector2(u_int));


   // If it starts at a boundary side, we must start with a velocity parallel to the side.
   if (curElem == nullptr) {
      // Make velocity parallel to the side
      const Vector2 u2parallel = dot(toVector2(u_int), t())*t();
      // u_int = {u2parallel.x, u2parallel.y, u_int.z};
      curElem = otherElem(curElem);
      
   }

// Nudges it towards the circumcenter
// printf("0{{jlev = %d}} side->idx = %d, x0 = %.15e,%.15e,%.15e\n", jlev, idx, x0.x, x0.y, x0.z);
Vector2 xcirc = {(curElem->nodes[0]->p.x + curElem->nodes[1]->p.x + curElem->nodes[2]->p.x)/3.0,
                 (curElem->nodes[0]->p.y + curElem->nodes[1]->p.y + curElem->nodes[2]->p.y)/3.0};
const double eps = 9.013e-3;                 
x0.x = xcirc.x*eps + (1 - eps)*x0.x;
x0.y = xcirc.y*eps + (1 - eps)*x0.y;



   const int MAX_ITERATIONS = 1000; 
   Element* nextElem = curElem;
   int iter = 0;
   while(true) {

      // Calculate timestep
      double dtb = min(dtbk, time_rm);

      // Next point in the backtracking curve
      Point3 xt = x0 - dtb*u_int;
      bool found = curElem->trackThroughElementsAndLevels(x0, xt, u_int, nextElem, jlev, dtb, this);
      assert(found && "Backtracking failed in btrack");

      if (time_rm == 0.0) break; // from the definition of dtb, time_rm will reach exactly 0.0

      time_rm -= dtb; // if trackThroughElementsAndLevels terminates early (if it intersects a boundary or a dry element),
                      // then time_rm needs to be increased by the unused time.

      // If the velocity is too small, we stop the backtracking
      const double vmag = length(toVector2(u_int));
      if (vmag <= velmin_btrack) break;

      // The current point becomes the next segment's starting point.
      x0 = xt;
 
      curElem = nextElem;

      iter++;
      assert(iter < MAX_ITERATIONS && "Max iterations reached in btrack");
   }


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


__device__ void Node::updateKb(void) {
// Updates the node's kb index. Must be called before compZCoord(), preferibly right after Depth and/or eta2 are updated

   if constexpr(ivcor == 1) return;


   if (isDry2()) return;

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
      }
   } 
}


__device__ __inline__ void Node::updateZCoordData(void) {
// XXX: To be called before compZCoord(), preferibly right after Depth and/or eta2 are updated

   if (isDry2()) {
      kb = -1;
      return;
   }

   if constexpr(ivcor == 2) {  // SZ coordinates
      
      // Compute bottom level index
      if (depth <= mesh->h_s) {
            kb = izlev;

      } else {
         kb = -1;
         for (int i = 0; i < nzlev; i++) {
            // XXX Use binary search
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


__device__ double Node::compZCoord(int k) const {
   // Computes the z coordinate at level k.
   // Depending on the vertical coordinate system, the z coordinate is computed differently.
   // When ivcor=1, the z coordinate is computed using the localized sigma levels.
   // When ivcor=2, the z coordinate is computed using the zeta/sigma levels. Z levels are
   // used 
   
   // NOTE: Assumes that updateZCoordData() has been called before this function.

   if (isDry2()) return 0.0;

   assert(htot2() > h0 && "ZCOOR: invalid depth\n");

   if constexpr(ivcor == 2) {

      const double hmod2 = sigmaTransitionDepth(); // Depth transition between z and sigma levels
      const double h_c = mesh->h_c;
      const double theta_f = mesh->theta_f;
      const double sinh_theta_f = mesh->sinh_theta_f;

      // TODO: Here we are assuming that the bottom level, kbpl, was found in a previous step.
      // The FORTAN code computes it every time this function is called, but it should be done only 
      // once per time step/node.
      // CHECK: updateZCoordData

      //  To simplify comparisons, I define index versions of count variables.
      if      (k == ivrt)           return eta2; 
      else if (k == izlev)          return -hmod2;
      else if (k == kb)             return -depth;
      else if (k < izlev && k >= 0) return mesh->ztot[k];
      else if (k > izlev) {
         const int ksig = k - izlev;

         if (hmod2 <= h_c) {

            return mesh->sigma[ksig] * (hmod2 + eta2) + eta2;
         } 
         else if (eta2 <= - h_c - (depth - h_c) * theta_f/sinh_theta_f) {
            printf("Please choose larger h_c: %f %f\n", eta2, h_c);
            assert(false && "ZCOOR: invalid h_c\n");
         } 
         else {
            return eta2*(1.0 + mesh->sigma[ksig]) + 
                         h_c * mesh->sigma[ksig] + 
                         (hmod2 - h_c) * mesh->cs[ksig];
         }
      } else { 
         assert(false && "ZCOOR: invalid level or coordinate system\n");
      }

   }
   else if constexpr(ivcor == 1) {
      if (k == kb) return -depth;
      else if (k == nvrt-1) return eta2;
      else if (k >= kb) return (eta2 + depth) * mesh->sigma_lcl[k] + eta2;
   }
   
   assert(false && "ZCOOR: invalid level or coordinate system\n");
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

      assert(false && "computeBottomDragVegetation: nchi == 1 not fully implemented. see below");
      // Cdmax=-1 !max. Cd at node for this process (info only)
      // Drag at nodes
/////////////////!$OMP   do reduction(max: iwbl_itmax) reduction(.or.: ltmp)   XXXXXXXXXXXX
// !        z0b_save(:) = 0.d0 ! Initialization of mixing length at bottom (z0b, T. Guerin) 
           
      assert(rough >= 0.0 && "computeBottomDragVegetation: rough < 0");
      if (rough == 0.0) {
         Cd = 0.0;

      } else { //roughness >0 
         
         const double bthick_ori = z[kb + 1] - z[kb];  // thickness of bottom bnd layer
         const double bthick = max(dzb_min, bthick_ori);

         assert(bthick > 0.0 && "computeBottomDragVegetation: bthick <= 0");

         double tmp = 1.0/(2.5*log(bthick/rough));
         Cd = tmp*tmp;

      } // if rough


      // RECOVER THIS: XXX
      // if (it == iths_main + 1) log("Cd min/max at 1st step= ", minval(Cdp), maxval(Cdp));
      // Output warning for WBL if iteration didn't converge


   } else if constexpr(nchi == -1) {
   // Manning 
      Cd = 0.0;

      if (isDry) return;

      const double htot = max(hmin_man, htot2()); 
      Cd = grav2*rmanning*rmanning/pow(htot, 0.333);

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
                                    double *drag, double *roughness, double *rmanning) {

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

// printf("NODE i = %d, node = %p, node.numNodes = %d\n", node.idx, &node, node.numNodes);      

   node.wind_angle = wind_rotate_angle[i];


   node.p.x   = node_coords[i * 3];
   node.p.y   = node_coords[i * 3 + 1];
   node.depth = node_coords[i * 3 + 2];

   if (use_ll) {
   node.lon = node_coords_ll[i * 2    ]*deg2rad;
   node.lat = node_coords_ll[i * 2 + 1]*deg2rad;
   }


   node.numSides = 0;

   // drag, rough, and manning are actually stored in the same memory location. This is just sugar to make the code more readable without using unnecessary memory.
   if      constexpr(nchi ==  0) node.drag     = drag[i];
   else if constexpr(nchi ==  1) node.rough    = roughness[i];
   else if constexpr(nchi == -1) node.rmanning = rmanning[i];
   else assert(false && "populateNodes: selected nchi not implemented.");
   
   for (int k = 0; k < nvrt; k++) node.rho[k] = rho0;

   
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
         
         node.nodeBnd->segment->ifltype = (isbnd[i]==0)?3:0;
         node.nodeBnd->segment->iettype = 3; // XXXXXXXXX

         // Initialize the length and area of the segment to nan to indicate that they are not yet computed
         node.nodeBnd->segment->clen  = nan("1");
         node.nodeBnd->segment->carea = nan("1");

      }
//          else {  XXXXXX
//             node.nodeBnd->segment = (BoundarySegment *)(isbnd[i]); // Hack, isbnd[i] will be -1 or -2 in this case
// printf("EEEEEPP isbnd[i] = %d, node.nodeBnd->segment = %p\n", isbnd[i], node.nodeBnd->segment);            
//          }
      

   }
   else {

      node.nodeBnd = nullptr; // For internal nodes

      if (isbnd[i] == landSegment) {
         node.nodeBnd = &landBoundaryNode;
         node.nodeBnd->segment = &landBoundarySegment;
      }
      // node.nodeBnd->node = &node;
      // node.nodeBnd->segment = (BoundarySegment *)(isbnd[i]); // Hack, isbnd[i] will be -1 or -2 in this case
      // // node.nodeBnd = nullptr;
      // node.nodeBnd->segment = (BoundarySegment *)(isbnd[i]); // Hack, isbnd[i] will be -1 or -2 in this case
      
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

   // Normal, unitary vector to the side
   side.n = { -(side.n1->p.y - side.n2->p.y), (side.n1->p.x - side.n2->p.x)};
   side.n = normalize(side.n);

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

   // Compute the center of the side
   side.center = 0.5*(side.n1->p + side.n2->p);

   assert(ishapiro == 1 && "Only shapiro = 1 is implemented");
   side.shapiro = shapiro0;

#ifdef GAUSSIAN_HILL_VELOCITIES
// XXXX
if (i == 0) printf("GAUSSIAN HILL VELOCITIES!!!\n");
constexpr double rot_per=3000.0;
constexpr double rot_f=0.0*2.0*pi/rot_per;
for (int k = 0; k < nvrt; k++) side.u2[k] = {-side.center.y*rot_f , side.center.x*rot_f};
// XXXXX      
#endif
   // Compute the z level at the center of the side
   for (int k = 0; k < nvrt; k++) side.zc[k] = 0.5*(side.n1->z[k] + side.n2->z[k]);
   
   // Compute the length of the side
   const Vector2 vLen = side.n2->p - side.n1->p;
   side.len = sqrt(dot(vLen, vLen));


   // Compute the Coriolis coefficient

   if constexpr(ncor == -1 || ncor == 0) {
      side.coriolis = coricoef;
   }
   else {
      const double sphi = side.lat();  // Latitude at the center of the side

// XXXXXX ERROR WHEN use_ll is false??????
      if constexpr(ics == 1) {
         // XXXXXXX
         side.coriolis = fc + beta*(sphi - sfea0);  // First order Taylor expansion of the Coriolis coefficient around the reference latitude sfea0
      }
      else {
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

   // Compute the centroid of the element and its depth
   Vector2 center = {0.0, 0.0};
   elem.depth = 0.0;
   for (int j = 0; j < i34[i]; j++) {

// printf("WWW i = %d, node = %p, elem.i34 = %d\n", i, elem.nodes[j], elem.i34);    


      center += elem.nodes[j]->p/i34[i];
      elem.depth  += elem.nodes[j]->depth/i34[i];
   }
   elem.center = center;

   
   // if constexpr(ics==2) elem.zc += elem.znd[j]/i34;   XXX Implement this?
   // if(dp(elnode(j,ie))<dpe(ie)) dpe(ie)=dp(elnode(j,ie))

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
   //  For ics=2, this is done inside element/ll frame XXX
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
         const double dPhi_j_dx = (elem.nodes[id2]->p.y - elem.nodes[id3]->p.y)/(2*area); 
         const double dPhi_j_dy = (elem.nodes[id3]->p.x - elem.nodes[id2]->p.x)/(2*area); 
         elem.gradPhi[j] = {dPhi_j_dx, dPhi_j_dy};
      }
      else  { // Quad; evaluate at centroid
         const double dPhi_j_dx = (elem.nodes[id2]->p.y - elem.nodes[id4]->p.y)/(2*area); 
         const double dPhi_j_dy = (elem.nodes[id4]->p.x - elem.nodes[id2]->p.x)/(2*area); 
         elem.gradPhi[j] = {dPhi_j_dx, dPhi_j_dy};
      } // if Tri

   } // for j   
   
   

}

__device__ void Mesh::populateBoundaryNodes(int *nond, int *iond,
                                            double *emo, double *efa, 
                                            double *umo, double *ufa, 
                                            double *vmo, double *vfa) {


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

   for (int jfr = 0; jfr < nbfr; jfr++) {

      const int idx2 = iseg*mnond*nbfr + i + jfr*mnond; //XXXXX
      nodeBnd->efa[jfr] = efa[idx2];
      nodeBnd->emo[jfr] = emo[idx2];
      nodeBnd->umo[jfr] = {umo[idx2], vmo[idx2]};
      nodeBnd->ufa[jfr] = {ufa[idx2], vfa[idx2]};

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

   // printf("Populating mesh for partition %d\n", partNum);

   // XXXXXX THIS function does almost nothing. 

   nbs = *sharedIdxSideBoundary + 1;
   nbp = *sharedIdxNodeBoundary + 1;

   nEntries = (nodeNeighNodesStart[nbp] - nodeNeighNodesStart[0]) + nbp;

}


extern "C" __global__ void setupCoefficientMatrixAndRHS(double time) {
   // Setup coefficient matrix, sparsem, for the shallow water wave equation
   // No elevation essential b.c. are imposed yet but other b.c. is imposed
   // 
   // This kernel is called for all nodes.


   NODE(node1);

if (node1->numElements > 7) printf("111 idx = %d, node1->numElements = %d\n", node1->idx, node1->numElements);

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

         const Point2 x1 = node1->p;
         const Point2 x2 = node2->p;
         const Point2 x3 = node3->p;

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

         // To clarify the notation in this section: 
         //   i stands for the node under consideration (node1)
         //   j (node2), k (node3), l (node4) are the other nodes of the quad.
         const int id2 = elem->nxq(1 + id1);
         const int id3 = elem->nxq(2 + id1);
         const int id4 = elem->nxq(3 + id1);
         Node *node2 = elem->nodes[id2];
         Node *node3 = elem->nodes[id3];
         Node *node4 = elem->nodes[id4];         

         // integrals of phi_i*phi_j, phi_i*phi_k, phi_i*phi_l
         const double ipp_ii = elem->intPhiPhi(id1, id1);
         const double ipp_ij = elem->intPhiPhi(id1, id2);
         const double ipp_ik = elem->intPhiPhi(id1, id3);
         const double ipp_il = elem->intPhiPhi(id1, id4);

         // integrals of (grad phi_i) * (grad phi_j), (grad phi_i) * (grad phi_k), (grad phi_i) * (grad phi_l)
         const double igpgp_ii = elem->intGradPhiGradPhi(id1, id1);
         const double igpgp_ij = elem->intGradPhiGradPhi(id1, id2);
         const double igpgp_ik = elem->intGradPhiGradPhi(id1, id3);
         const double igpgp_il = elem->intGradPhiGradPhi(id1, id4);

         // Computes the different components of integral I_1 of each node with respect to node i.
         const double I1_ii = ipp_ii + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ii;
         const double I1_ij = ipp_ij + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ij;
         const double I1_ik = ipp_ik + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_ik;
         const double I1_il = ipp_il + (grav3*thetai*thetai*dt*dt*hhatb)*igpgp_il;

         node1->addToSparseMatrix(node1, I1_ii);
         node1->addToSparseMatrix(node2, I1_ij); 
         node1->addToSparseMatrix(node3, I1_ik);
         node1->addToSparseMatrix(node4, I1_il);
  


      
      } 

 
      // Integral I_4
      // I_4 = \int_\Omega \left[ \phi_i\eta^n + \theta\Delta t \nabla\phi_i\cdot\pmb{E} + (1-\theta)\Delta t \nabla\phi_i \cdot\pmb{U}^n \right]d\Omega
      // Where
      // U is pmb{U}^{n+1} = \pmb{E} - g\theta\breve{H}\Delta t\nabla\eta^{n+1}

      // Since \nabla\phi_i is constant over the element, for the integral of U we can compute the mean value
      // of U over the element.
      Vector2 Umean = {0.0, 0.0};
      for (int m=0; m<elem->i34; m++) {      
         const Side *side = elem->sides[m];

         if constexpr(ics == 1) Umean += side->U; // XXXXX IS THIS JUST ONE VALUE OR TWO????
         else                   Umean += side->frame.project_hvec(side->U, elem->frame);
      }
      Umean = Umean/elem->i34;

      const Vector2 gradPhi_1 = elem->gradPhi[id1];

      // Integral of \nabla\phi_i \cdot \pmb{U}^n
      const double I_4_1 = dot(gradPhi_1, Umean);

      // Integral of \nabla\phi_i \cdot \pmb{E}
      const double I_4_2 = dot(gradPhi_1, elem->ghat1);
     
      node1->qel += dt*area*((1 - thetai)*I_4_1 + thetai*I_4_2);      

      // Additional terms 
      double vnorm = 0.0;
      if constexpr(imm == 2)  { // Pre-compute vnorm on the element
         Vector3 ubed;
         // update_bdef(time, elem->center, elem->depth, ubed); XXX

         Vector3 gradP = {0.0, 0.0, 1.0};
         for (int m=0; m < elem->i34; m++) {
           const Vector2 gradP_2D = elem->nodes[m]->depth*elem->gradPhi[m];
           gradP += toVector3(gradP_2D, 0.0);
         }   
         vnorm = dot(ubed, gradP)/sqrt(dot(gradP, gradP));
      } // imm==2

      if (elem->isTri()) {
         for (int l=0; l<3; l++) {
            double fac;

            Node *node = elem->nodes[l];

            fac = (node == node1) ? 2.0 : 1.0;
            if constexpr(imm==2) node1->qel += (1.0/12.0)*elem->area*fac*(node->eta2 + dt*vnorm);
            else                 node1->qel += (1.0/12.0)*elem->area*fac*(node->eta2 + node->bdef2 - node->bdef1);
           
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
   // I_2,6 are not needed if essential b.c. are enforced by elminating rows and columns
   NodeBoundary *nodeBnd = node1->nodeBnd;
   if (nodeBnd != nullptr && nodeBnd != &(mesh->landBoundaryNode))  { 
      // If it is an open boundary node...

      // At least one open bnd side, maybe two
      for (int l=0; l<2; l++) { 

         // Gets each of the sides sharing the node. 
         SideBoundary *sideBnd = nodeBnd->sideBnd(l);

         // if (sideBnd == nullptr) continue; // Skip if no boundary side
         

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
                  for (int k = side->kb; k<nvrt-1; k++) {
                     const double vn1 = dot(sideBnd->ubc[k  ], side->n); // outer normal
                     const double vn2 = dot(sideBnd->ubc[k+1], side->n);

                     Vn += 0.5*(side->zc[k+1] - side->zc[k])*(vn1 + vn2);
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

                  const double etam = 0.5*(node1->eta2 + node2->eta2);

                  const double tmp = polyDischarge->eval(etam)/sideSegment->clen; // XXXX IMPLEMENT Poly. It may need the grade (4 in this case)
                  assert(tmp >= 0.0 && "bar{f} < 0 at discharge boundary");

                  const double con0 = (1.0/6.0)*side->len*tmp*thetai*dt;
                  node1->addToSparseMatrix(node1, con0*2.0);
                  node1->addToSparseMatrix(node2, con0);

               } // if Flather or Discharge

// printf("node1 = %p, thetai = %f, dt = %f, ri3 = %f\n", node1, thetai, dt, ri3);
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
   __global__ void setupCGSolver(double* data, int* indices, int* indptr, double* b) {
      int off = 0;
      for (int i = 0; i < mesh->np; i++) {

         const Node *node = mesh->localNodeData + i;
         indptr[i] = off;

         // Diagonal first
         data[off]    = node->coeffDiag;
         indices[off] = i;
         ++off;
         // printf("i = %d, off = %d, node->nodes[i] = %p\n", i, off, node->nodes[i]);
         // Off-diagonals
         for (int k = 0; k < node->numNodes; ++k) {
            const Node* nb = node->nodes[k];
            // printf("i = %d, k = %d, nb = %p\n", i, k, nb);
            // Map pointer -> column index by pointer arithmetic
            int j = nb->idx;
            // (optionally: if (j == i) data[diag_off] += nd.coeffMatrix[k]; else ...)
            data[off]    = node->coeffMatrix[k];
            indices[off] = j;
            ++off;
            // printf("i = %d, j = %d, off = %d, data[off] = %f\n", i, j, off, data[off]);
         }

         // RHS
         b[i] = node->qel;
      }

      indptr[mesh->np] = off;
   }
   
}

extern "C" {
   __global__ void setupCGSolver2(double* diag, double* data, int* indices, int* indptr, double* b) {
      NODE(node);
      const int idx = node->idx;
      // Computes the index in the nzz array

      const int nzzIdx = (node->nodes - mesh->nodeNeighNodes) + idx;

      indptr[idx] = nzzIdx;
      data[nzzIdx]    = node->coeffDiag;
      indices[nzzIdx] = idx;

      for (int k = 0; k < node->numNodes; ++k) {
         const Node* nb = node->nodes[k];

         int j = nb->idx;
         // (optionally: if (j == i) data[diag_off] += nd.coeffMatrix[k]; else ...)
         data[nzzIdx + k + 1]    = node->coeffMatrix[k];
         indices[nzzIdx + k + 1] = j;
      }

      diag[idx] = 1.0/node->coeffDiag;
      b[idx] = node->qel;
   }
}


extern "C" {
   __global__ void setupCGSolver3(float* diag, float* data, int* indices, int* indptr, float* b) {
      NODE(node);
      const int idx = node->idx;
      // Computes the index in the nzz array

      const int nzzIdx = (node->nodes - mesh->nodeNeighNodes) + idx;

      indptr[idx] = nzzIdx;
      data[nzzIdx]    = (float)node->coeffDiag;
      indices[nzzIdx] = idx;

      for (int k = 0; k < node->numNodes; ++k) {
         const Node* nb = node->nodes[k];

         int j = nb->idx;
         // (optionally: if (j == i) data[diag_off] += nd.coeffMatrix[k]; else ...)
         data[nzzIdx + k + 1]    = (float)node->coeffMatrix[k];
         indices[nzzIdx + k + 1] = j;
      }

      diag[idx] = (float)(1.0/node->coeffDiag);
      b[idx] = (float)node->qel;
   }
}

#ifdef popo
extern "C" {
   __global__ void setupCGSolver4(double* diag_d, double* data_d, double* b_d, 
                                  float*  diag_f, float*  data_f, float*  b_f, 
                                  int* indices, int* indptr) {

      NODE(node);
      const int idx = node->idx;
      // Computes the index in the nzz array

      const int nzzIdx = (node->nodes - mesh->nodeNeighNodes) + idx;

      indptr[idx] = nzzIdx;
      data[nzzIdx]    = node->coeffDiag;
      indices[nzzIdx] = idx;

      for (int k = 0; k < node->numNodes; ++k) {
         const Node* nb = node->nodes[k];

         int j = nb->idx;
         // (optionally: if (j == i) data[diag_off] += nd.coeffMatrix[k]; else ...)
         data[nzzIdx + k + 1]    = node->coeffMatrix[k];
         indices[nzzIdx + k + 1] = j;
      }

      diag[idx] = 1.0/node->coeffDiag;
      b[idx] = node->qel;
   }
}


extern "C" {
   __global__ void setupCGSolver3(float* diag, float* data, int* indices, int* indptr, float* b) {
      NODE(node);
      const int idx = node->idx;
      // Computes the index in the nzz array

      const int nzzIdx = (node->nodes - mesh->nodeNeighNodes) + idx;

      indptr[idx] = nzzIdx;
      data[nzzIdx]    = (float)node->coeffDiag;
      indices[nzzIdx] = idx;

      for (int k = 0; k < node->numNodes; ++k) {
         const Node* nb = node->nodes[k];

         int j = nb->idx;
         // (optionally: if (j == i) data[diag_off] += nd.coeffMatrix[k]; else ...)
         data[nzzIdx + k + 1]    = (float)node->coeffMatrix[k];
         indices[nzzIdx + k + 1] = j;
      }

      diag[idx] = (float)(1.0/node->coeffDiag);
      b[idx] = (float)node->qel;
   }
}
#endif


extern "C" {
   __global__ void solverCopyResult(double* x) {
      NODE(node);
      const int idx = node->idx;

      node->mX = x[idx];
   }
}

extern "C" {
   __global__ void rollUpElevation(void) {

      NODE(node);

      node->eta1 = node->eta2;

   }


   __global__ void backtrackAllNodes(void) {

      // max/min # of sub-divisions
      const int ndelt_max = max(1.0, dt/dtb_min);
      const int ndelt_min = max(1.0, dt/dtb_max); 

      // Removed all references to swild98 XXX
      // swild98=0.d0 !init


      SIDE(side);

      const int kb = side->kb;


      // Initialize btrack velocities to the actual velocities at the side 
      for (int k = 0; k < nvrt; k++) side->btu[k] = side->u2[k];

      if (side->isDry) return;


      // Checks that at least one element is not dry
      const Element *elem1 = side->elems[0];
      const Element *elem2 = side->elems[1];


      assert ((elem1 != nullptr && !(elem1->isDry)) || (elem2 != nullptr && !(elem2->isDry)) && "backtrackAllNodes: non dry side without two dry elements.");
      
      // XXXXXX
      // !swild_tmp to store global coord. of the starting side
      // swild_tmp(1)=xcj(isd0); swild_tmp(2)=ycj(isd0); swild_tmp(3)=zcj(isd0)
      // !swild10_tmp to store frame at starting pt for ics=2 (not used for ics=1)
      // !Use sframe2 
      // swild10_tmp(1:3,1:3)=sframe2(:,:,i) !pframe(:,:,isidenode(1,isd0)) 


      for (int k = kb; k < nvrt; k++) {
      //  Initialize (xt,yt,zt), nnel and vel.
      //  For ics=2, the coordinates are in local frames
      //  Caution! nnel must be initialized inside this loop as it is updated inside.


         const int iadvf = side->iadv();  // Advection flag

         Point3 xt;
         if constexpr(ics==1) {
            xt = side->center3D(k); 
         }
         else // lat/lon; in side lat/lon frame
         {
            assert(false && "lat/lon not implemented in backtrackAllNodes");
         // XXX IMPLEMENT THIS
      //          xt=0.d0
      //          yt=0.d0
      //          !centroid coord. for nudging
      //          call project_pt('g2l',xctr(ie0),yctr(ie0),zctr(ie0), &
      //   &(/xcj(isd0),ycj(isd0),zcj(isd0)/),swild10_tmp,xctr2,yctr2,tmp)
         }  // ics



         // If on land or open boundary...
         const bool ifl_bnd =  (side->sideBnd != nullptr);
         
         const double vmag = length(side->u2[k]);

         if (vmag <= velmin_btrack) {
            // No activity 

            side->btu[k] = side->u2[k];
            
         } else {  // do btrack

            // Compute number of sub-divisions based on local gradients 
            // XXX IS THIS NECESSARY?
            double suma = 0.0;
            int icount = 0;
            for (int ii=0; ii<2; ii++) {

               const Element *elem = side->elems[ii];
                 
               if (elem == nullptr) continue;

               icount++;

               // Not strictly along z; in ll frame for ics=2

               double utmp[4];
               double vtmp[4];
               for (int jj = 0; jj < elem->i34; jj++) {

                  const Node *node = elem->nodes[jj];

                  if constexpr(ics==1) {
                     utmp[jj] = node->u2[k].x;
                     vtmp[jj] = node->u2[k].y;  
                  } else {
                     assert(false && "lat/lon not implemented in backtrackAllNodes");
                     // call project_hvec(uu2(j,nd0),vv2(j,nd0),pframe(:,:,nd0),eframe(:,:,ie),utmp0(jj),vtmp0(jj))
                  }
               } // for jj

               const Vector2 grad_u = elem->grad2(utmp); 
               const Vector2 grad_v = elem->grad2(vtmp);

               suma += dt*sqrt(dot(grad_u, grad_u) + dot(grad_v, grad_v));
            } // for ii

            assert(icount > 0 && "Impossible 77");

            // '4' is somewhat arbitrary
            const int ndelt = clamp(int(suma/icount)*4, ndelt_min, ndelt_max); 

            const double dtbk = dt/ndelt; // Target btrack step; may be smaller sometimes
         

            Vector3 u_int = {0.0, 0.0, 0.0};

            side->btrack(k, dt, dtbk, u_int);
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
        node->grav2=grav;

      }

   }


   __global__ void computeTidalPotential(double time, double* tfreq, double* tamp, double* tnf, double* tear, int* jspc) {

      NODE_BND(nodeBnd);

      Node *node = nodeBnd->node;

      node->ETP = 0.0;

      for (int j=0; j<ntip; j++) {

         const double ncyc = int(tfreq[j]*time/(2.0*pi));
         const double arg = tfreq[j]*time - ncyc*2.0*pi + jspc[j]*node->lon + tear[j];
         node->ETP += 0.69*ramp*tamp[j]*tnf[j]*node->fun_lat[jspc[j]]*cos(arg);

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

      assert(ics == 1 && "computeHorizontalViscosity: implemented only for ics=1");

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
               const Side *side1 = elem->sides[elem->nxq(idxSide + 1)];
               const Side *side2 = elem->sides[elem->nxq(idxSide - 1)];

               sumU += side1->u2[k] + side2->u2[k] - 2*side->u2[k];   // XXXX Is this correct?

            }

            side->nuGradU[k] = hvis_coef0*sumU/dt;
   
         }
      }
      else if constexpr(ihorcon == 2) {
         assert(false && "computeHorizontalViscosity: implemented only for ihorcon == 1");
      }


   }


   __global__ void setEssentialBC(void) {
      // Set essential boundary conditions 
      NODE(node); // Each boundary node
      NodeBoundary *nodeBnd = node->nodeBnd;
      
// nodeBnd->isEssentialBC = true;      
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

      const Node *node1 = side->n1;
      const Node *node2 = side->n2;

      // Vertically integrated velocity U (in ll if ics=2)
      for (int k = side->kb; k < nvrt - 1; k++) { 
         side->U += side->Dz(k + 1)*0.5*(side->u2[k] + side->u2[k + 1]);  
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
      const double vmag1 = length(side->btu[kb + 1]); //BACKTRACKED VELOCITY
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
         if (ihhat == 1) side->hhat = max(0.0, side->hhat);
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

      int k1 = max(kb, nvrt / 2);
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
            // assert(zrat >= 0.0 && zrat <= 1.0 && "computeSideDensityGradient: zrat is not in [0, 1]");
printf("###### k1 = %d, side->idx = %d, zc[%d] = %.16f, z[%d] = %.16f, z[%d] = %.16f, z[%d] = %.16f, zrat = %.16f\n",
      k1, side->idx, k, side->zc[k], k1, e1->z[k1], k1 + 1, e1->z[k1 + 1], e1->gradRho[k1].x, e1->gradRho[k1].y, zrat);            
            gradRho += (e1->gradRho[k1]*(1.0 - zrat) + e1->gradRho[k1 + 1]*zrat)*w1;
         }      
         if (e2 != nullptr && k > kb) {
            k2 = e2->findLevel(side->zc[k], k2); // finds the index at the half level starting the search at k2
            const double zrat = (side->zc[k] - e2->z[k2])/(e2->z[k2 + 1] - e2->z[k2]);
            
            gradRho += (e2->gradRho[k2]*(1.0 - zrat) + e2->gradRho[k2 + 1]*zrat)*w2;
         }

         if constexpr(ics==2) {
            assert(false && "Not implemented: ics == 2");
            // IMPLEMENT THIS XXXXXX
            // call project_hvec(swild(1),swild(2),eframe(:,:,ie),sframe2(:,:,i),tmp1,tmp2)
            // swild(1)=tmp1
            // swild(2)=tmp2
         }
         
         // Integrate (using a central Gauss point) the density gradient from the surface to the bottom 
         // multiplied by -rampbc*(g/rho0)  (In ll if ics=2)
         // CHECK INDICES!!! XXX
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

// printf("side->isDry: %d\n", side->isDry);
      if (side->isDry) return;

      // Computes the gradient of the variables at the side center.
      auto isWet1 = [](const Element* e) { return !(e->isDry); };
      auto isWet2 = [](const Element* e) { return !(e->isDry2()); };
      auto isWet3 = [](const Element* e) { return !(e->depth >= tip_dp); };

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
         for (int k=0; k<nvrt; k++) side->u2[k] = {0.0, 0.0};         

         return;
      }

      const Node *n1 = side->n1;
      const Node *n2 = side->n2;

      const int kb = side->kb;

      const double veg_h = side->veg_h();
      const double grav3 = side->grav3();

      if (side->is2D()) {
         // Warning: Do not use htot2, which is updated.
         const double htot = side->htot1();
         const double hhat = side->hhat;

         // assert(hhat > 0.0 && htot > 0.0 && "solveMomentumEquations: Impossible");

         Vector2 tau2 = {0.0, 0.0};
          
         if constexpr(ics==1) {
            tau2 = side->tau();

         } else {
            assert(false && "Not implemented: ics == 2");
            // call project_hvec(tau(1,node1),tau(2,node1),pframe(:,:,node1),sframe2(:,:,j),vn1,vn2)
            // call project_hvec(tau(1,node2),tau(2,node2),pframe(:,:,node2),sframe2(:,:,j),tt1,ss1)
            // taux2=(vn1+tt1)*0.5d0
            // tauy2=(vn2+ss1)*0.5d0
         } // if ics         

         // hat_gam has units of m/s

         const Vector2 coriolisTerm = side->coriolis*rot90(side->u2[nvrt-1]);

         Vector2 hat_gam = side->btu[nvrt-1] + dt*(coriolisTerm - side->gradPr/rho0 + grav3*side->gradETP +
                           side->bcc[kb] + tau2/htot - grav3*((1 - thetai)*side->gradEta1 + thetai*side->gradEta2));

         hat_gam += dt*side->nuGradU[nvrt-1];  // XXXX compute nuGradU


         const double ux = clamp(hat_gam.x*hhat/htot, -rmaxvelx, rmaxvelx);
         const double uy = clamp(hat_gam.y*hhat/htot, -rmaxvely, rmaxvely);

         for (int k = 0; k < nvrt; k++) side->u2[k] = {ux, uy};

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
            const double dfz = 0.5*(side->dfv(k - 1) + side->dfv(k)); ///XXXXX
 
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
            rrhs[kin] = - 0.5*dz*dt*grad2D;

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

            } else { // k == nvrt - 1 SURFACE
               Vector2 tau2 = {0.0, 0.0};
               if constexpr(ics==1) {
                  tau2 = side->tau();
      
               } else {
                  assert(false && "Not implemented: ics == 2");
                  // call project_hvec(tau(1,node1),tau(2,node1),pframe(:,:,node1),sframe2(:,:,j),vn1,vn2)
                  // call project_hvec(tau(1,node2),tau(2,node2),pframe(:,:,node2),sframe2(:,:,j),tt1,ss1)
                  // taux2=(vn1+tt1)*0.5d0
                  // tauy2=(vn2+ss1)*0.5d0
               } // if ics

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

            } // if k>kbs + 1

            // Baroclinic
            Vector2 *bcc = side->bcc;
            if constexpr(ibc == 0) {
               if (k < nvrt  ) rrhs[kin] += (dz_1*dt/6.0)*(2.0*bcc[k] + bcc[k+1]);
               if (k > kb + 1) rrhs[kin] += (dz_0*dt/6.0)*(2.0*bcc[k] + bcc[k-1]);

            } // if ibc==0

         } // k=kb + 1,nvrt


         tridag_sch<2>(nvrt - kb - 1, alow, bdia, cupp, (double *)rrhs, (double *)soln, (double *)gam);

         for (int k = kb + 1; k < nvrt; k++) {  // YYYYYYYYY can be kb 
            const int kin = k - kb - 1;

            // Velocity is always clamped to maximum values
            side->u2[k] = make_double2(clamp(soln[kin].x, -rmaxvelx, rmaxvelx), clamp(soln[kin].y, -rmaxvely, rmaxvely));
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
      for (int k = 0; k <= kb - 1; k++) side->u2[k] = {0.0, 0.0};


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
  
         for (int k = kb; k < nvrt; k++) { // YYYYYYYYY can be kb + 1
            
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

               const Side *nSide = side->neigh(j); // Every internal side has 4 neighbors

               int kin = k;
               if (!(nSide->isDry)) kin = max(k, kb + 1);
               
               if constexpr(ics==1) {
                  suru += nSide->u2[kin];

               } else {
                  assert(false && "Not implemented: ics == 2");
                  // Vector2 vn1, vn2;
                  // project_hvec(side->u2[kin], side->v2[kin], nSide->sframe, side->sframe, vn1, vn2);
                  // suru=suru+vn1;
                  // surv=surv+vn2;
               }
            }


            side->bcc[k] = side->u2[k] + side->shapiro*(0.25*suru - side->u2[k]); // sframe2 if ics=2

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
         side->u2[kb    ] = side->u2 [nvrt - 1];

      } else {
         // 3D
         for (int k = kb + 1; k < nvrt; k++) side->u2[k] = side->bcc[k];
         
      }

   }


   __global__ void ul0UpdateZCoords(void) {
   // For wet nodes updates the nodes z-coordinates and bottom index, kb.
      
      NODE(node);      


      if (node->isDry2()) { 

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

   __global__ void ul0ComputeDryElementsAndNodes(void)
   // Updates isDry based on eta2 for the nodes. Also computes the element's kb and z.
   {
      ELEM(elem);

      elem->isDry = elem->isDry2(); // Element is dry if any of its nodes are dry (based on eta2, not on isDry)

      // All the nodes of a wet element are wet.
      if (!elem->isDry) {
         elem->nodes[0]->isDry = false;
         elem->nodes[1]->isDry = false;
         elem->nodes[2]->isDry = false;
         if (elem->isQuad()) elem->nodes[3]->isDry = false;

      }


      elem->kb = 0;

      if (elem->isDry) return;

      // Wet
      const Node *n1 = elem->nodes[0];
      const Node *n2 = elem->nodes[1];
      const Node *n3 = elem->nodes[2];  
      const Node *n4 = elem->nodes[3]; // nullptr for triangles

   //    if(elem->isDry2()) { // XXXX
   //        write(errmsg,*)'level0: Element-node inconsistency (0):',ielg(i),idry_e2(i), &
   //   &iplg(elnode(1:i34(i),i)),idry2(elnode(1:i34(i),i)),idry(elnode(1:i34(i),i))
   //        call parallel_abort(errmsg)
   //    }


      // Computes element's kb as the minimum of the nodes' kb.
      const int kb1 = n1->kb;
      const int kb2 = n2->kb;
      const int kb3 = n3->kb;
      const int kb4 = (elem->isQuad()) ? n4->kb : INT_MAX;


      elem->kb = min(min(kb1, kb2), min(kb3, kb4));


      // Computes the element z as the mean of the nodes z.
      for (int k = elem->kb; k < nvrt; k++) {

         // Compute the element z as the mean of the nodes z.
         elem->z[k] = (n1->z[max(k, kb1)] + n2->z[max(k, kb2)] + n3->z[max(k, kb3)]);
         if (elem->isQuad()) elem->z[k] += n4->z[max(k, kb4)];          
         elem->z[k] /= elem->i34;

         assert(((k < elem->kb + 1)  || (elem->z[k] - elem->z[k-1] > 0.0)) && "Weird element z");
      }
                  
   }

   __global__ void ul0UpdateRewettedNodes(int it, int iths)
   {
      NODE(node);
      // Compute velocity, S, T for re-wetted nodes (q2 and xl are fine)

      if (it == iths) return;


      if(!node->isDry && node->tmp.b) { // rewetted
         for (int k=0; k < nvrt; k++) {

            Vector3 utmp={0.0, 0.0, 0.0};

            double ttmp=0.0;
            double stmp=0.0;
            int icount=0;
            for (int j = 0; j < node->numNodes; j++) {


               // Wet nbrs not affected by this part and so each sub-domain should use same values
               if (!node->isDry) { // all indices extended
                  icount++;

                  // Assume small element size in wet/dry zone so pframes are close to each other
                  utmp += node->uvw2(k);
                  
                  ttmp += node->tracers[k].tracer[0]; 
                  stmp += node->tracers[k].tracer[1]; 
               }
            }

            if(icount > 0) {
               const double inv_icount = 1.0/icount;
               node->u2[k]     = {inv_icount*utmp.x, inv_icount*utmp.y};
               node->w2[k]     = inv_icount*utmp.z;
               node->tracers[k].tracer[0]   = inv_icount*ttmp;
               node->tracers[k].tracer[1]   = inv_icount*stmp;
            }
         }
      }
   }


   __global__ void ul0UpdateSides(int it, int iths) {

      SIDE(side);

      // A side is wet if and only if at least one of its elements is wet.
      const bool isWet = !side->elems[0]->isDry || (side->isInternal() && !side->elems[1]->isDry);
 
      const Node *n1 = side->n1;
      const Node *n2 = side->n2;
     
      const int kb1 = n1->kb;
      const int kb2 = n2->kb;

      side->kb = 0; // dry

      if (isWet) {

         assert(~(n1->isDry || n2->isDry) && "Side-node inconsistency. At least one node of a wet side is dry.");

         assert(side->htot2() > h0 && "Weird side. Total depth is smaller than h0.");

         side->kb = min(kb1, kb2);

         for (int k = 0; k < side->kb; k++) side->zc[k] = 0.0;  // XXXX NECESSARY???

         for (int k = side->kb; k < nvrt; k++) {
            side->zc[k] = 0.5*(n1->z[max(k, kb1)] + n2->z[max(k, kb2)]);

         
            assert(((k < side->kb + 1) || (side->zc[k] - side->zc[k-1] > 0.0)) && "Weird side (3)");

            // Compute velocity for re-wetted sides 
            if ((it != iths) && (side->isDry && isWet)) {
            
               Vector2 utmp={0.0, 0.0};

               int icount=0;

               // For all the sides on the two elements that share this side
               for (int j=0; j < 2; j++) {
                  const Element *elem = side->elems[j];
                  if (elem == nullptr) continue;

                  for (int jj=0; jj < elem->i34; jj++) {
                     const Side *nSide = elem->sides[jj];
                     if (nSide->isDry) continue;

                     icount++;
                     utmp += nSide->u2[k];                     
                  }
               }
                  
               assert(icount > 0 && "No wet sides around this side");
               if (icount > 0)  side->u2[k] = utmp/icount;
                  
            } // if ((it != iths) && (side->isDry && isWet))

         } // for k

      } // if isWet
   
      side->isDry = !isWet;

   }


   __global__ void zeroNodalVelocities(void) {
      NODE(node);

      for (int k = 0; k < nvrt; k++) {
         node->u2[k]   = {0.0, 0.0};
         node->w2[k]   = 0.0;
         node->utmp[k] = {0.0, 0.0};
         node->tmp.d = 0.0;  // reset the temporary real value. It will be used to compute the area-averaged vertical velocity
      }
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
   }


   __global__ void computeNodalHorizontalVelocities1(void) {
      assert(indvel == 0 && "computeNodalVelocities only works for indvel <= 0. IMPLEMENT THIS FOR INDVEL = 1");
      assert(ics == 1 && "computeNodalVelocities only works for ics = 1. IMPLEMENT THIS FOR other ics");

      ELEM(elem);

      const int kb = elem->kb;
      const bool isTri = elem->isTri();

      if constexpr(indvel <= 0) {
         // Assuming discontinoius hvel.


         if (isTri) { // For triangles
         // Interpolates the side velocities at the nodes. A magnitude defined at the sides and interpolated (extrapolated) at one node can be
         // computed as the sum of the adjacent side magnitudes minus the side magnitude that is not adjacent to the node.

            for (int j = 0; j < 3; j++) {
               // For every node of the triangle...
               Node *node = elem->nodes[j];

               const Side *s1 = elem->sides[j];
               const Side *s2 = elem->sides[(j + 1) % 3];
               const Side *s3 = elem->sides[(j + 2) % 3];

               const int numElements = node->numElements;
               const double inv_numElements = 1.0/numElements;

               for (int k = kb; k < nvrt; k++)  {
                  Vector2 u = s2->u2[k] + s3->u2[k] - s1->u2[k]; // u interpolated at the node   
                  u.x = clamp(u.x, -rmaxvelx, rmaxvelx);
                  u.y = clamp(u.y, -rmaxvely, rmaxvely);
                  // node->utmp[k] = node->utmp[k] + (1.0/numElements) * u; 
                  atomicAdd(&(node->utmp[k].x), inv_numElements * u.x);
                  atomicAdd(&(node->utmp[k].y), inv_numElements * u.y);

               }

            }
         }
         else { // For quads
            assert(false && "computeNodalHorizontalVelocities1: not implemented for quads");

         }

      }
      
   }


   __global__ void computeNodalHorizontalVelocities2(void) {
      NODE(node);

      for (int k = 0; k < nvrt; k++) {
         node->u2[k] = node->utmp[k];
         node->w2[k] = node->w2[k];
      }
      
   }   

   __global__ void computeNodalHorizontalVelocitiesOLD(void) {
   // XXXX I'm not sure where I got this function from. 
   // Compute node's horizonal velocities from side velocities at WHOLE levels.
   // It goes for each side and adds to the nodes the velocity divided by the number of sides.
   // It does it for all nodes, including boundary nodes. Afterwards, the open boundary nodes are set to a different value.
   // It uses the utmp array as a temporary array to store the velocities.

      assert(indvel == 0 && "computeNodalVelocities only works with indvel = 0. IMPLEMENT THIS FOR INDVEL <= 0");

      SIDE(side);

      const int kb = side->kb;

      Node *n1 = side->n1;
      Node *n2 = side->n2;

      for (int k = kb; k < nvrt; k++) {

         // Cast the u2 vector into a Vector2 
         Vector2& n1_u2 = *reinterpret_cast<Vector2*>(&(n1->utmp[k]));
         Vector2& n2_u2 = *reinterpret_cast<Vector2*>(&(n2->utmp[k]));

         n1_u2 = side->u2[k]/n1->numSides;
         n2_u2 = side->u2[k]/n2->numSides;

      }


   }

   __global__ void resetVelocityOpenBndNodes(void) {
   // Set velocity to zero for open boundary nodes where the velocity is imposed.

      NODE_BND(nodebnd);

      if (nodebnd->segment != nullptr && nodebnd->segment->ifltype != 0) {
         for (int k = 0; k < nvrt; k++) { nodebnd->node->u2[k] = {0.0, 0.0}; nodebnd->node->w2[k] = 0.0; }
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

      if (elem->isDry) return;

      const int kb = elem->kb;
      const bool isTri = elem->isTri();

      assert(kb >= 0 && "computeVerticalVelocity: wrong kb");

      // Compute the average (over the element) of the bottom deformation at the previous and current timesteps
      const double av_bdef1 = elem->mean(&Node::bdef1);
      const double av_bdef2 = elem->mean(&Node::bdef2);

      Vector3 elem_n1 = elem->n(kb); // normal vector to the bottom face of the "volumetric" element (3D)

      double area3D;

      // Bottom face flux
      double bottomFlux = 0.0;
      Vector2 u2meanBottom = {0.0, 0.0};
      const Vector2 gradh = elem->grad2(&Node::depth);
      if constexpr(imm==2) {
         Vector3 ubed;

         assert(false && "computeVerticalVelocity: not implemented for imm==2");
         // updateBDEF(ubed);
     
         // const double bflux0 = dot(ubed, elem_n1);

         // elem->w[kb] += ubed.z;
      }
      else {
         Vector2 umean;
         if constexpr(ics == 1) {
            u2meanBottom = elem->mean(kb, &Side::u2);
         }
         else {
            assert(false && "computeVerticalVelocity: ics==1 not implemented for imm==0 or 1");
         }

         // Error: /=0 for 2D (but OK b/cos fluxes are 0 below for transport)
         elem->w[kb] = (av_bdef2 - av_bdef1)/dt - dot(gradh, u2meanBottom); // XXXX Is the dot prod always 0?

         bottomFlux = - dot(gradh, u2meanBottom);
      } // if imm

      // Compute w for kb+1 to nvrt-1
      for (int l = kb; l < nvrt - 1; l++) {

         // normal vector to the top face of the "volumetric" element (3D)
         Vector3 elem_n2 = elem->n(l + 1, area3D);  // Get the normal vector and the area of the top face (as a 3D surface)

         // Compute the total flux at the sides
         double sideFlux = 0.0;
         Vector2 u2top = {0.0, 0.0};
         for (int j = 0; j < elem->i34; j++) {
            const Side *side = elem->sides[j];

            const double sideArea = side->area(l);
            const Vector2 u2c = 0.5*(side->u2[l] + side->u2[l + 1]); // velocity at the center of the side
            const double sign = elem->ssign[j]; 

            sideFlux += sign*dot(u2c, side->n)*sideArea;
         }
         
         // Compute the horizontal velocity at the center of the top face as an average of the 
         //velocities at the centers of the sides
         const Vector2 u2meanTop = elem->mean(l + 1, &Side::u2);

         // The flx at the top ace can be computed from the other flows.
         const double topFlux = - sideFlux - bottomFlux;
         const double topHorizontalFlux = dot(u2meanTop, toVector2(elem_n2))*area3D; // Part of the flux related to the horizontal velocity

         elem->w[l + 1] = (-sideFlux - topHorizontalFlux - bottomFlux)/(elem_n2.z*area3D);

         elem_n1 = elem_n2;

         bottomFlux = - topFlux;
      }
   }  






//    __device__ void updateLevels1(void) {
//       // Routine to update level indices and wetting and drying.
//       // Used when resolution is fine enough.
//       // ONLY WORKS WITH PURE TRIs


//       ELEM(elem);

//       const int i34 = elem->i34;

//       // Flag for comm timing
//       const bool cwtime = (it != iths);

//       // An element is wet if and only if depths at all nodes >h0 
//       // A node is wet if and only if at least one surrounding element is wet
//       // A side is wet if and only if at least one surrounding element is wet

//       // Initialize element flags for first step

//       if (it == iths) {
//          elem->isDry = false;

//          for (int j = 0; j < i34; j++) {
//             // For all nodes in the element...

//             const Node *node = elem->node(j);
//             const double htot = node->htot();
               
//             // For an element to be dry, at least one node must be dry
//             if (htot <= h0) elem->isDry = true;
//             break;
//          }
//       }

//       //  Wetting/drying algorithm
//       const double isDry2 = elem->isDry; // starting from step n's indices

//       if (it != iths) {

//          // XXXX removed a fragment of code that that recomputed elem->isDry2, but didn't do anything different

//          // istop: 1- ready for final extrap. stage; 2- ready for final
//          // checks and exit loop15
//          int istop = 0;
//          int itr = 0;
//          while (true) {
//             itr = itr + 1;
//             assert(itr <= 100 && "LEVELS1: Too many iterations in wet/dry");

//             int icolor2 = 0; //interface sides

          
//             const Element *elem1 = side->elem[1];
//             const Element *elem2 = side->elem[2];
//             if (elem2 == nullptr && elem1->isDry2 + elem2->isDry2)) {
//                side->icolor2 = 1;
//             }
          


//             // Aug. shoreline sides (must be internal sides)
//             int nsdf = 0;
//             if (side->icolor2==1) {
//                nsdf=nsdf+1
//                isdf(nsdf)=i
//             }
          

//          REDUCE 
//          call mpi_allreduce(nsdf,nsdf_gb,1,itype,MPI_SUM,comm,ierr)
//          if (nsdf_gb==0) return // all wet

//          // Final extrapolation
//          srwt_xchng(1) = (istop == 1)
//          call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
//          if(srwt_xchng_gb(1)) then // all ranks ready
//             if(myrank==0) write(16,*)'doing final extrapolation in levels1...'

//             icolor=0 // frontier nodes for extrapolation
//             inew=0 // for initializing and counting su2 sv2
//             for (int i=0; i<nsdf; i++) {
//                const int isd = isdf[i];

//                if (isdel(1,isd)==0.or.isdel(2,isd)==0) then
//                   write(errmsg,*) 'LEVELS1: bnd side (2):',isdel(:,isd),iplg(isidenode(1:2,isd))
//                   call parallel_abort(errmsg)
//                endif


//                // Try to find a dry elem (to take care of some odd cases where
//                // nodeA is interface btw sub-domains)
//                for (int m = 0; m <= 1; m++) {
                  
//                   Element *elem = side->elem[m];
//                   if (elem != nullptr) {
//                      if (elem->isDry2) {

//                         break;
//                      }
//                   }

//                } //m

//                if(elem == nullptr) return;

//                Node *n1 = side->node[0];
//                Node *n2 = side->node[1];

//                // TRICK: finds the node of elem that is not n1 or n2. To do so, adds the pointer values of all the nodes (as ints)
//                // and subtracts the pointer values of n1 and n2.
//                Node *nodeA = (Node *)((int)elem->node[0] + (int)elem->node[1] + (int)elem->node[2] - (int)n1 - (int)n2);

//                if(nodeA->icolor == 1) return;  //this node is done

//                nodeA->icolor = 1; // this node will be done


//                inun=0 // inundation flag
//                for (int j=0; j<nodeA->nne; j++) {
//                   Element *elem2 = nodeA->elem[j];
//                   Side *side2 = nodeA->side[j];

//                   // int id = nodeA->elem[j]->iself;
//                   // int isd2 = nodeA->elem[j]->elside;
//                   if(elem2->icolor2 == 1) {

//                      double tmp = dot(side2->u2[isd2], side2->n);
//                      double flux_t = -tmp*ssign(id,ie2) // inward normal
//                      if(flux_t > 0.0) {

//                      const double etm = max(side2->n1->eta2(), side2->n2->eta2());
//                      if (etm + dp(nodeA) > h0) {}
//                      WHY depth at node????
//                         inun=1;
//                         exit;
//                      }
//                   } // flux_t>0

//                } //icolor2(isd2)==1
               
//             } //j

//             if(inun==1) then
//                eta2(nodeA) = max(eta2(nodeA),-dp(nodeA)+2._rkind*h0)

//                 do j=1,nne(nodeA)
//                   ie2=indel(j,nodeA)
//                   id=iself(j,nodeA)
//                   isd2=elside(id,ie2)
//                   if(icolor2(isd2)==1) then
//                     do l=1,3
//                       nd=elnode(l,ie2)
//                       if(eta2(nd)+dp(nd)<=h0) then 
//                         write(errmsg,*)'LEVELS1: Failed to wet element:',ielg(ie2),iplg(nodeA)
//                         call parallel_abort(errmsg)
//                       endif
//                     enddo !l=1,3
//                     idry_e2(ie2)=0
//                     do l=1,2 !sides sharing nodeA
//                       id1=elside(nx(id,l),ie2)

//                       swild2(1,1:nvrt)=su2(1:nvrt,isd2)
//                       swild2(2,1:nvrt)=sv2(1:nvrt,isd2)

//                       if(inew(id1)==0) then
//                         su2(1:nvrt,id1)=swild2(1,1:nvrt)
//                         sv2(1:nvrt,id1)=swild2(2,1:nvrt)
//                         inew(id1)=1
//                       else
//                         su2(1:nvrt,id1)=su2(1:nvrt,id1)+swild2(1,1:nvrt)
//                         sv2(1:nvrt,id1)=sv2(1:nvrt,id1)+swild2(2,1:nvrt)
//                         inew(id1)=inew(id1)+1
//                       endif
//                     enddo !l=1,2
//                   endif !icolor2(isd2)==1
//                 enddo !j=1,nne(nodeA)
//               endif !inun==1
//             enddo !i=1,nsdf

//             call exchange_e2di(idry_e2)
//             call exchange_p2d(eta2)


//             ltmp=.false. // flag for wetting occurring

//             do i=1,ns
//               ltmp=ltmp.or.inew(i)/=0
//               if(inew(i)/=0) then

//                 su2(1:nvrt,i)=su2(1:nvrt,i)/dble(inew(i))
//                 sv2(1:nvrt,i)=sv2(1:nvrt,i)/dble(inew(i))
//               endif
//             enddo !i

//             srwt_xchng(1)=ltmp

//             istop=2

//           endif // srwt_xchng_gb; final extrapolation

//           if(istop/=2) then

//             istop=1 // stop iteration and go to extrapolation stage; initialize first
//             do i=1,nsdf // aug.
//               isd=isdf(i)
//               do j=1,2
//                 nd=isidenode(j,isd)
//                 if (eta2(nd)+dp(nd)<=h0) then


//                   istop=0
//                   do l=1,nne(nd)
//                     ie=indel(l,nd)
//                     if(ie>0) idry_e2(ie)=1
//                   enddo !l
//                 endif
//               enddo !j=1,2 nodes
//             enddo !i=1,nsdf
//             call exchange_e2di(idry_e2)

//           //  Wetting
//             inew=0 !for initializing and counting su2 sv2
//             srwt_xchng(1)=.false. !flag for wetting occurring
//             do i=1,nsdf !aug. domain for updating vel. at interfacial sides (between 2 sub-domains)
//               isd=isdf(i) !must be internal side
//               if(isdel(1,isd)<0.or.isdel(2,isd)<0) cycle !neither element can have interfacial sides
//               if(isdel(1,isd)==0.or.isdel(2,isd)==0) then
//                 write(errmsg,*)'LEVELS1: bnd side:',isdel(:,isd),iplg(isidenode(1:2,isd))
//                 call parallel_abort(errmsg)
//               endif
//               if(idry_e2(isdel(1,isd))+idry_e2(isdel(2,isd))/=1) cycle
//               // 2 end nodes have total depths > h0

//               if(idry_e2(isdel(1,isd))==1) then
//                 ie=isdel(1,isd) !>0
//               else
//                 ie=isdel(2,isd) !>0
//               endif
//               n1=isidenode(1,isd)
//               n2=isidenode(2,isd)
//               nodeA=elnode(1,ie)+elnode(2,ie)+elnode(3,ie)-n1-n2   ! eli: is the 2,ie one right?
//               l0=lindex(nodeA,ie)

//               if(l0==0.or.nodeA==n1.or.nodeA==n2) then
//                 write(errmsg,*)'Frontier node outside, or on the interface:', &
//        &l0,iplg(nodeA),iplg(n1),iplg(n2),itr,it,iths !icolor(nodeA)

//                 write(12,*)'LEVELS1: fatal error message'
//                 do l=1,ns
//                   if(icolor2(l)==1) then
//                     write(12,*)l,iplg(isidenode(1:2,l))
//                     write(12,*)l,ielg(isdel(1:2,l)),idry_e2(isdel(1:2,l)),idry_e(isdel(1:2,l))
//                   endif
//                 enddo //l
//                 do l=1,nea
//                   write(12,*)l,idry_e2(l),idry_e(l)
//                 enddo //l
//                 call parallel_abort(errmsg)
//               endif // end fatal

//               if(eta2(nodeA)+dp(nodeA)>h0) then !all 3 nodes have depths > h0
//                // Check
//                 do j=1,3
//                   nd=elnode(j,ie)
//                   if(eta2(nd)+dp(nd)<=h0) then
//                     write(errmsg,*)'Failed to wet element (13):',ielg(ie),iplg(nd),iplg(nodeA)
//                     call parallel_abort(errmsg)
//                   endif
//                 enddo !j


//                 srwt_xchng(1)=.true.
//                 istop=0
//                 idry_e2(ie)=0

//                 do j=1,2 !sides sharing nodeA
//                   id1=elside(nx(l0,j),ie)
//                   if(icolor2(id1)==0) then

//                   swild2(1,1:nvrt)=su2(1:nvrt,isd)
//                   swild2(2,1:nvrt)=sv2(1:nvrt,isd)

//                   if(inew(id1)==0) then
//                       !vel. only accurate in resident domain
//                       su2(1:nvrt,id1)=swild2(1,1:nvrt) !su2(1:nvrt,isd)
//                       sv2(1:nvrt,id1)=swild2(2,1:nvrt) !sv2(1:nvrt,isd)
//                       inew(id1)=1
//                   else
//                       su2(1:nvrt,id1)=su2(1:nvrt,id1)+swild2(1,1:nvrt)
//                       sv2(1:nvrt,id1)=sv2(1:nvrt,id1)+swild2(2,1:nvrt)
//                       inew(id1)=inew(id1)+1
//                     endif
//                   endif // icolor2(id)==0
//                 enddo // j=1,2
//               endif // eta2(nodeA)+dp(nodeA)>h0
//             enddo // i=1,nsdf; shoreline sides

//             // Compute average vel. for rewetted sides
//             do i=1,ns
//               if(inew(i)/=0) then
//                 su2(1:nvrt,i)=su2(1:nvrt,i)/real(inew(i),rkind)
//                 sv2(1:nvrt,i)=sv2(1:nvrt,i)/real(inew(i),rkind)
//               endif !inew(i)/=0
//             enddo // i=1,ns

//           endif //istop/=2


//           ltmp = false; // for vel. exchange


//           // Enforce wet/dry flag consistency between nodes and elements due to added wet elements

//           idry2=1


//           do i=1,np
//             do j=1,nne(i)
//               if(idry_e2(indel(j,i))==0) then
//                 idry2(i)=0; exit
//               endif
//             enddo !j
//           enddo !i

//           do i=1,nea
//             inew(i)=0 !use for temp. storage of new element wet/dry flags
//             do j=1,3
//               if(idry2(elnode(j,i))==1) inew(i)=1
//             enddo !j
//           enddo !i=1,nea


//           do i=1,ns
//             if(.not.(idry_e2(isdel(1,i))==1.and.(isdel(2,i)==0.or.isdel(2,i)>0.and.idry_e2(max(1,isdel(2,i)))==1))) cycle
//             // Dry side that may need new vel.

//             iwet=0 // flag
//             do j=1,2
//               ie=isdel(j,i)
//               if(ie>0.and.idry_e2(max(1,ie))==1.and.inew(max(1,ie))==0) iwet=1
//             enddo !j

//             if(iwet==1) then !vel. as average
//               sutmp=0._rkind; svtmp=0._rkind; icount=0
//               do m=1,2 !2 elements
//                 ie=isdel(m,i)
//                 if(ie<=0) cycle

//                 do jj=1,3 !3 sides
//                   !Find wet side
//                   isd2=elside(jj,ie)
//                   if(isdel(1,isd2)>0.and.idry_e2(max(1,isdel(1,isd2)))==0.or. &
//      &isdel(2,isd2)>0.and.idry_e2(max(1,isdel(2,isd2)))==0) then !at least one wet element
//                     icount=icount+1


//                     sutmp(1:nvrt)=sutmp(1:nvrt)+su2(1:nvrt,isd2)
//                     svtmp(1:nvrt)=svtmp(1:nvrt)+sv2(1:nvrt,isd2)
//                   endif
//                 enddo !jj
//               enddo !m=1,2; 2 elements

//               ltmp=ltmp.or.icount/=0
//               if(icount/=0) then
// !                srwt_xchng(1)=.true.
//                 su2(1:nvrt,i)=sutmp(1:nvrt)/real(icount,rkind)
//                 sv2(1:nvrt,i)=svtmp(1:nvrt)/real(icount,rkind)
//               endif
//             endif !iwet
//           enddo !i=1,ns


//           idry_e2(1:nea)=inew(1:nea)
//           srwt_xchng(1)=ltmp

//           call mpi_allreduce(srwt_xchng,srwt_xchng_gb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
//           if(srwt_xchng_gb(1)) then
//             allocate(swild(2,nvrt,nsa),stat=istat)
//             if(istat/=0) call parallel_abort('Levels1: fail to allocate (8)')

//             swild(1,:,:)=su2(:,:)
//             swild(2,:,:)=sv2(:,:)


//             su2(:,:)=swild(1,:,:)
//             sv2(:,:)=swild(2,:,:)
//             deallocate(swild)
          

//          // Sync
//          call parallel_barrier

//          if(istop==2) exit loop15

//       end do loop15


//       endif !it/=iths



//       do i=1,ne
//         do j=1,i34(i)
//           nd=elnode(j,i)
//           if(eta2(nd)+dp(nd)<=h0) then
//             idry_e2(i)=1
//             exit
//           endif
//         enddo !j
//       enddo !i


//       // Wet/dry flags for nodes/sides
//       idry2=1; 
//       idry_s2=1
//       do i=1,nea
//         if(idry_e2(i)==0) then
//           idry2(elnode(1:3,i))=0
//         endif
//       enddo !i

//       do i=1,np
//         do j=1,nne(i)
//           if(idry_e2(indel(j,i))==0) then
//             idry2(i)=0; exit
//           endif
//         enddo !j
//       enddo !i



//       do i=1,ns
//         do j=1,2
//           if(isdel(j,i)>0) then; if(idry_e2(isdel(j,i))==0) then
//             idry_s2(i)=0; exit
//           endif; endif
//         enddo !j
//       enddo !i



//       // Reset vel. at dry sides

//       do i=1,nsa
//         if(idry_s2(i)==1) then
//           su2(1:nvrt,i)=0._rkind
//           sv2(1:nvrt,i)=0._rkind
//         endif
//       enddo !i


//       // Limit elevation at dry nodes
//       do i=1,npa
//         if(idry2(i)==1) then
//           !eta2(i)=min(0.d0,-dp(i))
//           eta2(i)=min(eta2(i),-dp(i))
//         endif
//       enddo !i


//       // z-coordinates for nodes

//       do i=1,npa
//         if(ivcor==2) then; if(eta2(i)<=h0-h_s) then
//           write(errmsg,*)'Deep depth dry:',iplg(i)
//           call parallel_abort(errmsg)
//         endif; endif

//         if (idry2(i)==1) then
//           if(ivcor/=1) kbp(i)=0
//         else !wet
//           call zcoor(1,i,kbp(i),z(:,i))
//         endif !wet ot dry
//       enddo !i=1,npa

//       // Compute element bottom index
//       elem->kb = 0;


//       do i=1,nea
//         if(idry_e2(i)/=0) cycle

// !       Wet
//         n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i)
//         if(idry2(n1)/=0.or.idry2(n2)/=0.or.idry2(n3)/=0) then
//           write(errmsg,*)'level1: Element-node inconsistency (0):',ielg(i),idry_e(i), &
//      &iplg(elnode(1:3,i)),idry2(elnode(1:3,i))
//           call parallel_abort(errmsg)
//         endif
//         kbe(i)=min(kbp(n1),kbp(n2),kbp(n3))
//         do k=kbe(i),nvrt
//           ze(k,i)=(z(max(k,kbp(n1)),n1)+z(max(k,kbp(n2)),n2)+z(max(k,kbp(n3)),n3))/3._rkind
//           if(k>=kbe(i)+1) then; if(ze(k,i)-ze(k-1,i)<=0._rkind) then
//             write(errmsg,*)'Weird element (1):',k,i,ze(k,i),ze(k-1,i)
//             call parallel_abort(errmsg)
//           endif; endif
//         enddo !k
//       enddo !i


//       // Compute side bottom index. For wet side and its wet adjacent element,
//       // kbs>=kbe
//       do i=1,nsa
//         kbs(i)=0 !dry
//         if(idry_s2(i)==0) then !wet side with 2 wet nodes
//           n1=isidenode(1,i)
//           n2=isidenode(2,i)
//           if(idry2(n1)/=0.or.idry2(n2)/=0) then
//             write(errmsg,*)'Side-node inconsistency:',it,islg(i),'node:',iplg(n1),iplg(n2), &
//      &eta2(n1),eta2(n2),idry2(n1),idry2(n2),';element:', &
//      &(isdel(j,i),ielg(isdel(j,i)),idry_e2(isdel(j,i)),j=1,2)
//             call parallel_abort(errmsg)
//           endif
//           if(dps(i)+(eta2(n1)+eta2(n2))/2._rkind<=h0) then
//             write(errmsg,*)'Weird side (0):',islg(i),iplg(n1),iplg(n2),eta2(n1),eta2(n2)
//             call parallel_abort(errmsg)
//           endif
//           kbs(i)=min(kbp(n1),kbp(n2))
//           do k=kbs(i),nvrt
//             zs(k,i)=(z(max(k,kbp(n1)),n1)+z(max(k,kbp(n2)),n2))/2
//             if(k>=kbs(i)+1) then; if(zs(k,i)-zs(k-1,i)<=0._rkind) then
//               write(errmsg,*)'Weird side (1):',k,iplg(n1),iplg(n2),z(max(k,kbp(n1)),n1), &
//      &z(max(k,kbp(n2)),n2),z(max(k-1,kbp(n1)),n1),z(max(k-1,kbp(n2)),n2)
//               call parallel_abort(errmsg)
//             endif; endif
//           enddo !k
//         endif !wet side
//       enddo !i=1,nsa


//       // Compute vel., S,T for re-wetted nodes (q2 and xl are fine)
//       if(it/=iths) then
//         do i=1,np 
//           if(idry(i)==1.and.idry2(i)==0) then
//             do k=1,nvrt
//               uu2(k,i)=0._rkind
//               vv2(k,i)=0._rkind
//               ttmp=0._rkind
//               stmp=0._rkind
//               icount=0
//               do j=1,nnp(i)
//                 nd=indnd(j,i) !must be inside the aug. domain
// !               Wet nbrs not affected by this part and so each sub-domain should use same values
//                 if(idry(nd)==0) then !all indices extended
//                   icount=icount+1
//                   uu2(k,i)=uu2(k,i)+uu2(k,nd)
//                   vv2(k,i)=vv2(k,i)+vv2(k,nd)
//                   ttmp=ttmp+tr_nd(1,k,nd) !tnd(k,nd)
//                   stmp=stmp+tr_nd(2,k,nd) !snd(k,nd)
//                 endif
//               enddo !j
//               if(icount==0) then
//                 !Use last wet value
//               else
//                 uu2(k,i)=uu2(k,i)/real(icount,rkind)
//                 vv2(k,i)=vv2(k,i)/real(icount,rkind)
//                 tr_nd(1,k,i)=ttmp/real(icount,rkind)
//                 tr_nd(2,k,i)=stmp/real(icount,rkind)
//               endif
//             enddo !k=1,nvrt
//           endif !rewetted
//         enddo !i=1,np

//       endif !it/=iths


//       // Check wet/dry in ghost zone
//       prwt_xchng(1)=.false.
//       if(it/=iths) then
//         do i=np+1,npa !check ghosts wet/dry
//           if(idry(i)==1.and.idry2(i)==0) then
//             prwt_xchng(1)=.true. !ghost rewetted; need exchange
//             exit
//           endif


//         enddo !i
//       endif !it/=iths

 

//       // update ghost nodes
//       if(prwt_xchng_gb(1)) then
//           allocate(swild(4,nvrt,nsa),stat=istat)
//           if(istat/=0) call parallel_abort('Levels0: fail to allocate swild')

//           swild(1,:,1:npa)=uu2(:,:)
//           swild(2,:,1:npa)=vv2(:,:)
//           swild(3,:,1:npa)=tr_nd(1,:,:) !tnd(:,:)
//           swild(4,:,1:npa)=tr_nd(2,:,:) !snd(:,:)

 
//           uu2(:,:)=swild(1,:,1:npa)
//           vv2(:,:)=swild(2,:,1:npa)
//           tr_nd(1,:,:)=swild(3,:,1:npa)
//           tr_nd(2,:,:)=swild(4,:,1:npa)

//         endif !prwt_xchng_gb
//       endif !nproc>1

//       close(10)

//       // Update wet/dry flags
//       idry=idry2
//       idry_s=idry_s2
//       idry_e=idry_e2



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

      if constexpr(ics==1) {
         sdbt(1,k,i) = side->u2[k] - dt*(node1->bcc[1][k] + node2->bcc[1][k])/2.d0 - dt*side->u2p[k]
         sdbt(2,k,i) = side->v2[k] - dt*(node1->bcc[2][k] + node2->bcc[2][k])/2.d0 - dt*side->v2p[k]
      } else {
         // XXX IMPLEMENT THIS
         // call project_hvec(bcc(1,k,n1),bcc(2,k,n1),pframe(:,:,n1),sframe2(:,:,i),vn1,vn2)
         // call project_hvec(bcc(1,k,n2),bcc(2,k,n2),pframe(:,:,n2),sframe2(:,:,i),tt1,ss1)
         // sdbt(1,k,i)=su2(k,i)-dt*(vn1+tt1)/2.d0-dt*swild2(k,1)
         // sdbt(2,k,i)=sv2(k,i)-dt*(vn2+ss1)/2.d0-dt*swild2(k,2)
      }  // if ics


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

      if (elem->isDry) return;

      // Wet elements.
      // Warning: ghat1=elem average of \breve{G} must include all: Coriolis, atmo. pressure, 
      // tidal potential, horizontal diffusion, and baroclinic etc
      // Remember to update both fb and F 
      // If ics=2, ghat1 are in eframe (ll)
      // cff1-3: coefficients for 2D/3D
      elem->gradEta = {0.0, 0.0};  // \nabla \eta
      elem->gradPr  = {0.0, 0.0};  // \nabla pressure
      elem->gradETP = {0.0, 0.0};  // \nabla etp

      double av_cff1 = 0.0;  // av. for cff1
      double av_cff2 = 0.0;  // av. for cff2
      double av_cff3 = 0.0;  // av. for cff3
      double av_df = 0.0;  // av. for H^\alpha / H
      double av_cff2_chi = 0.0;  // av. for cff2*chi
      double av_cff3_chi = 0.0;  // av. for cff3*chi

      const int i34 = elem->i34;
      const int kb = elem->kb;

      for (int j=0; j<i34; j++) {
         // For each node and side in the element...

         const Node *node = elem->nodes[j];
         const Side *side = elem->sides[j];

         const double chi      = side->chi;
         const double veg_c2   = side->veg_c2;
         const double veg_beta = side->veg_beta;
         const double veg_h    = side->veg_h();
         const double veg_c    = veg_c2/(1.0 + veg_c2);  // veg_c2 = \alpha*|u|*dt

         // Node - const in elem
         elem->gradEta += node->eta2*elem->gradPhi[j];  // in eframe if ics=2
         elem->gradPr  += node->pr  *elem->gradPhi[j];

         if (elem->depth >= tip_dp) {
           elem->gradETP += node->ETP*elem->gradPhi[j];
         }

         
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
                  cff1=1;
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
            for (int k = side->kb; k < nvrt; k++) {
            
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
         

         // Add btrack & surface stress terms to \breve{G} first
         if constexpr(ics==1) {

            elem->ghat1 += cff1*Ustar - cff3*Ustar_alpha;

         } else {
            // XXX IMPLEMENT THIS
            // call project_hvec(xtmp,ytmp,sframe2(:,:,isd),eframe(:,:,i),tt1,ss1)
            // elem->ghat1[0] += cff1*vn1 - cff3*tt1;
            // elem->ghat1[1] += cff1*vn2 - cff3*ss1;
         }  // if ics

         // Surface and bottom stress terms
         if constexpr(ics==1) {
            const Vector2 tau    = side->tau();
            const Vector2 ubstar = side->btu[side->kb + 1]; // XXX CHECK INDICES

            // Term for surface and bottom stress
            elem->ghat1 -= chi*dt*cff2*ubstar + dt*cff1*tau;

         } else {
            assert(false && "computeElemAverageGhat: Not implemented");
            // XXX IMPLEMENT THIS
            // n1=isidenode(1,isd); n2=isidenode(2,isd)
            // call project_hvec(tau(1,n1),tau(2,n1),pframe(:,:,n1),eframe(:,:,i),vn1,vn2)
            // call project_hvec(tau(1,n2),tau(2,n2),pframe(:,:,n2),eframe(:,:,i),tt1,ss1)
            // tau_x=(vn1+tt1)/2.d0
            // tau_y=(vn2+ss1)/2.d0
            // call project_hvec(sdbt(1,kbs(isd)+1,isd),sdbt(2,kbs(isd)+1,isd),sframe2(:,:,isd),eframe(:,:,i),ubstar,vbstar)
         }  // if ics

         

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

         fb += side->nuGradU[side->kb];

         // Add to ghat1
         if constexpr(ics==2) {
            // XXX IMPLEMENT THIS
            // call project_hvec(bigf1,bigf2,sframe2(:,:,isd),eframe(:,:,i),vn1,vn2)
            // bigf1=vn1; bigf2=vn2
            // call project_hvec(botf1,botf2,sframe2(:,:,isd),eframe(:,:,i),vn1,vn2)
            // botf1=vn1; botf2=vn2
            // call project_hvec(bigfa1,bigfa2,sframe2(:,:,isd),eframe(:,:,i),vn1,vn2)
            // bigfa1=vn1; bigfa2=vn2
         } // if ics

         elem->ghat1 += cff1*dt*F - cff2*chi*dt*dt*fb - cff3*dt*Fa;
      } // for j (nodes and sides)

      elem->ghat1 = (1.0/i34)*elem->ghat1;
         
      // Finish off terms in F, F^\alpha and f_b
      const double grav3 = elem->mean(&Node::grav2);
      const Vector2 fb = grav3*elem->gradETP - elem->gradPr/rho0; //const in each elem


      Vector2 tmp = {0.0, 0.0}; //elem average of all terms; into ghat1
      for (int j=1; j<i34; j++) { 
         // For each side in the element...
         const Side *side = elem->sides[j];
         const double chi = side->chi;

         const double htot = side->htot2();
         const double veg_h = side->veg_h();
         const Vector2 F  = htot*fb;
         const Vector2 Fa = min(htot, veg_h)*fb;

         tmp += dt*(av_cff1*F - av_cff2*chi*dt*fb - av_cff3*Fa);
      } // for j

      //botf1 etc are already in elem frame
      elem->ghat1 += tmp/i34;

      // Baroclinic force
      if constexpr(ibc==0) {
         assert(elem->nodes[0]->rho[0] >= -98.0 && elem->nodes[1]->rho[0] >= -98.0 && elem->nodes[2]->rho[0] >= -98.0 && "computeElemAverageGhat1: Impossible dry");

         Vector2 Fa = {0.0, 0.0};
         Vector2 intGradRho = {0.0, 0.0};
         
         // \sum_{l=k}^N dr*dz; whole level (and eframe if ics=2)
         for (int k=nvrt-2; k > kb + 1; k--) { /// XXXX CHECK INDICES
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
                           double h_c, double h_s, double theta_b, double theta_f, double *sigma_lcl, double *sigma, 
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
                   h_c, h_s, theta_b, theta_f, sigma_lcl, sigma, ztot, cs, dcs);

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
   double *drag, double *roughness, double *rmanning) {

   mesh->populateNodes(i34, isbnd, sharedIdxNodeBoundary, 
                       node_coords, node_coords_ll, use_ll, wind_rotate_angle, 
                       nodeNeighNodesStart, nodeNeighNodesIndices,
                       nodeNeighElemsStart, nodeNeighElemsIndices,
                       localNodeIndices, localElemIndices, localSideIndices,
                       localNodePart, localElemPart, localSidePart,
                       drag, roughness, rmanning);
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
                                                double *vmo, double *vfa) {
   mesh->populateBoundaryNodes(nond, iond, emo, efa, umo, ufa, vmo, vfa);

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
      case (int)VarCode::lonlat:
         *offset = (int)((char*)&elems[0].xlon - (char*)elems);
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
      case (int)VarCode::w:
         *offset = (int)((char*)&sides[0].w2 - (char*)sides);
         *sizeType = sizeof(double);
         break;
      case (int)VarCode::u:
         *offset = (int)((char*)&sides[0].u2[0].x - (char*)sides);
         *sizeType = sizeof(Vector2);
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
         *offset = (int)((char*)&nodes[0].z[0] - (char*)nodes);
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
      default:
         *offset = -1;    
   }  
   
}

extern "C" __global__ void getElementVar(int varCode, bool is3D, double *var1, double *var2, int *ivar1, int *ivar2, int k = -1, int j = -1) {
   // Get the value of a variable from an element
   // varName is the name of the variable to get
   // var is the array to store the value
   // k is the vertical level index of the variable to get

   ELEM(elem);
   
   int idx = elem->idx;
   int i = threadIdx.x + blockIdx.x * blockDim.x;

   // printf("getElementVar: varCode = %d, k = %d, idx = %d, i = %d\n", varCode, k, idx, i);
   if      (varCode == (int)VarCode::w)       { var1[idx] = elem->w[k];}
   else if (varCode == (int)VarCode::z)       { var1[idx] = elem->z[k];}
   else if (varCode == (int)VarCode::x)       { var1[idx] = elem->center.x; var2[idx] = elem->center.x;}
   else if (varCode == (int)VarCode::area)    { var1[idx] = elem->area;}
   else if (varCode == (int)VarCode::depth)   { var1[idx] = elem->depth;}
   else if (varCode == (int)VarCode::radius)  { var1[idx] = elem->radius;}
   else if (varCode == (int)VarCode::lonlat)  { var1[idx] = elem->xlon;     var2[idx] = elem->ylat;}
   else if (varCode == (int)VarCode::ghat)    { var1[idx] = elem->ghat1.x;  var2[idx] = elem->ghat1.y;}
   else if (varCode == (int)VarCode::kb)      { ivar1[idx] = elem->kb;}
   else if (varCode == (int)VarCode::ssign && j > -1)  {ivar1[idx] = elem->ssign[j];}
   else if (i==0) printf("getElementVar: varCode = %d, k = %d not found\n", varCode, k);
      
   
}


extern "C" __global__ void getNodeVar(int varCode, double *var, int k = -1) {
   // Get the value of a variable from an node
   // varName is the name of the variable to get
   // var is the array to store the value
   // k is the vertical level index of the variable to get

   NODE(node);
   
   int idx = node->idx;
   int i = threadIdx.x + blockIdx.x * blockDim.x;

   // printf("getNodeVar: varCode = %d, k = %d, idx = %d, i = %d\n", varCode, k, idx, i);
   
   if      (varCode == (int)VarCode::z && k >= 0)       var[idx] = node->z[k];
   else if (varCode == (int)VarCode::x && k == -1)      var[idx] = node->p.x;
   else if (varCode == (int)VarCode::depth && k == -1)  var[idx] = node->depth;
   else if (varCode == (int)VarCode::htot && k == -1)   var[idx] = node->htot2();
   else if (varCode == (int)VarCode::kb && k == -1)     var[idx] = (double)(node->kb);
   else if (varCode == (int)VarCode::elbc && k == -1)   var[idx] = node->nodeBnd!=nullptr?node->nodeBnd->elbc:0.0;
   else if (varCode == (int)VarCode::eta && k == -1)    var[idx] = node->eta2;  
   else if (varCode == (int)VarCode::isdry && k == -1)  var[idx] = node->isDry?1.0:0.0;
   else if (varCode == (int)VarCode::qel && k == -1)    var[idx] = node->qel;
   else if (varCode == (int)VarCode::u && k >= 0)       var[idx] = node->u2[k].x;
   else if (varCode == (int)VarCode::w && k >= 0)       var[idx] = node->w2[k];
   
   else if (i==0) printf("getNodeVar: varCode = %d, k = %d not found\n", varCode, k);
      
   
}


extern "C" __global__ void getSideVar(int varCode, double *var, int k = -1) {
   // Get the value of a variable from an side
   // varName is the name of the variable to get
   // var is the array to store the value
   // k is the vertical level index of the variable to get
   SIDE(side);
   
   int idx = side->idx;
   int i = threadIdx.x + blockIdx.x * blockDim.x;

   if      (varCode == (int)VarCode::z && k >= 0)        var[idx] = side->zc[k];
   else if (varCode == (int)VarCode::x && k == -1)       var[idx] = side->center.x;
   // else if (varCode == (int)VarCode::y && k == -1)       var[idx] = side->center.y;
   else if (varCode == (int)VarCode::isdry && k == -1)   var[idx] = side->isDry?1.0:0.0;
   else if (varCode == (int)VarCode::nx && k == -1)      var[idx] = side->n.x;
   else if (varCode == (int)VarCode::ny && k == -1)      var[idx] = side->n.y;
   else if (varCode == (int)VarCode::u && k >= 0)        var[idx] = side->u2[k].x;
   // else if (varCode == (int)VarCode::v && k >= 0)        var[idx] = side->u2[k].y;
   else if (varCode == (int)VarCode::U && k == -1)       var[idx] = side->U.x;
   // else if (varCode == (int)VarCode::V && k == -1)       var[idx] = side->U.y;   
   else if (varCode == (int)VarCode::depth && k == -1)   var[idx] = side->depth();
   else if (varCode == (int)VarCode::kb && k == -1)      var[idx] = (double)(side->kb);
   else if (varCode == (int)VarCode::n1 && k == -1)      var[idx] = side->n1->idx;
   else if (varCode == (int)VarCode::n2 && k == -1)      var[idx] = side->n2->idx;
   else if (varCode == (int)VarCode::eta && k == -1)     var[idx] = side->eta2();
   else if (varCode == (int)VarCode::btu && k >= 0)      var[idx] = side->btu[k].x;
   // else if (varCode == (int)VarCode::btv && k >= 0)      var[idx] = side->btu[k].y;
   else if (varCode == (int)VarCode::bcc && k >= 0)     var[idx] = side->bcc[k].x;
   // else if (varCode == (int)VarCode::bccy && k >= 0)     var[idx] = side->bcc[k].y;
   else if (varCode == (int)VarCode::ubc && k >= 0)      var[idx] = side->sideBnd!=nullptr?side->sideBnd->ubc[k].x:0.0;
   // else if (varCode == (int)VarCode::vbc && k >= 0)      var[idx] = side->sideBnd!=nullptr?side->sideBnd->ubc[k].y:0.0;
   else if (varCode == (int)VarCode::uth && k >= 0)      var[idx] = side->sideBnd!=nullptr?side->sideBnd->uth(k).x:0.0;
   // else if (varCode == (int)VarCode::vth && k >= 0)      var[idx] = side->sideBnd!=nullptr?side->sideBnd->uth(k).y:0.0;
   else if (varCode == (int)VarCode::ifltype && k == -1) var[idx] = side->sideBnd!=nullptr?side->sideBnd->segment->ifltype:0.0;
   else if (varCode == (int)VarCode::island && k == -1)  var[idx] = side->sideBnd!=nullptr?(double)(side->sideBnd->segment->isLand):7.7;
   // else if (varCode == (int)VarCode::island && k == -1)  var[idx] = side->sideBnd != nullptr ? (double)((unsigned long long)(side->sideBnd->segment) & 0xFFFFFFFF) : 7.7;
   else if (varCode == (int)VarCode::hhat && k == -1)    var[idx] = side->hhat;
   // else if (varCode == (int)VarCode::gradEta1x   && k == -1)    var[idx] = side->gradEta1.x;
   // else if (varCode == (int)VarCode::gradEta1y   && k == -1)    var[idx] = side->gradEta1.y;
   // else if (varCode == (int)VarCode::gradEta2x   && k == -1)    var[idx] = side->gradEta2.x;
   // else if (varCode == (int)VarCode::gradEta2y   && k == -1)    var[idx] = side->gradEta2.y;
   // else if (varCode == (int)VarCode::gradETPx   && k == -1)    var[idx] = side->gradETP.x;
   // else if (varCode == (int)VarCode::gradETPy   && k == -1)    var[idx] = side->gradETP.y;
   // else if (varCode == (int)VarCode::gradPrx   && k == -1)    var[idx] = side->gradPr.x;
   // else if (varCode == (int)VarCode::gradPry   && k == -1)    var[idx] = side->gradPr.y;
   else if (i==0) printf("getSideVar: varCode = %d, k = %d not found\n", varCode, k);
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