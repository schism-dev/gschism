


// Each device has a global mesh object
__device__ Mesh* mesh = nullptr;



Polynomial<4> *polyDischarge; // XXX read from file

//===============================================================================
//     
//===============================================================================
__device__ __forceinline__ void update_bdef(double time, double x0, double y0, double &dep, Vector3 &vel) {
   // For imm=2, user needs to update bottom depth and velocity

   dep = min(1.0, 7.0 - (x0 + time));
   vel = { -1.0, 0.0, 0.0 };

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
      const Point2 x3 = side2->n1->point2();
      const Point2 x4 = side2->n2->point2();

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

            nextElement = nextSide->otherElement(this);
            tIntersect = max(small, min(tt1, 1.0-small));
            t = tt2;

            found_intersection = true;
         }
      
      }

   }


   // found_intersection should always be true, since we checked all the sides.  
   side = nextSide;
// if (found_intersection) {
// if (nextElement == nullptr)  printf("$$$  0, %d, %.5f, %.5f\n",  side->idx, side->x(t), side->y(t));
// else   printf("$$$  %d, %d, %.5f, %.5f\n", nextElement->idx, side->idx, side->x(t), side->y(t));
// }   
   return found_intersection;
}




__device__ bool Element::findNextIntersection(Point2 x1, Point2 x2, Side *&side, Element *&nextElement,
                                              double &t, double& tIntersect) {

   double closest_intersect_t = 1e10;
   bool found_intersection = false;

   for (int i = 0; i < i34; i++) {

      Side *side2 = sides[i];

      // (x3,y3)-(x4,y4) is the element side line
      const Point2 x3 = side2->n1->point2();
      const Point2 x4 = side2->n2->point2();

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
            side = side2;
            nextElement = side->otherElement(this);
            t = tt2;
            tIntersect = max(small, min(tt1, 1.0-small));
            found_intersection = true;
         }
      
      }

   }

// if (found_intersection) {
// if (nextElement == nullptr)  printf("$$$  0, %d, %.5f, %.5f\n",  side->idx, side->x(t), side->y(t));
// else   printf("$$$  %d, %d, %.5f, %.5f\n", nextElement->idx, side->idx, side->x(t), side->y(t));
// }

   // found_intersection should always be true, since we checked all the sides.  
   return found_intersection;
}


__device__ void Element::trackGoesOutOfWater(Side *nextSide, Element *curElem, double tt, int jlev,
                                             double &dt, Point3 &xt, double zint, Vector3 &u_int) {

   printf("Tracking segment goes out of water at element %d\n", curElem->idx);

   // When this happens, we stop our algorithm that usually assumes that we move to the next element.
   // Instead, we nudge the intersection point to the barycenter of the element, and nextElement becomes
   // curElement.

   // Then it finds the tangential velocity at the element side (the side where the intersection was found),
   // And uses this velocity instead to continue the tracking.

   // computes the intersection point in cartesian coordinates
   double xint = nextSide->x(tt);
   double yint = nextSide->y(tt);
  
   // Nudges the point to the barycenter of the element. 
   // This ensures that the point is inside the element and hence, we will find an intersection if the
   // next segment intersects the element sides.  
   curElem->nudgeToBarycenter(xint, yint);

   // Computes the tangential velocity at the side. Note that we are not using the x,y,z location
   // of the intersection point. We are approximating it by the velocity at the side center at jlev.
   // NOTE: With the new code it will be easy to compute the tangential velocity at the intersection point.
   double vtan = dot(nextSide->t(), nextSide->u2[jlev]);

   // If the side is an open boundary, the tangential velocity is zero.
   if (nextSide->isBnd()) vtan = 0.0;

   // The new velocity components are computed from the tangential velocity and the tangent vector.
   const double us = vtan*nextSide->t().x;
   const double vs = vtan*nextSide->t().y;
   const double ws = 0.5*(nextSide->n1->u2[jlev].z + nextSide->n2->u2[jlev].z);

   // Compute the new "tracked" time based on the remaining segment length and the original velocity.
   const double len2 = (xint - xt.x)*(xint - xt.x) + (yint - xt.y)*(yint - xt.y);
   const double vel2 = u_int.x*u_int.x + u_int.y*u_int.y;        
   assert(vel2 > 0.0 && "Zero velocity found");
   dt = dt - sqrt(len2/vel2); 

   // The end of this tracking segment is now the intersection point.
   xt = {xint, yint, zint};
   
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
                                                       Element *&nextElem, int &jlev, double &dt) {
// printf("trackThroughElementsAndLevels: x0=%f, y0=%f, xt=%f, yt=%f\n", x0, y0, xt, yt);
   // First thing is to check if the target point is inside the element. If it is there is nothing
   // else to do.
   if (containsPoint(toPoint2(xt))) {
      nextElem->interpolateVelocities(xt, u_int);
      return true;
   }

   Side *side = nullptr;

   // 1. Initial checks
   // Verify starting element is valid (wet)
   assert(!nextElem->idry && "Starting element is dry");

   // Checks if track has a zero length.
   const double len2 = (xt.x - x0.x)*(xt.x - x0.x) + (xt.y - x0.y)*(xt.y - x0.y);   
   assert(len2 > 0.0 && "Zero length track"); 
   

   const int MAX_ITERATIONS = 1000;
   
   double t;
   Side *nextSide = nullptr;
   double tIntersect = 1e10;
   double zint = x0.z;

   // Horizontal tracking loop
   // Finds the first intersection with the element sides. It is guaranteed to exist because
   // (x0,y0) is inside the element but the target point is not.
   bool found_intersection = findNextIntersection(toPoint2(x0), toPoint2(xt), nextSide, nextElem, t, tIntersect);
   assert(found_intersection && "No intersection found");

   zint = (1 - tIntersect)*zint + tIntersect*xt.z;

   // if the tracking segment goes out of water (either dry or edge boundary)
   if (isEdgeBoundary(nextSide) || (nextElem->idry)) {
      trackGoesOutOfWater(nextSide, this, t, jlev, dt, xt, zint, u_int);
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
      zint = (1 - tIntersect)*zint + tIntersect*xt.z;

      // If we are already in the target element, we are done.
      // Remember that, if we are here, the point (xt,yt) is inside curElem.
      if (!found_intersection) break;


  
      // if the tracking segment goes out of water (either dry or edge boundary).
      if (curElem->isEdgeBoundary(nextSide) || (nextElem->idry)) {

         trackGoesOutOfWater(nextSide, curElem, t, jlev, dt, xt, zint, u_int);

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
   nextElem->interpolateVelocities(xt, u_int);


   return true;
}




__inline__ __device__ double Element::computeSignedArea(void) const {
   // Compute the area of the element

   double res = signedArea(nodes[0]->point2(), nodes[1]->point2(), nodes[2]->point2());

   // In the quad case we need to add the area of the other triangle.
   if (isQuad()) {
      double ar1 = signedArea(nodes[0]->point2(), nodes[2]->point2(), nodes[3]->point2());
      double ar2 = signedArea(nodes[0]->point2(), nodes[1]->point2(), nodes[3]->point2());
      double ar3 = signedArea(nodes[1]->point2(), nodes[2]->point2(), nodes[3]->point2());
      if (ar1 < 0.0 || ar2 < 0.0 || ar3 < 0.0) {
         // XXX Aren't we computing a signed area? why is it not allowed to be negative?
         Error("Negative area found in quad element");
      }
      res += ar1;
   }

   return res;
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
__device__ bool Element::computeShapeCoords(Point2 xt, double shpCoords[4],
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
         x0 = nodes[0]->point2();
         x1 = nodes[1]->point2();
         x2 = nodes[2]->point2();
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
         const double a = - tmpmin + 1.0e-2;  // XXX it used to be 1.0e-4 in the comments of the original code
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

}


/**  Interpolate velocities at a given point (xt,yt,zt)

   @param xt, yt, zt: coordinates of the point to interpolate
   @param uuint: velocity interpolated at the final position (xt,yt,zt)

   This is done in four steps:

   1. Find the shape coordinates for the final position in the element.
   2. Interpolate the z level of the top and bottom levels at xt, yt. 
   3. Perform a *binary search* to find the vertical level containing zt. It will have to interpolate the 
      z levels at xt, yt at each step.   
   4. Interpolate the velocities at the final position (xt,yt,zt). This is done by interpolating the velocities
      horizontally and vertically.
**/
__device__ void Element::interpolateVelocities(Point3 xt, Vector3& u_int) {
   // Interpolate the velocities at the point xt. 
   // This function works the same for triangular and quadrilateral elements.
   double shpCoord[4];


   // First step is to compute the shape coordinates of (xt,yt)
   computeShapeCoords(toVector2(xt), shpCoord, false);
   assert((min4(shpCoord[0], shpCoord[1], shpCoord[2], shpCoord[3]) >= -small) && "Invalid shape coordinates");

   // Second step is to find the level that contains zt
   // Uses a binary search to find the level
   // Note that z coordinates are interpolated from their values at the nodes of the element.
   int L = kbe;   // I use L(eft) and R(ight) as usual, but here L is the bottom level and R is the top level
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
// printf("$$$ wght1 = %f, wght2 = %f, zt = %f, zL = %f, zR = %f\n", wght1, wght2, zt, zL, zR);
// printf("$$$ shpCoord = %f, %f, %f, %f\n", shpCoord[0], shpCoord[1], shpCoord[2], shpCoord[3]);
   // Interpolate horizontally and vertically the velocities at the two levels L, R
   u_int = wght1*interpolate(L, shpCoord, &Node::u2) + wght2*interpolate(R, shpCoord, &Node::u2);

   // Neither of the interpolated velocities should be NaN
   assert(!isnan(u_int.x) && !isnan(u_int.y) && !isnan(u_int.z) && "NaN velocities in interpolateVelocities");

}



__inline__ __device__ double Side::x(double t) const { assert(t >= 0.0 && t <= 1.0); return (1-t)*n1->x + t*n2->x; }
__inline__ __device__ double Side::y(double t) const { assert(t >= 0.0 && t <= 1.0); return (1-t)*n1->y + t*n2->y; }

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
 * @param sclr: Array for backtracked values
 * @return: True if the backtracking was successful, false otherwise
 */
__inline__ __device__ bool Side::btrack(
      int jlev,
      double dt,            // Tracking time
      double dtbk,          // Tracking timestep
      double* sclr,         // Array for backtracked values
      Vector3 &u_int
) const {

   // Start from side center
   Point3 x0 = {xc, yc, zc[jlev]};
 

   double time_rm = dt; // Remaining time


   // // Initialize outputs
   // for(int i = 0; i < 4 + N_TRACERS; i++) {
   //    sclr[i] = 0.0;
   // }


   // Initial velocities
   u_int = {u2[jlev].x, u2[jlev].y, 0.5*(n1->u2[jlev].z + n2->u2[jlev].z)};

   // Find the element that is in the direction opposite to the velocity
   Element* curElem = findElementTowards({u_int.x, u_int.y});

printf("$$$  %d, %d, %.5f, %.5f\n",  curElem->idx,idx, x0.x, x0.y)   ;
   const int MAX_ITERATIONS = 100000; // XXX 1000
   Element* nextElem = curElem;
   int iter = 0;
   while(true) {
      // Calculate timestep
      double dtb = min(dtbk, time_rm);

      // Next point in the backtracking curve
      Point3 xt = x0 - dtb*u_int;

      bool found = curElem->trackThroughElementsAndLevels(x0, xt, u_int, nextElem, jlev, dtb);
      assert(found && "Backtracking failed in btrack");
// printf("<>T = %f, %f %f %f\n", time_rm - dtb, xt, yt, zt);
      time_rm -= dtb; // if trackThroughElementsAndLevels terminates early (if it intersects a boundary or a dry element),
                      // then time_rm needs to be increased by the unused time.


      if (time_rm == 0.0) break; // from the definition of dtb, time_rm will reach exactly 0.0

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








// extern "C" {
// __global__ void backtrackKernel(
//     const Mesh* mesh,
//     const double dt,
//     BacktrackResults* results  // Structure to store results
// ) {
//     // Get side index
//     int sIdx = blockIdx.x * blockDim.x + threadIdx.x;
//     if (sIdx >= mesh->ns) return;  // Check bounds for local sides
    
//     // Get local side data
//     const Side* side = mesh->localSideData + sIdx;
//     if (side->isDry()) return;  // Skip dry sides
    
//     // Get starting position (side center)
//     double x0 = side->xc;
//     double y0 = side->yc;
//     double z0 = side->zc;
    
//     // Get starting element
//     Element* startElement = side->getStartElement();
    
//     // Do backtracking
//     double xt, yt, zt;  // Final positions
//     bool exitedDomain = startElement->backtrack(
//         x0, y0, z0,
//         dt,
//         xt, yt, zt
//     );
    
//     // Store results
//     results[sid].finalPos = {xt, yt, zt};
//     results[sid].exitedDomain = exitedDomain;
// }

__host__ __device__ double Node::sigmaTransitionDepth(void) const { return min(depth, mesh->h_s); }  

__device__ double Node::compZCoord(int k) const {
   // Computes the z coordinate at level k.
   // Depending on the vertical coordinate system, the z coordinate is computed differently.
   // When ivcor=1, the z coordinate is computed using the localized sigma levels.
   // When ivcor=2, the z coordinate is computed using the zeta/sigma levels. Z levels are
   // used 
   
   // NOTE: Assumes that updateZCoordData() has been called before this function.
      
   if (idry) {
      return 0.0;
   }

   if constexpr(ivcor == 2) {

      const double hmod2 = sigmaTransitionDepth(); // Depth transition between z and sigma levels
      const int kz = mesh->kz;

      // TODO: Here we are assuming that the bottom level, kbpl, was found in a previous step.
      // The FORTAN code computes it every time this function is called, but it should be done only 
      // once per time step/node.
      // CHECK: updateZCoordData
      if (k == kz)               return -hmod2;
      else if (k == nvrt-1)      return eta2;
      else if (k == kb)          return -depth;
      else if (k < kz && k >= 0) return mesh->ztot[k];
      else if (k > kz) {
         const int ksig = k - kz + 1;

         if (hmod2 <= mesh->h_c) {
            return mesh->sigma[ksig] * (hmod2 + eta2) + eta2;
         } 
         else if (eta2 <= - mesh->h_c - (depth - mesh->h_c) * mesh->theta_f/mesh->sinh_theta_f) {
            printf("Please choose larger h_c: %f %f\n", eta2, mesh->h_c);
            Error("ZCOOR: invalid h_c\n");
         } 
         else {
            return eta2*(1.0 + mesh->sigma[ksig]) + 
                         mesh->h_c * mesh->sigma[ksig] + 
                         (hmod2 - mesh->h_c) * mesh->cs[ksig];
         }
      } else { 
         Error("ZCOOR: invalid level or coordinate system\n");
      }

   }
   else if constexpr(ivcor == 1) {
      if (k == kb) return -depth;
      else if (k == nvrt-1) return eta2;
      else if (k >= kb) {
            return (eta2 + depth) * mesh->sigma_lcl[k] + eta2;
      }
   }
   
   printf("ZCOOR: invalid level or coordinate system\n");
   return 0.0;
}

#ifdef TESTING
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
      const int kbs = side->kb;
      
      // Process each vertical level
      for (int k = kbs + 1; k < nvrt; k++) {
          if (isLand) {
              side->flux_adv_hface[k] = 0.0;

          } else {
              // Calculate normal component of the velocity at current and previous levels
              const double vnor1 = dot(side->u2[k  ], side->n);
              const double vnor2 = dot(side->u2[k-1], side->n);
              
              // Calculate flux (normal * area in local x-direction)
              const double height = side->zc[k] - side->zc[k-1]; // Vertical distance between the two side centers
              const double vmean = 0.5*(vnor1 + vnor2);
              const double length = side->length;
              
              side->flux_adv_hface[k] = height*length*vmean;
          }
      }
   }
}

#endif 

__device__ void Mesh::populateMesh(int *i34, int *elnode, int *elside,
                double *node_coords, double *wind_rotate_angle, 
                int *isdel, int *isidenode, int *isbnd,
                int *nond, int *nlnd, int *iond, int *ilnd,
                int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                int *localNodePart,     int *localElemPart,     int *localSidePart) {

   printf("Populating mesh for partition %d, %d, %d\n", partNum, nond, nlnd);

   // Horizontal grid
   printf("Processing %d nodes\n", np);      
   for (int i = 0; i < np; i++)
   {
      Node &node = localNodeData[i];

      const int idxStart = nodeNeighNodesStart[i];
      const int idxEnd   = nodeNeighNodesStart[i+1];
      node.numNodes = idxEnd - idxStart;
      node.nodes = (Node **)(nodeNeighNodes + idxStart*sizeof(Node *));
      node.idx = i;

      for (int idx = idxStart; idx < idxEnd; idx++) {
         int idxNode = nodeNeighNodesIndices[idx];
         Node *node_ptr = (Node *)&(nodes[localNodePart[idxNode]][localNodeIndices[idxNode]]);

         nodeNeighNodes[idxNode] = node_ptr;
         node.nodes = &(nodeNeighNodes[idxNode]);
      }

      // Element **elements; // pointer to the elements connected to this node
      // int numElements;   // number of elements shared by this node

      node.x = node_coords[i * 3];
      node.y = node_coords[i * 3 + 1];
      node.z = node_coords[i * 3 + 2];  // XXX this is probably wrong
      node.depth = node_coords[i * 3 + 2];

      // Links the node to the boundary information
      int idxNodeBoundary = 0;
      if (isbnd[i] >= 0) {
         node.nodeBnd = &(localNodeBoundaryMem[idxNodeBoundary]);
         if (isbnd[i]>0) {

            node.nodeBnd->segment = &(localSegmentMem[isbnd[i]]);
            node.nodeBnd->segment->ifltype = isbnd[i];

            // Initialize the length and area of the segment to nan to indicate that they are not yet computed
            node.nodeBnd->segment->clen  = nan("1");
            node.nodeBnd->segment->carea = nan("1");

         }
         else {
            node.nodeBnd->segment = (BoundarySegment *)(isbnd[i]); // Hack, isbnd[i] will be -1 or -2 in this case
         }
         
         idxNodeBoundary++;
      }
      else {
         node.nodeBnd = nullptr;
      }

// XXXX
constexpr double rot_per=3000.0;
constexpr double rot_f=2.0*pi/rot_per;
for (int k = 0; k < nvrt; k++) {node.u2[k] = {0.0, -node.y*rot_f, node.x*rot_f};}
// XXXXX
   }

   printf("Processing %d sides\n", ns);
   int idxSideBoundary = 0;
   for (int i = 0; i < ns; i++)
   {
      Side &side = localSideData[i];
      side.idx = i;
      // finding the pointer to the nodes that form the side
      const int idx_n1 = isidenode[globalSideIndices[i]*2];
      const int idx_n2 = isidenode[globalSideIndices[i]*2 + 1];

      const int locidx_n1 = localNodeIndices[idx_n1];
      const int locidx_n2 = localNodeIndices[idx_n2];
      const int part_n1 =   localNodePart[idx_n1];
      const int part_n2 =   localNodePart[idx_n2];

      side.n1 = (Node *)&(mesh->nodes[part_n1][locidx_n1]);
      side.n2 = (Node *)&(mesh->nodes[part_n2][locidx_n2]);

      // Normal vector to the side
      side.n = {- (side.n1->y - side.n2->y), (side.n1->x - side.n2->x)};
      const double norm = sqrt(dot(side.n, side.n));
      side.n = side.n/norm;

      // Side boundary type
      const NodeBoundary *n1bnd = side.n1->nodeBnd;
      const NodeBoundary *n2bnd = side.n2->nodeBnd;
      // XXXX NOTE: I'm chosing the boundary type of the first node. This may not be correct.
      
      if (n1bnd != nullptr || n2bnd != nullptr) 
      {
         side.sideBnd = &(localSideBoundaryMem[idxSideBoundary]);
         if (n1bnd != nullptr) side.sideBnd->segment = n1bnd->segment;
         else                  side.sideBnd->segment = n2bnd->segment;
         side.sideBnd->side = &side;
         idxSideBoundary++;
      }
      else side.sideBnd = nullptr;

      // finding the pointer to the elements that share the side
      const int idx_e1 = isdel[globalSideIndices[i]*2];
      const int idx_e2 = isdel[globalSideIndices[i]*2 + 1];

      if (idx_e1 > -1) {
         const int locidx_e1 = localElemIndices[idx_e1];
         const int part_e1 =  localElemPart[idx_e1];
         side.elems[0] = (Element *)&(elems[part_e1][locidx_e1]);
      }
      else {
         side.elems[0] = nullptr;
      }

      if (idx_e2 > -1) {
         const int locidx_e2 = localElemIndices[idx_e2];
         const int part_e2 =  localElemPart[idx_e2];
         side.elems[1] = (Element *)&(elems[part_e2][locidx_e2]);
      }
      else {
         side.elems[1] = nullptr;
      }



      // Compute the center of the side
      side.xc = 0.5*(side.n1->x + side.n2->x);
      side.yc = 0.5*(side.n1->y + side.n2->y);

// XXXX
constexpr double rot_per=3000.0;
constexpr double rot_f=2.0*pi/rot_per;
for (int k = 0; k < nvrt; k++) side.u2[k] = {-side.yc*rot_f , side.xc*rot_f};
// XXXXX      

      // Compute the z level at the center of the side
      for (int k = 0; k < nvrt; k++) side.zc[k] = 0.5*(side.n1->znl[k] + side.n2->znl[k]);

      // Compute the length of the side
      const Vector2 vLen = side.n2->point2() - side.n1->point2();
      side.length = sqrt(dot(vLen, vLen));

   }

   printf("Processing %d elements\n", ne);
   for (int i = 0; i < ne; i++)
   {
      Element &elem = localElemData[i];
      elem.idx = i;
      for (int j = 0; j < i34[i]; j++) {
         const int idxNode = elnode[i*4 + j];  
         const int idxSide = elside[i*4 + j];
         
         elem.nodes[j] = (Node *)&(nodes[localElemPart[idxNode]][localElemIndices[idxNode]]);
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
      elem.wind_angle = wind_rotate_angle[i];

      // Compute the centroid of the element and its depth
      elem.xc = 0.0;
      elem.yc = 0.0;
      elem.depth = 0.0;
      for (int j = 0; j < i34[i]; j++) {
         elem.xc    += elem.nodes[j]->x/i34[i];
         elem.yc    += elem.nodes[j]->y/i34[i];
         elem.depth += elem.nodes[j]->z/i34[i];
      }
 
      
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
      const int i34 = elem.i34;
      for (int j = 0; j < i34; j++) {
         const double area = elem.area;
         const int id2 = elem.nxq(1 + j);
         const int id3 = elem.nxq(2 + j);
         const int id4 = elem.nxq(3 + j);
         if (elem.isTri()) {
            // This expression can be derived from the derivatives of the shape functions in area coordinates.
            // a_1 = [(x2 - x)(y1 - y) - (x1 - x)(y2 - y)]/(2*area)
            // da1/dx = (y1 - y2)/(2*area)
            const double dPhi_j_dx = (elem.nodes[id2]->y - elem.nodes[id3]->y)/(2*area); 
            const double dPhi_j_dy = (elem.nodes[id3]->x - elem.nodes[id2]->x)/(2*area); 
            elem.gradPhi[j] = {dPhi_j_dx, dPhi_j_dy};
         }
         else  { // Quad; evaluate at centroid
            const double dPhi_j_dx = (elem.nodes[id2]->y - elem.nodes[id4]->y)/(2*area); 
            const double dPhi_j_dy = (elem.nodes[id4]->x - elem.nodes[id2]->x)/(2*area); 
            elem.gradPhi[j] = {dPhi_j_dx, dPhi_j_dy};
         } // if Tri

      } // for j

   }



   
}


__global__ void setupCoefficientMatrixAndRHS(double time) {
// Setup coefficient matrix, sparsem, for the shallow water wave equation
// No elevation essential b.c. are imposed yet but other b.c. is imposed
// 
// This kernel is called for all nodes.

   int in = 1; //XXX blockIdx.x * blockDim.x + threadIdx.x;

   Node *node1 = &(mesh->localNodeData[in]);


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
      const double hhatb = elem->mean(&Node::hhat);        

      assert((hhatb >= 0.0 || ihh != 1) && "Invalid hhat");

      const int id1 = elem->nodeLocalIdx(node1); 
      const double area = elem->area;

      if (elem->isTri()) {

         // To clarify the notation in this section: 
         //   i standtimes for the node under consideration (node1)
         //   j (node2) and k (node3) are the other nodes of the triangle.
         const int id2 = elem->nxq(1 + id1);
         const int id3 = elem->nxq(2 + id1);
         Node *node2 = elem->nodes[id2];
         Node *node3 = elem->nodes[id3];

         const Point2 x1 = node1->point2();
         const Point2 x2 = node2->point2();
         const Point2 x3 = node3->point2();

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
      if (imm == 2)  { // Pre-compute vnorm on the element
         Vector3 ubed;
         update_bdef(time, elem->xc, elem->yc, elem->depth, ubed);

         Vector3 gradP(0.0, 0.0, 1.0);
         for (int m=0; m < elem->i34; m++) {
           const Vector2 gradP_2D = elem->nodes[m]->depth*elem->gradPhi[m];
           gradP.x += gradP_2D.x;
           gradP.y += gradP_2D.y;
         }   
         vnorm = dot(ubed, gradP)/sqrt(dot(gradP, gradP));
      } // imm==2

      if (elem->isTri()) {
         for (thenint l=0; l<3; l++) {
           double fac;

           Node *node = elem->nodes[l];

           (node == node1)? fac=2.0 : fac=1.0;

           if (imm==2) node->qel += (1.0/12.0)*elem->area*fac*(node->eta2 + dt*vnorm);
           else        node->qel += (1.0/12.0)*elem->area*fac*(node->eta2 + node->bdef2 - node->bdef1);
           
         }

         if (elem->idry == 0) {
           const Vector2 gradEta = elem->nodes[0]->eta2*elem->gradPhi[0] +
                                   elem->nodes[1]->eta2*elem->gradPhi[1] +
                                   elem->nodes[2]->eta2*elem->gradPhi[2];

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
   if (nodeBnd != nullptr)  { 
      // If it is an open boundary node...

      // At least one open bnd side, maybe two
      for (int l=0; l<2; l++) { 

         // Gets each of the sides sharing the node. 
         SideBoundary *sideBnd = nodeBnd->sideBnd(l);

         if (sideBnd == nullptr) continue; // Skip if no boundary side

         const Node *node2 = sideBnd->side->otherNode(node1); // Other node of the side


         if (sideBnd != nullptr) {
            // If the side is on a boundary...

            // Get the boundary segment
            const BoundarySegment *sideSegment = sideBnd->segment;
            const Side *side = sideBnd->side;

            // Integral I_3
            if (sideSegment->ifltype != 0) {
            // Natural or Flather b.c.

               // Calculate I_3 even if i is on essential b.c. so as to check symmetry later
               // especially for Flather b.c.
               double ri3 = 0.0;
               if (side->idry != 1) {
                  // If the side is not dry
                  
                  // Computes the side flux of (uth, vth)
                  double Vn = 0.0;
                  for (int k = side->kb; k<nvrt-1; k++) {
                     // uth, vth in lat/lon frame if ics=2
                     // XXX Are uth,vth really necessary? Why not just store the normal component?
                     const double vn1 = dot(side->uth[k  ], side->n); // outer normal
                     const double vn2 = dot(side->uth[k+1], side->n);

                     Vn += 0.5*(side->zc[k+1] - side->zc[k])*(vn1 + vn2);
                  }

                  ri3 = 0.5*side->length*Vn;
               }

               if (sideSegment->ifltype == -1) { // Flather 1

                  assert(node1->eta2 >= -98.0 && node2->eta2 >= -98.0 && "Mismatch 1") ;
                  assert(side->depth > 0.0 && "Negative depth at Flather boundary");

                  const double con0 = (1.0/6.0)*side->length*sqrt(grav*side->depth); // For coefficient matrix
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

                  const double con0 = (1.0/6.0)*side->length*tmp*thetai*dt;
                  node1->addToSparseMatrix(node1, con0*2.0);
                  node1->addToSparseMatrix(node2, con0);

               } // if Flather or Discharge

               node1->qel -= thetai*dt*ri3;

            } // I_3
          

            // Integral I_5
            if(side->idry != 1) {

               const double Unbar = dot(side->U, side->n);

               node1->qel -= 0.5*(1.0 - thetai)*dt*side->length*Unbar;
               
            } // I_5

         } // if is boundary side

      } // for l

   } // if isbnd

}



__device__ void backtrackAllNodes(void) {

   // max/min # of sub-divisions
   const double ndelt_max = max(1.0, dt/dtb_min);
   const double ndelt_min = max(1.0, dt/dtb_max); 


   swild98=0.d0 !init

   do i=1,ns
        sdbt(1,:,i)=su2(:,i)
        sdbt(2,:,i)=sv2(:,i)
   enddo !i

   side = FROM kernel

   if (side->idry == 1) return;

   isd0=i
   jmin=kbs(isd0)
   ie0=0

   // Checks that at least one element is not dry
   const Element *elem1 = side->elements[0];
   const Element *elem2 = side->elements[1];
   assert ((elem1 != nullptr && !(elem1->idry)) || (elem2 != nullptr && !(elem2->idry)) && "backtrackAllNodes: non dry side without two dry elements.");


   !swild_tmp to store global coord. of the starting side
   swild_tmp(1)=xcj(isd0); swild_tmp(2)=ycj(isd0); swild_tmp(3)=zcj(isd0)
   !swild10_tmp to store frame at starting pt for ics=2 (not used for ics=1)
   !Use sframe2 
   swild10_tmp(1:3,1:3)=sframe2(:,:,i) !pframe(:,:,isidenode(1,isd0)) 


   for (int j = jmin; j < nvrt; j++) {
   //  Initialize (xt,yt,zt), nnel and vel.
   //  For ics=2, the coordinates are in local frames
   //  Caution! nnel must be initialized inside this loop as it is updated inside.

      ipsgb=islg(isd0)

      Node *node1 = side->nodes[0];
      Node *node2 = side->nodes[1];

      iadvf = min(iadv(n1),iadv(n2))

      Point3 xt;
      if constexpr(ics==1) {
         xt = side->center(j); 
      }
      else // lat/lon; in side lat/lon frame
      {
      // XXX IMPLEMENT THIS
   //          xt=0.d0
   //          yt=0.d0
   //          !centroid coord. for nudging
   //          call project_pt('g2l',xctr(ie0),yctr(ie0),zctr(ie0), &
   //   &(/xcj(isd0),ycj(isd0),zcj(isd0)/),swild10_tmp,xctr2,yctr2,tmp)
      }  // ics


      // Velocity at the side (in ll frame for ics, same vertical direction)
      // XXX IS THIS NECESSARY?
      Vector3 u_int = side->velocity(j); //{side->u2[j].x, side->u2[j].y, 0.5*(side->w2[j] + side->w2[j+1])};

      // If on land or open boundary...
      const bool ifl_bnd =  (side->sideBnd != nullptr);
      
      const double vmag = length(toVector2(u_int));

      nnel = ie0;

      vis_coe = 0.0;

      if (vmag <= velmin_btrack) {
         // No activity 
         // XXX IS THIS NECESSARY? 

         sdbt(1,j,isd0)=su2(j,isd0)
         sdbt(2,j,isd0)=sv2(j,isd0)
         if(ielm_transport/=0) sdbt(3:2+ntracers,j,isd0)=(tr_nd(1:ntracers,j,n1)+tr_nd(1:ntracers,j,n2))*0.5d0
         swild98(1:2,j,isd0)=su2(j,isd0) !max/min
         swild98(3:4,j,isd0)=sv2(j,isd0)


      } else {  // do btrack

         // Compute number of sub-divisions based on local gradients 
         // XXX IS THIS NECESSARY?
         double suma = 0.0;
         int icount = 0;
         for (int ii=0; ii<2; ii++) {

            const Element *elem = side->elements[ii];
              
            if (elem == nullptr) continue;

            icount=icount+1

            // Not strictly along z; in ll frame for ics=2

            for (int jj = 0; jj < elem->i34; jj++) {

               const Node *node = elem->nodes[jj];

               if constexpr(ics==1) {
                  Vector2 u_node = {node->u2[j].x, node->u2[j].y};
               } else {
                  // XXX IMPLEMENT THIS
                  // call project_hvec(uu2(j,nd0),vv2(j,nd0),pframe(:,:,nd0),eframe(:,:,ie),utmp0(jj),vtmp0(jj))
               }
            } // for jj

            const Vector2 grad_u = elem->grad(u_int.x); // XXX IMPLEMENT THIS
            const Vector2 grad_v = elem->grad(u_int.y);
            // dudx = dot_product(utmp0(1:i34(ie)),dldxy(1:i34(ie),1,ie))
            // dudy = dot_product(utmp0(1:i34(ie)),dldxy(1:i34(ie),2,ie))
            // dvdx = dot_product(vtmp0(1:i34(ie)),dldxy(1:i34(ie),1,ie))
            // dvdy = dot_product(vtmp0(1:i34(ie)),dldxy(1:i34(ie),2,ie))

            suma += dt*sqrt(dot(grad_u, grad_u) + dot(grad_v, grad_v));
         } // for ii

         assert(icount > 0 && "Impossible 77");

         // '4' is somewhat arbitrary
         ndelt = max0(ndelt_min,min0(ndelt_max,int(suma/icount)*4)) // >=1
         dtbk=dt/ndelt // Target btrack step; may be smaller sometimes

         // Perturb starting pt to avoid underflow 
         // XXX I don't think this is needed with the current implementation
         // eps=btrack_nudge
         // if constexpr(ics==1) {
         //    xt=(1.d0-eps)*xt+eps*xctr(nnel)
         //    yt=(1.d0-eps)*yt+eps*yctr(nnel)
         // } else { //lat/lon
         //    xt=(1.d0-eps)*xt+eps*xctr2
         //    yt=(1.d0-eps)*yt+eps*yctr2
         // } // ics

         time_rm=dt
         
         // call btrack(l,ipsgb,ifl_bnd,j,iadvf,swild_tmp,swild10_tmp, dtbk,vis_coe,time_rm,time_rm2,uuint,vvint,wwint,nnel,jlev,xt,yt,zt,swild3,ltmp)

            
      } // if vmag (do backtrack)

   } // for j=jmin,nvrt-1

}

template<int nc>
__device__ void tridag_sch(int n, double *a, double *b, double *c, double *rhs, double *u, double *gam) {

// Solves a tridiagonal system. It was adapted from "Numerical Recipes in FORTRAN (pp.43 ).
//
// a,b,c,r: input vectors and are not modified by this program.
// b is the main diagonal, a below and c above. a(1) and c(n) are not used.
// r is the r.h.s.
// n: actual rank of the system.
// rows in u() and r();
// u: output with nc columns (transposed to rows)
// gam: a working array.
//
// NOTES:
// This is a template function that allows for different numbers of columns in the rhs vector.
// [n, nc] is the shape of Fortran-type arrays rhs and u.

   assert(n >= 1       && "TRIDAG: n must be >= 1");
   assert(nc <= nvec   && "TRIDAG: Increase number of columns");
   assert(b[0] != 0.d0 && "TRIDAG: b[0]=0");

   double bet = b[0];
   if constexpr (nc == 1) {
      u[0] = rhs[0]/bet;
   } else {
      for (int k = 0; k < nc; k++) u[k*nc] = rhs[k*nc]/bet;
   }

   for (int j = 1; j < n; j++) {
      gam[j] = c[j-1]/bet;
      bet = b[j] - a[j]*gam[j];

      assert (bet!=0.d0 && "TRIDIAGONAL SOLVER: bet=0");

      if constexpr (nc == 1) {
         u[j] = (rhs[j] - a[j]*u[j - 1])/bet;
      } else {
         for (int k = 0; k < nc; k++) {
            u[j + k*nc] = (rhs[j + k*nc] - a[j]*u[j - 1 + k*nc])/bet;
         }
      }
   }

   // Back-substitution
   for (int j=n-1; j>=1; j--) {
      if constexpr (nc == 1) {
         u[j] = u[j] - gam[j+1]*u[j-1];
      } else {
         for (int k = 0; k < nc; k++) {
            u[j + k*nc] = u[j + k*nc] - gam[j + 1]*u[j-1 + k*nc];
         }
      }
   }
}


__device__ void cubicSplineDerivatives(int n, double *x, double *y, double *y_top, double *y_bot, double *yp_out, double *ypp_out) {
// Compute the 1st and 2nd derivatives of y = f(x) using continuous cubic splines.
// 
// This function uses the General expression for a C2 interpolating cubic spline (search for "C2 interpolating cubic spline" in Wikipedia)
// This method uses the 2nd derivatives and the y values at the endpoints to generate the polynomial coefficients for each level.
//
// NOTE: (x,y) here are just any function y=f(x). In fact we use these splines for vertical interpolation.
// 
// Inputs: 
//            n: dimesion for x etc; # of data pts (>=2);
//            x[n]: x coordinates; must be in ascending order (and distinctive);
//            y[n]: functional values; 
//            y_top and y_bot: 1st derivatives at x[0] and x[n-1];
//     Output: 
//            ypp[n]: 2nd derivatives used in interpolation.
//             yp[n]: 1st derivatives 

   assert(n >= 2 && "CUBIC SPLINE: n < 2");

   double a[n];   // Lower diagonal
   double b[n];   // Diagonal
   double c[n];   // Upper diagonal
   double gam[n]; // Gamma, used by the tridiagonal solver
   double rhs[n]; // Right-hand side

   for (int k=0; k<n; k++) {

      if (k == 0) {
         // Bottom level
         assert(x[k+1] - x[k] != 0.0 && "CUBIC SPLINE: bottom problem");

         b[k] = (x[k+1] - x[k])/3.0;
         c[k] = b[k]/2.0;
         rhs[k] = (y[k+1] - y[k])/(x[k+1] - x[k]) - yp1;
         
      } else if (k == npts - 1) {
         // Top level (surface)
         assert(x[k] - x[k-1] != 0.0 && "CUBIC SPLINE: surface problem");

         b[k] = (x[k] - x[k-1])/3.0;
         a[k] = b[k]/2.0;
         rhs[k] = -(y[k] - y[k-1])/(x[k] - x[k-1]) + yp2;
        
      } else {
         // Intermediate levels
         assert(x[k] - x[k-1] != 0.0 && x[k+1] - x[k] != 0.0 && "CUBIC SPLINE: middle problem");

         b[k] = (x[k+1] - x[k-1])/3.0;
         a[k] = (x[k  ] - x[k-1])/6.0;
         c[k] = (x[k+1] - x[k  ])/6.0;

         rhs[k] = (y[k+1] - y[k])/(x[k+1] - x[k]) - (y[k] - y[k-1])/(x[k] - x[k-1]);

      }
   }

   // Solve the tridiagonal system
   tridag_sch<1>(n, 1, a, b, c, rhs, ypp, gam)
    
   yp[    0] = yp1; 
   yp[n - 1] = yp2;

   // Compute rest of 1st derivatives using the splines
   for (int k=1; k < n - 1; k++) {
      yp[k] = (y[k+1] - y[k])/(x[k+1] - x[k]) - (x[k+1] - x[k])/6.0*(2.0*ypp[k] + ypp[k+1]);
   } 
}



__device__ void computeContinuityEquationRHS(void) {
// Continuity equation: preparation of matrix

   // Compute elevation essential boundary conditions
   // In case of border node on >1 bnd with imposed elevation, the bnd with largest segment # prevails.
   double elbc[nope_global]; // flags
   for (int i=0; i < nope_global; i++) {
      // For all boundary segments...

      for (int j=0; j < nond_global(i); j++) {
         // For all nodes in the segment...

         const Node *node = iond_global(i,j); //global

         if (ipgl(nd)%rank==myrank) then

            ip=ipgl(nd)%id
            if(nramp_elev==1) then
              eic=etaic(ip)
            else
              eic=0.d0
            endif

            // Prep tide
            if (iettype(i)==3 || iettype(i)==5) then
              eta1_bar=(1.d0-ramp)*eic !etaic(ip) !initialize
              do jfr=1,nbfr
                ncyc=int(amig(jfr)*time/2.d0/pi)
                arg=amig(jfr)*time-real(ncyc,rkind)*2.d0*pi+face(jfr)-efa(i,j,jfr)
                eta1_bar=eta1_bar+ramp*ff(jfr)*emo(i,j,jfr)*cos(arg)
              enddo !jfr=1,nbfr
            endif

            if(iettype(i)==1.or.iettype(i)==2) then
              elbc(ip)=ramp*eth(1,i)+(1.d0-ramp)*eic
            else if(iettype(i)==3) then
              elbc(ip)=eta1_bar
            else if(iettype(i)==4) then
              elbc(ip)=ramp*eth(j,i)+(1.d0-ramp)*eic !etaic(ip)
            else if(iettype(i)==5) then
              elbc(ip)=ramp*eth(j,i)+eta1_bar
            endif

            !Add inverse barometric effects
            if(iettype(i)/=0.and.inv_atm_bnd==1) elbc(ip)=elbc(ip)+ramp*(prmsl_ref-pr(ip))/grav/rho0

          endif !ipgl(nd)%rank==myrank
        enddo !j
   } // for i


   side->U=0.d0; 

   // Computes the vertically integrated velocity U.
   Side *side = FROM kernel;
      
   if(side->idry) return;


   const Node *node1 = side->n1;
   const Node *node2 = side->n2;


   // Vertically integrated velocity U (in ll if ics=2)
   side->U = {0.0, 0.0};
   for (int k=side->kbs; k<nvrt-1; k++) { /// XXX CHECK INDICES
      side->U += side->Dz(k+1)*0.5*(side->u2[k] + side->u2[k+1]);  // CHECK INDICES!!!!
   }

}


__device__ void computeVegetationModelParameters(void) {
   
   Side *side = FROM kernel;

      
   chi=0.d0; chi2=0.d0; hhat=0.d0; 
   hhat2=0
   veg_c2=0.d0; veg_beta=0.d0


   const double htot = side->htot();

   assert(htot > 0.0 && "computeVegetationModelParameters: htot <= 0");

   const double veg_h  = side->veg_h();
   const double veg_nv = side->veg_nv();


   // Compute the drag coefficients chi2 and chi (\tilde\chi in notes)
   const double vmag1 = module(side->u2[side->kbs + 1])
   side->chi2 = side->Cd*vmag1 
   side->chi  = side->chi2
   if constexpr(iveg != 0) {
      veg_alpha_bot = (veg_alpha3D(kbp(n1),n1)+veg_alpha3D(kbp(n2),n2))/2.d0
      side->chi = side->chi/(1.0 + veg_alpha_sd_bot*vmag1*dt)
   }

   // Calc consts in SAV model; make sure veg_[c2,beta]=0 at 2D, dry side, and emergent side
   // uuint=0.00 !\int_{-h}^{z_v} |u|dz / H^\alpha = \bar{|u|}^\alpha
         if (iveg != 0 && nvrt - kbs(i) > 1 && veg_h_sd > 0.0) {
            // 3D wet side with veg
            const double zctr2 = side->zTopCanopy(); // top of canopy
            for (int k = side->kbs; k < nvrt - 1; k++) {
               if(zctr2>zs(k,i)) then !partial or full SAV layer
               zz1=min(zctr2,zs(k+1,i))
               zrat=(zz1-zs(k,i))/(zs(k+1,i)-zs(k,i)) !\in (0,1]
               !if(zrat<=0.or.zrat>1) call parallel_abort('STEP: WOW2')
               ub2=(1.d0-zrat)*su2(k,i)+zrat*su2(k+1,i) !@upper level
               vb2=(1.d0-zrat)*sv2(k,i)+zrat*sv2(k+1,i) 
   !              bigu1=bigu1+(zz1-zs(k,i))*(su2(k,i)+ub2)/2
   !              bigv1=bigv1+(zz1-zs(k,i))*(sv2(k,i)+vb2)/2
               ubar1=sqrt(su2(k,i)*su2(k,i)+sv2(k,i)*sv2(k,i))
               ubar2=sqrt(ub2*ub2+vb2*vb2)
               uuint=uuint+(zz1-zs(k,i))*(ubar1+ubar2)/2.d0
               } // if
            } // for
      uuint=uuint/min(htot,veg_h_sd) !>=0

      veg_alpha     = sum(veg_alpha_vert_mean(isidenode(1:2,i)))/2.d0
      side->veg_c2 = side->veg_alpha*dt*uuint; // >=0
      if(side->veg_h < 0.99*htot) { // 3D wet submergent side with vegetation
      tmpx2 = max(1.0e-5, chi2(i)*vmag1/grav/htot) // Energy gradient
      // \beta_2; arguments checked
      tmpy2=sqrt(sqrt(veg_nv_sd)/veg_h_sd)*(-0.32d0-0.85d0*log10((htot-veg_h_sd)/veg_h_sd*tmpx2))
      veg_beta(i)=exp(tmpy2*(zctr2-zs(kbs(i)+1,i)))-1.d0 !\beta
      veg_beta(i)=min(10.d0,max(veg_beta(i),0.d0))
      endif !veg_h_sd
   endif !iveg

      // hhat is \breve{H} in notes
      if (side->is2D()) { // 2D

         double tmp = htot + side->chi2*dt;
         if (iveg != 0) tmp = tmp + veg_alpha_sd*vmag1*htot*dt;

         assert(tmp > 0.0 && "COMPUTE CONTINUITY EQUATION: Impossible dry 53");

         hhat(i) = htot*htot/tmp; 

      else { // 3D
         side->hhat = htot - side->chi*dt;
!        hhat(i)=hhat2(i)

         if constexpr(iveg != 0) {
            const double veg_h = side->veg_h();
            if (veg_h < 0.99*htot) { // Submergent
               const double veg_c    = side->veg_c();
               const double veg_beta = side->veg_beta();
               side->hhat -= veg_c*(veg_h + veg_beta*side->chi*dt);

            } else { // Emergent
               side->hhat = side->hhat/(1.0 + veg_c2); // veg_alpha_sd*uuint*dt)
            }
         } // if iveg

         // Enforce positivity for 3D model
         if (ihhat == 1) side->hhat = max(0.0, side->hhat)
      } // if 3D

        
      // Baroclinic force at side and whole levels

      bcc=0 !init for 2D cases

}

#ifdef DONTCOMPILE
__device__ void computeVerticalAdvection(void) {
   // Compute the vertical advection term d(uw)/dz and d(vw)/dz, using cubic splines.

   sdbt(1:2,:,:)=0.d0

   const Side *side = FROM kernel;

   if(side->idry) return;

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
   cubicSplineDerivatives(kbs, nvrt, side, ws, up)

   double top_dvw_dz = (ws[nvrt - 1]*side->v2[nvrt] - ws[nvrt - 1]*side->v2[nvrt - 1])/(side->z[nvrt - 1] - side->z[nvrt - 2])
   double bot_dvw_dz = (ws[kbs     ]*side->v2[kbs ] - ws[kbs     ]*side->v2[kbs  - 1])/(side->z[kbs ] - side->z[kbs  - 1])
   cubicSplineDerivatives(kbs, nvrt, side, ws, vp)

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
!      call exchange_s3d_4(sdbt)
      allocate(swild96(2,nvrt,nsa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: fail to allocate swild96 (3.4)')
      swild96=sdbt(1:2,:,:)
      call exchange_s3d_2(swild96)
      sdbt(1:2,:,:)=swild96
      deallocate(swild96)

#ifdef INCLUDE_TIMING
      wtimer(4,2)=wtimer(4,2)+mpi_wtime()-cwtmp
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      endif !ELM or upwind

#ifdef TIMER2
      tmp=mpi_wtime()
      write(12,*)'Time taken for mom advection=',tmp-cwtmp3,it
      cwtmp3=tmp !reset
#endif

      if(itransport_only==0) then

}



__device__ void computeElemAverageGhat1(void) {
// Element average ghat1 (in eframe if ics=2; \breve{G} in notes). Dimensions: m^2/s.

   Element *elem = FROM kernel;
   
   elem->ghat1= {0.0, 0.0};  

   if constexpr(ihydraulics != 0) {
      if (nhtblocks > 0 && (elem->isblock > 0)) return;
   }

   if (elem->idry == 1) return;

   // Wet elements.
   // Warning: ghat1=elem average of \breve{G} must include all: Coriolis, atmo. pressure, 
   // tidal potential, horizontal diffusion, and baroclinic etc
   // Remember to update both fb and F 
   // If ics=2, ghat1 are in eframe (ll)
   // cff1-3: coefficients for 2D/3D
   Vector2 gradEta = {0.0, 0.0};  // \nabla \eta
   Vector2 gradPr  = {0.0, 0.0};  // \nabla pressure
   Vector2 gradEtp = {0.0, 0.0};  // \nabla etp

   double av_cff1 = 0.0;  // av. for cff1
   double av_cff2 = 0.0;  // av. for cff2
   double av_cff3 = 0.0;  // av. for cff3
   double av_df = 0.0;  // av. for H^\alpha / H
   double av_cff2_chi = 0.0;  // av. for cff2*chi
   double av_cff3_chi = 0.0;  // av. for cff3*chi

   const int i34 = elem->i34;

   for (int j=1; j<i34; j++) {
      // For each node and side in the element...

      const Node *node = elem->nodes[j];
      const Side *side = elem->sides[j];

      const double chi      = side->chi;
      const double veg_c2   = side->veg_c2;
      const double veg_beta = side->veg_beta;
      const double veg_h    = side->veg_h();
      const double veg_c    = veg_c2/(1.0 + veg_c2);  // veg_c2 = \alpha*|u|*dt

      // Node - const in elem
      gradEta += node->eta2*elem->gradPhi[j];  // in eframe if ics=2
      gradPr  += node->pr  *elem->gradPhi[j];

      if (elem->depth >= tip_dp) {
        gradEtp += node->etp*elem->gradPhi[j];
      }

      
      const double htot  = side->htot();        // Total depth at the side
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
      for (int k = side->kb; k < nvrt; k++) {
         if (side->is2D()) {  //2D; no need to compute U^{\star\alpha}
            Ustar += side->bt[nvrt]*side->Dz(k);

         } else { //3D
            Vector2 w2 = (side->bt[k] + side->bt[k-1])/2.0*side->Dz(k);
            Ustar += w2;

            if constexpr(iveg != 0) {
                if (zctr2 > side->zc[k-1]) {  // XXX CHECK INDICES!!! 
                    const double zrat = min(1.0, (zctr2 - side->zc[k-1])/side->Dz(k)); // For the top of the canopy, ratio of the level that is below the canopy
                    Ustar_alpha += zrat*w2; // XXXX
                }
            }  // if iveg
         }  //2/3D
      }  //k
      

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
         const Vector2 ubstar = side->sdbt[side->kb+1]; // XXX CHECK INDICES
      } else {
         // XXX IMPLEMENT THIS
         // n1=isidenode(1,isd); n2=isidenode(2,isd)
         // call project_hvec(tau(1,n1),tau(2,n1),pframe(:,:,n1),eframe(:,:,i),vn1,vn2)
         // call project_hvec(tau(1,n2),tau(2,n2),pframe(:,:,n2),eframe(:,:,i),tt1,ss1)
         // tau_x=(vn1+tt1)/2.d0
         // tau_y=(vn2+ss1)/2.d0
         // call project_hvec(sdbt(1,kbs(isd)+1,isd),sdbt(2,kbs(isd)+1,isd),sframe2(:,:,isd),eframe(:,:,i),ubstar,vbstar)
      }  // if ics

      // Term for surface and bottom stress
      elem->ghat1 -= chi*dt*cff2*ubstar + dt*cff1*tau;

      // All terms in F,F^\alpha,f_b except baroclinic
      // F, f_b: all explicit terms excluding vegetation 
      // F^\alpha=\int_{-h}^z_v f dz (f ncludes all explicit terms except for vegetation)
      // Only need to calculate F^\alpha for 3D sides (but init=0) 
      const double coriolis = side->cori;
      const Vector2 F  =  coriolis*side->U;  // F_x
      const Vector2 fb =  coriolis*side->u2;  // f_b_x
      Vector2 Fa = {0.0, 0.0};  // F^\alpha, vertically forcing terms baroclinic gradient, horizontal viscosity, Coriolis, earth tidal potential, atmospheric pressure, radiation stress)

      if (iveg != 0 && nvrt > side->kb) {  // 3D
         
         for (int k = side->kb; k < nvrt; k++) {
            if (zctr2 > side->z[k-1]) {
                const Vector2 w2 = side->Dz[k]*(side->u2[k] + side->u2[k-1])/2.0;
                const double zrat = min(1.0, (zctr2-side->z[k-1])/side->Dz[k]);
                Fa +=  coriolis*w2*zrat;
            } // if zctr2
         } // k
      } // iveg

      // Integrate horizontal viscosity vertically
      for (int k=side->kb; k<nvrt; k++) {
         const Vector2 levelHorVisc = side->Dz(k)*(side->d2uv[k] + side->d2uv[k-1])/2.d0 // (d2uv(1,k,isd)+d2uv(1,k-1,isd))/2

         F += w2;

         // Vegetation term
         if constexpr(iveg!=0) {
            if (nvrt > side->kbs + 1 && zctr2 > side->z[k-1]) { // 3D
               const double zrat = min(1.0, (zctr2 - side->z[k-1])/side->Dz(k));
               Fa += w2*zrat;
            } 
         } // if iveg
      } // for k

      fb += side->d2uv[side->kb];


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
   const double grav3 = elem->mean(Element::grav2)
   const Vector2 fb = grav3*gradP - gradPr/rho0 //const in each elem


   tmp = {0.0, 0.0}; //elem average of all terms; into ghat1
   for (int j=1; j<i34; j++) { 
      // For each side in the element...
      Side *side = elem->sides[j];

      const double htot = side->htot();
      const double F  = htot*fb;
      const double Fa = min(htot, veg_h)*fb;

      tmp += av_cff1*dt*F - av_cff2*chi*dt*dt*fb - av_cff3*dt*Fa;
   } // for j

   //botf1 etc are already in elem frame
   elem->ghat1 += tmp/i34;

   // Baroclinic force
   if constexpr(ibc==0) {
      assert(prho(1,elnode(1,i))>= -98.d0 && prho(1,elnode(2,i))>= -98.d0 && prho(1,elnode(3,i))>= -98.d0 && "Impossible dry 5");
      
      // swild2(k,:) = \sum_{l=k}^N dr*dz; whole level (and eframe if ics=2)
      swild2[nvrt-1,1:2]=0.d0
      for (int k=nvrt-2; k>elem->kbe; k--) { /// XXXX CHECK INDICES
         // dr_dxy is grad(rho)
         swild2(k,1:2)=swild2(k+1,1:2)+dr_dxy(1:2,k+1,i)*(ze(k+1,i)-ze(k,i)) 
      }
      swild(1:2)=0.d0 !\in F [m^2/s/s]
      for (int k=side->kb; k<nvrt; k++) {
         swild(1:2) -= grav3/rho0*(side->Dz(k)*(2.d0*swild2(k,1:2)+dr_dxy(1:2,k,i)*side->Dz(k))
      }
      fb = -grav3/rho0*swild2(kbe(i)+1,1:2) //[m/s/s]; elem average
      Fa = av_df*swild(1:2) //approx
      elem->ghat1 += rampbc*dt*av_cff1*swild(1) - rampbc*dt*dt*av_cff2_chi*fb- rampbc*av_cff3*dt*Fa;

   } // if ibc

}



#ifdef DONTCOMPILE
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
                           int np, int ne, int ns, int nvrt, int kz,
                           int max_np, int max_ne, int max_ns,
                           const Node **nodesPtr, const Element **elemsPtr, const Side **sidesPtr,
                           Node *nodes, Element *elems, Side *sides,
                           Node **partNodeNeighNodesMem, Element **partNodeNeighElemsMem,
                           NodeBoundary *localNodeBoundaryMem, SideBoundary *localSideBoundaryMem, BoundarySegment *localSegmentMem,
                           double h_s, double theta_b, double theta_f, double *sigma_lcl, double *sigma, 
                           double *ztot, double *cs, double *dcs) {

   if (mesh != nullptr) 
   {
      printf("ERROR: Mesh is already initialized\n");
      return;
   };

   printf("Creating mesh for partition %d\n", partNum);
   mesh = new Mesh(partNum, npgb, negb, nsgb, nvrt, kz, np, ne, ns, nseg, max_np, max_ne, max_ns,
                   nodesPtr, elemsPtr, sidesPtr, nodes, elems, sides,
                   partNodeNeighNodesMem, partNodeNeighElemsMem, 
                   localNodeBoundaryMem, localSideBoundaryMem, localSegmentMem,
                   h_s, theta_b, theta_f, sigma_lcl, sigma, ztot, cs, dcs);

   // XXX read from file
   polyDischarge = new Polynomial<4>(0.0, 0.0, 0.0, 0.0, 0.0);                   
   
}

extern "C"
  __global__ void populateMesh(int *i34, int *elnode, int *elside,
                   double *node_coords, double *wind_rotate_angle, 
                   int *isdel, int *isidenode, int *isbnd,
                   int *nond, int *nlnd, int *iond, int *ilnd,
                   int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                   int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                   int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                   int *localNodeIndices, int *localElemIndices, int *localSideIndices,
                   int *localNodePart,    int *localElemPart,    int *localSidePart) {

    mesh->populateMesh(i34, elnode, elside, node_coords, wind_rotate_angle, 
                       isdel, isidenode, isbnd,
                       nond, nlnd, iond, ilnd,
                       nodeNeighNodesStart, nodeNeighNodesIndices,
                       nodeNeighElemsStart, nodeNeighElemsIndices,
                       globalNodeIndices, globalElemIndices, globalSideIndices,
                       localNodeIndices, localElemIndices, localSideIndices,
                       localNodePart, localElemPart, localSidePart);
  }

extern "C" 
__global__ void destroyMeshCPU(void) {
   delete mesh;
}

extern "C"
__global__ void test(double x, double y, double z) {
   mesh->test(x, y, z);
}