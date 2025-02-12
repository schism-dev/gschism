


// Each device has a global mesh object
__device__ Mesh* mesh = nullptr;





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
   @param nextElement: element that contains the intersection point (input and output)
   @param tIntersect: New intersection point tIntersect*side.n1 + (1-tIntersect)*side.n2

   The equations for the two segments (the one we are considering and one of the sides of the element) are:
   P1(tt1) = (x1,y1) + tt1*(x2-x1, y2-y1)  where 0 <= tt1 <= 1
   P2(tt2) = (x3,y3) + tt2*(x4-x3, y4-y3)  where 0 <= tt2 <= 1

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
__device__ bool Element::findNextIntersection(double x2, double y2, 
                                              Side *&side, Element *&nextElement,
                                              double &t, double& tIntersect) {


   double closest_intersect_t = 1e10;
   bool found_intersection = false;

   // (x1,y1)-(x2,y2) is the track line segment. Every step we move (x1,y1) to the next
   // intersection point, at position t on "side". The final position is (x2,y2)=(xt,yt).
   const double x1 = (1-t)*side->n1->x + t*side->n2->x;
   const double y1 = (1-t)*side->n1->y + t*side->n2->y;
   printf(">findNextIntersection: x1=%f, y1=%f, x2=%f, y2=%f\n", x1, y1, x2, y2);

   Side *nextSide = nullptr;

   for (int i = 0; i < i34; i++) {

      // Element's i-th side 
      Side *side2 = sides[i];

      printf("SIDE: side = %p, side2 = %p\n", side, side2);
      // Skip the side that contains the initial point
      if (side2 == side) continue;

      printf("### i = %d, side2->n1->x=%f, side2->n1->y=%f, side2->n2->x=%f, side2->n2->y=%f\n", i, side2->n1->x, side2->n1->y, side2->n2->x, side2->n2->y);

      // (x3,y3)-(x4,y4) is the element side line
      const double x3 = side2->n1->x;
      const double y3 = side2->n1->y;
      const double x4 = side2->n2->x;
      const double y4 = side2->n2->y;

      // delta is 0 if the two lines are parallel. In this case, they do not intersect.
      const double delta = (x2-x1)*(y3-y4) - (y2-y1)*(x3-x4);

      if (delta != 0.0) {

         const double delta1 = (x3-x1)*(y3-y4) - (y3-y1)*(x3-x4);
         const double delta2 = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1); 

         // Solves the two linear equations for tt1 and tt2
         const double tt1 = delta1/delta; // delta is not zero, so we can divide
         const double tt2 = delta2/delta;

         if (tt1 > 0.0 && tt2 >= -small && tt2 <= 1.0 + small) {
            closest_intersect_t = min(closest_intersect_t, tt1);
            nextSide = side2;
            nextElement = nextSide->otherElement(this);
            tIntersect = max(small, min(tt1, 1.0-small));
            t = tt2;
            printf("### t = %f, tIntersect=%f\n", t , tIntersect);
            found_intersection = true;
         }
      
      }

   }

   if (!found_intersection) {
      printf("### findNextIntersection: No intersection found, tIntersect=%f\n", tIntersect);
   }

   
   // found_intersection should always be true, since we checked all the sides.  
   side = nextSide;
   return found_intersection;
}


__device__ bool Element::findNextIntersection(double x1, double y1, double x2, double y2, 
                                              Side *&side, Element *&nextElement,
                                              double &t, double& tIntersect) {

   double closest_intersect_t = 1e10;
   bool found_intersection = false;

   printf("findNextIntersection: x1=%f, y1=%f, x2=%f, y2=%f\n", x1, y1, x2, y2);

   for (int i = 0; i < i34; i++) {

      Side *side2 = sides[i];

      // (x3,y3)-(x4,y4) is the element side line
      const double x3 = side2->n1->x;
      const double y3 = side2->n1->y;
      const double x4 = side2->n2->x;
      const double y4 = side2->n2->y;

      // delta is 0 if the two lines are parallel. In this case, they do not intersect.
      const double delta = (x2-x1)*(y3-y4) - (y2-y1)*(x3-x4);

      if (delta != 0.0) {

         const double delta1 = (x3-x1)*(y3-y4) - (y3-y1)*(x3-x4);
         const double delta2 = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1); 

         // Solves the two linear equations for tt1 and tt2
         const double tt1 = delta1/delta; // delta is not zero, so we can divide
         const double tt2 = delta2/delta;

         if (tt1 > 0.0 && tt2 >= -small && tt2 <= 1.0 + small) {
            closest_intersect_t = min(closest_intersect_t, tt1);
            side = side2;
            nextElement = side->otherElement(this);
            t = tt2;
            tIntersect = max(small, min(tt1, 1.0-small));
            found_intersection = true;
            
         }
      
      }

   }

   // found_intersection should always be true, since we checked all the sides.  
   return found_intersection;
}


__inline__ __device__ double Element::computeSignedArea(void) const {
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



/**
   Compute the shape coordinates of the point (xt,yt) in the element.
   @param xt, yt: coordinates of the point
   @param shpCoords: shape coordinates of the point
   @param fixNegShapeCoords: if true, fix negative shape coordinates
   @return true if the point is inside the element, false otherwise
**/
__device__ bool Element::computeShapeCoords(double &xt, double &yt, double shpCoords[4],
                                            bool fixNegShapeCoords) const {

   // Compute the shape coordinates of the point (xt,yt) in the element.
   if (i34==3) {
      // Compute area coordinates for a given pt w.r.t. to a triangular element
      // If ifl=1, will fix 0 or negative area coord. (assuming it's not too negative)
      // and in this case, the pt will be nudged into the element
      // It's not reliable to use the area coord for ics=2 when the pt is far away
      // from the local frame, so use additional check.

      double x0, y0, x1, y1, x2, y2;

      if constexpr(ics==1) {
         x0 = nodes[0]->x;
         y0 = nodes[0]->y;
         x1 = nodes[1]->x;
         y1 = nodes[1]->y;
         x2 = nodes[2]->x;
         y2 = nodes[2]->y;
      }
      else // lat/lon
      {
#ifdef ICS2
         globalToLocal(nodes[0]->x, nodes[0]->y, nodes[0]->z, x0, y0, z0);
         globalToLocal(nodes[1]->x, nodes[1]->y, nodes[1]->z, x1, y1, z1);
         globalToLocal(nodes[2]->x, nodes[2]->y, nodes[2]->z, x2, y2, z2);
#endif
      }

      shpCoords[0] = signedArea(xt,x1,x2,yt,y1,y2)/area;
      shpCoords[1] = signedArea(x0,xt,x2,y0,yt,y2)/area; 
      shpCoords[2] = 1.0 - shpCoords[0] - shpCoords[1];
      double tmpmin  = min3(shpCoords[0], shpCoords[1], shpCoords[2]);

      // If at least one of the shape coordinates is negative, and we have the flag fixNegShapeCoords set,
      // then we fix the negative shape coordinates.
      if (fixNegShapeCoords && tmpmin<=0.0) {

         // XXX I have changed this algorithm. Now it scales and offsets all the shape coordinates by the same factor,
         // keeping the sum of the shape coordinates equal to 1.
         // I have tested and compared the new algorithm with the old one, and they give
         // very similar results in term of how close the "fixed" point is to the original point in cartesian coordinates.
         // The new algorithm is faster since it is branchless.
         const double a = - tmpmin + 1.0e-2;  // XXX it used to be 1.0e-4 in the comments of the original code
         const double b = 1.0/(shpCoords[0] + shpCoords[1] + shpCoords[2] + 3.0*a);
         
         shpCoords[0] = (shpCoords[0] + a)*b;
         shpCoords[1] = (shpCoords[1] + a)*b;
         shpCoords[2] = (shpCoords[2] + a)*b;          

         if (shpCoords[0] < 0.0 || shpCoords[1] < 0.0 || shpCoords[2] < 0.0 ||
             shpCoords[0] > 1.0 || shpCoords[1] > 1.0 || shpCoords[2] > 1.0) {
            Error("Negative shape coordinates in computeShapeCoords, when fixNegShapeCoords is true");
         }

         // Update the point coordinates for the new shape coordinates
         double xcoords[3] = {x0, x1, x2};
         double ycoords[3] = {y0, y1, y2};
         xt = dotprod3(xcoords, shpCoords);
         yt = dotprod3(ycoords, shpCoords);

      }

      // Check if the point is inside the element
      // XXX printf("Shape coordinates: %f %f %f %f\n", shpCoords[0], shpCoords[1], shpCoords[2], shpCoords[3]);
      return min(min(shpCoords[0], shpCoords[1]), shpCoords[2]) > -small;      
   }
   else {
      // XXX Not implemented yet
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
   @param uuint, vvint, wwint: velocities interpolated at the final position (xt,yt,zt)

   This is done in four steps:

   1. Find the shape coordinates for the final position in the element.
   2. Interpolate the z level of the top and bottom levels at xt, yt. 
   3. Perform a *binary search* to find the vertical level containing zt. It will have to interpolate the 
      z levels at xt, yt at each step.   
   4. Interpolate the velocities at the final position (xt,yt,zt). This is done by interpolating the velocities
      horizontally and vertically.
**/
__device__ bool Element::interpolateVelocities(double xt, double yt, double zt, 
                                               double& uuint, double& vvint, double& wwint) {
   // Interpolate the velocities at the point (xt,yt,zt). 
   // This function works the same for triangular and quadrilateral elements.
   double shapeCoords[4];
   printf("### interpolateVelocities not implemented yet\n");
   // Error("interpolateVelocities Not implemented yet");
   return false;
}
#ifdef BACKTRACKING   

   auto n0 = nodes[0];
   auto n1 = nodes[1];
   auto n2 = nodes[2];
   auto n3 = nodes[3];
   
   // First step is to compute the shape coordinates of (xt,yt)
   if (!computeShapeCoords(xt, yt, shapeCoords)) return tec_failed_to_compute_area_coords;


   // Second step is to find the level that contains zt
   // Uses a binary search to find the level
   // Note that the z coordinates need to be computed, they are not stored.
   // Also, note that interpolate only takes into account 3 nodes for triangular elements and 4 nodes for quadrilateral elements.
   int L = kbe;   // I use L(eft) and R(ight) as usual, but here L is the bottom level and R is the top level
   int R = nvrt;
   double zLAtNodes[4] = {n0->compZCoord(L), n1->compZCoord(L), n2->compZCoord(L), n3->compZCoord(L)};
   double zRAtNodes[4] = {n0->compZCoord(R), n1->compZCoord(R), n2->compZCoord(R), n3->compZCoord(R)};
   double zL = interpolate(shapeCoords, zLAtNodes);
   double zR = interpolate(shapeCoords, zRAtNodes);
   while (R - L > 1) {
      int M = (L + R) / 2;

      double zMAtNodes[4] = {nodes[0]->compZCoord(M), nodes[1]->compZCoord(M), nodes[2]->compZCoord(M), nodes[3]->compZCoord(M)};
      double zM = interpolate(shapeCoords, zMAtNodes);
      if (zM < zt) 
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
   // At this point, L and R are the two levels that contain zt, and zL and zR are the z coordinates at the nodes of the element at these levels

   // Final step is the vertical interpolation

   // Interpolate horizontally the velocities at the two levels L, R
   const double uLAtNodes[4] = {nodes[0]->u(L), nodes[1]->u(L), nodes[2]->u(L), nodes[3]->u(L)};
   const double vLAtNodes[4] = {nodes[0]->v(L), nodes[1]->v(L), nodes[2]->v(L), nodes[3]->v(L)};
   const double wLAtNodes[4] = {nodes[0]->w(L), nodes[1]->w(L), nodes[2]->w(L), nodes[3]->w(L)};
   const double uLint = interpolate(shapeCoords, uLAtNodes);
   const double vLint = interpolate(shapeCoords, vLAtNodes);
   const double wLint = interpolate(shapeCoords, wLAtNodes);

   const double uRAtNodes[4] = {nodes[0]->u(R), nodes[1]->u(R), nodes[2]->u(R), nodes[3]->u(R)};
   const double vRAtNodes[4] = {nodes[0]->v(R), nodes[1]->v(R), nodes[2]->v(R), nodes[3]->v(R)};
   const double wRAtNodes[4] = {nodes[0]->w(R), nodes[1]->w(R), nodes[2]->w(R), nodes[3]->w(R)};
   const double uRint = interpolate(shapeCoords, uRAtNodes);
   const double vRint = interpolate(shapeCoords, vRAtNodes);
   const double wRint = interpolate(shapeCoords, wRAtNodes);

   // Interpolate vertically the velocities at the level zt
   const double ratio = (zt - zL)/(zR - zL);
   uuint = (1-ratio)*uLint + ratio*uRint;
   vvint = (1-ratio)*vLint + ratio*vRint;
   wwint = (1-ratio)*wLint + ratio*wRint;

}
#endif

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
//       const double x4 = nodes[nqx(i+1)]->x;
//       const double y4 = nodes[nqx(i+1)]->y;

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
              const double vnor1 = side->u2[k  ]*side->nx + side->v2[k  ]*side->ny;
              const double vnor2 = side->u2[k-1]*side->nx + side->v2[k-1]*side->ny;
              
              // Calculate flux (normal * area in local x-direction)
              const double height = side->zc[k] - side->zc[k-1]; // Vertical distance between the two side centers
              const double vmean = 0.5*(vnor1 + vnor2);
              const double length = side->length;
              
              side->flux_adv_hface[k] = height*length*vmean;
          }
      }
   }
}



__device__ void Mesh::populateMesh(int *i34, int *elnode, int *elside,
                double *node_coords, double *wind_rotate_angle, 
                int *isdel, int *isidenode, int *isbnd,
                int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                int *localNodeIndices,  int *localElemIndices,  int *localSideIndices,
                int *localNodePart,     int *localElemPart,     int *localSidePart) {

   printf("Populating mesh for partition %d\n", partNum);


   printf("Processing %d nodes\n", np);      
   for (int i = 0; i < np; i++)
   {
      Node &node = localNodeData[i];

      const int idxStart = nodeNeighNodesStart[i];
      const int idxEnd   = nodeNeighNodesStart[i+1];
      node.numNodes = idxEnd - idxStart;
      node.nodes = (Node **)(nodeNeighNodes + idxStart*sizeof(Node *));

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
      node.z = node_coords[i * 3 + 2];
   }


   printf("Processing %d sides\n", ns);
   for (int i = 0; i < ns; i++)
   {
      Side &side = localSideData[i];

      // finding the pointer to the nodes that form the side
      const int idx_n1 = isidenode[globalSideIndices[i]*2];
      const int idx_n2 = isidenode[globalSideIndices[i]*2 + 1];

      const int locidx_n1 = localNodeIndices[idx_n1];
      const int locidx_n2 = localNodeIndices[idx_n2];
      const int part_n1 =   localNodePart[idx_n1];
      const int part_n2 =   localNodePart[idx_n2];

      side.n1 = (Node *)&(mesh->nodes[part_n1][locidx_n1]);
      side.n2 = (Node *)&(mesh->nodes[part_n2][locidx_n2]);

      // finding the pointer to the elements that share the side
      const int idx_e1 = isdel[globalSideIndices[i]*2];
      const int idx_e2 = isdel[globalSideIndices[i]*2 + 1];

      const int locidx_e1 = localElemIndices[idx_e1];
      const int locidx_e2 = localElemIndices[idx_e2];
      const int part_e1 =  localElemPart[idx_e1];
      const int part_e2 =  localElemPart[idx_e2];

      side.elems[0] = (Element *)&(elems[part_e1][locidx_e1]);
      side.elems[1] = (Element *)&(elems[part_e2][locidx_e2]);

      // Compute the center of the side
      side.xc = 0.5*(side.n1->x + side.n2->x);
      side.yc = 0.5*(side.n1->y + side.n2->y);

      // Compute the z level at the center of the side
      for (int k = 0; k < nvrt; k++) side.zc[k] = 0.5*(side.n1->znl[k] + side.n2->znl[k]);

      // Compute the length of the side
      side.length = sqrt((side.n2->x - side.n1->x)*(side.n2->x - side.n1->x) + 
                         (side.n2->y - side.n1->y)*(side.n2->y - side.n1->y));

   }

   printf("Processing %d elements\n", ne);
   for (int i = 0; i < ne; i++)
   {
      Element &elem = localElemData[i];

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
                           Node *nodes, Element *elems, Side *sides,
                           Node    **partNodeNeighNodesMem, Element **partNodeNeighElemsMem) {

   if (mesh != nullptr) 
   {
      printf("ERROR: Mesh is already initialized\n");
      return;
   };

   mesh = new Mesh(partNum, npgb, negb, nsgb, nvrt, np, ne, ns, max_np, max_ne, max_ns,
                   nodesPtr, elemsPtr, sidesPtr, nodes, elems, sides,
                   partNodeNeighNodesMem, partNodeNeighElemsMem);
   
}

extern "C"
  __global__ void populateMesh(int *i34, int *elnode, int *elside,
                   double *node_coords, double *wind_rotate_angle, 
                   int *isdel, int *isidenode, int *isbnd,
                   int *nodeNeighNodesStart, int *nodeNeighNodesIndices,
                   int *nodeNeighElemsStart, int *nodeNeighElemsIndices,
                   int *globalNodeIndices, int *globalElemIndices, int *globalSideIndices,
                   int *localNodeIndices, int *localElemIndices, int *localSideIndices,
                   int *localNodePart,    int *localElemPart,    int *localSidePart) {

    mesh->populateMesh(i34, elnode, elside, node_coords, wind_rotate_angle, 
                       isdel, isidenode, isbnd,
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