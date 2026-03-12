
    Element *elem;

    elem->setFrame();



    void Element::setFrame(void) {

        double lat, lon;
        computeLatLon(xc, yc, zc, lat, lon);

        frame.update(lat, lon);

    }



    # I1
        if (!isQuad()) {

            // To clarify the notation in this section: i stands for the node under consideration, j and k are the other nodes of the triangle.
            // Node i is related to local index id1, nodes j and k are related to local indices id2 and id3.
            
            id2=nxq(1,id1)
            id3=nxq(2,id1)

            // integral of phi_i*phi_j
            const double ipp_ii = (1.0/6.0 )*area;
            const double ipp_ij = (1.0/12.0)*area;
            const double ipp_ik = (1.0/12.0)*area;

            // integral of [(grad phi_i) * (grad phi_j)]
            const double igpgp_ii = (1.0/4.0)/area*((x[id3]-x[id2])*(x[id3]-x[id2]) + (y[id3]-y[id2])*(y[id3]-y[id2]));
            const double igpgp_ij = (1.0/4.0)/area*((x[id3]-x[id2])*(x[id1]-x[id3]) + (y[id3]-y[id2])*(y[id1]-y[id3]));
            const double igpgp_ik = - igpgp_ii - igpgp_ij; // This comes from the fact that if a,b,c are vectors forming a triangle (in order), then a+b+c = 0, hence: -a*c = a*(a+b+c-c) = a*a + a*b.

            // Computes the different components of integral I_1 of each node with respect to node i.
            const double I1_ii = ipp_ii + grav3*thetai*thetai*dt*dt*hhatb*igpgp_ii;
            const double I1_ij = ipp_ij + grav3*thetai*thetai*dt*dt*hhatb*igpgp_ij;
            const double I1_ik = ipp_ik + grav3*thetai*thetai*dt*dt*hhatb*igpgp_ik;


            sparsem(0,i)=sparsem(0,i)+tmp0


            do jj=1,2 !other 2 nodes
              nd=elnode(nxq(jj,id,i34(ie)),ie)
              indx=0
              do m=1,nnp(i)
                if(indnd(m,i)==nd) then
                  indx=m; exit
                endif
              enddo !m
              if(indx==0) call parallel_abort('STEP: failed to find (9)')
              sparsem(indx,i)=sparsem(indx,i)+swild(jj)
            enddo !jj
        }
        else { // Quad

        }
