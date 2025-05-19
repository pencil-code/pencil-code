input int idir
input int3 dir, stop
input real3 unit_vec

Kernel Qintrinsic(int idir, int3 dir, int3 stop, real3 unit_vec){

      int3 boundary = get_boundary(get_normal())
      int3 sign=sgn(dir)
      stop *= sign

      l,m,n = boundary.x, boundary.y, boundary.z

      int3 sign=sgn(dir)
      stop *= sign

      while (l*sign.x < stop.x && m*sign.y < stop.y && n*sign.z < stop.z)
      {
        idx      = DEVICE_VTXBUF_IDX[l,m,n]
        idx_for  = DEVICE_VTXBUF_IDX[l+dir.x,m+dir.y,n+dir.z)
        idx_back = DEVICE_VTXBUF_IDX[l-dir.x,m-dir.y,n-dir.z)

        dtau_m=sqrt(KAPPARHO[idx_back]*KAPPARHO[idx])* &
                    0.5*(dlength[idir][n-dir.z]+dlength[idir][n])
        dtau_p=sqrt(KAPPARHO[idx_for]*KAPPARHO[idx])* &
                    0.5*(dlength[idir][n]+dlength[idir][n+dir.z])

        dtau_m = max(dtau_m,epsi)
        dtau_p = max(dtau_p,epsi)

        dSdtau_m = (SRAD[idx]-SRAD[idx_back])/dtau_m
        dSdtau_p = (SRAD[idx_for]-SRAD[idx])/dtau_p

        Srad1st=(dSdtau_p*dtau_m+dSdtau_m*dtau_p)/(dtau_m+dtau_p)
        Srad2nd=2.*(dSdtau_p-dSdtau_m)/(dtau_m+dtau_p)

        if (dtau_m>dtau_thresh_max) {
          emdtau=0.0
          emdtau1=1.0
          emdtau2=-1.0
        } else {
          if (dtau_m<dtau_thresh_min) {
            emdtau1=dtau_m*(1.-0.5*dtau_m*(1.-dtau_m/3.))
            emdtau=1-emdtau1
            emdtau2=-dtau_m**2*(0.5-dtau_m/3.)
          } else {
            emdtau=exp(-dtau_m)
            emdtau1=1.-emdtau
            emdtau2=emdtau*(1.+dtau_m)-1.
          }
        }

        TAU[idx]=TAU[idx_back]+dtau_m
        QRAD[idx]=QRAD[idx_back]*emdtau-Srad1st*emdtau1-Srad2nd*emdtau2

        if (ldoppler_rad) {

          dQdtau_m=(QRAD[idx]-QRAD[idx_back])/dtau_m
          dQdtau_p=(QRAD[idx_for]-QRAD[idx])/dtau_p
          QRAD1st=(dQdtau_p*dtau_m+dQdtau_m*dtau_p)/(dtau_m+dtau_p)

          u_dot_n= UUX[idx]*unit_vec.x)  // dot(UU[idx],unit_vec)
                  +UUY[idx]*unit_vec.y)
                  +UUZ[idx]*unit_vec.z)

          if (unit_vec.z)!=0.) {
            hemisign=sgn(unit_vec.z)
          } else {
            if (unit_vec.y)!=0.) {
              hemisign=sign(unit_vec.y)
            } else {
              hemisign=sign(unit_vec.x)
            }
          }

          u_dot_n_p= UUX[idx_for]*unit_vec.x) 
                    +UUY[idx_for]*unit_vec.y) 
                    +UUZ[idx_for]*unit_vec.z)  // dot(UU[idx_for],unit_vec)

          u_dot_n_m= UUX[idx_back]*unit_vec.x) 
                    +UUY[idx_back]*unit_vec.y) 
                    +UUZ[idx_back]*unit_vec.z)  // dot(UU[idx_back],unit_vec)

          dudtau_m=(u_dot_n-u_dot_n_m)/dtau_m
          dudtau_p=(u_dot_n_p-u_dot_n)/dtau_p

          u_dot_n_1st=(dudtau_p*dtau_m+dudtau_m*dtau_p)/(dtau_m+dtau_p)

          if (ldoppler_rad_includeQ) {

            QRAD[idx] +=    emdtau1* u_dot_n*(4.*SRAD[idx]+hemisign*QRAD[idx])*Q2fact
                        +4.*emdtau2*(u_dot_n*Srad1st+u_dot_n_1st*SRAD[idx])*Qfact
                        +   emdtau2*(u_dot_n*QRAD1st+u_dot_n_1st*QRAD[idx])*Qderfact

          } else {

            QRAD[idx] +=    emdtau1*u_dot_n*4.*SRAD[idx]
                        +4.*emdtau2*(u_dot_n*Srad1st+u_dot_n_1st*SRAD[idx])
          }
        }   //  if (ldoppler_rad)
        l += dir.x
        m += dir.y
        n += dir.z
      }  //  while ...
}
