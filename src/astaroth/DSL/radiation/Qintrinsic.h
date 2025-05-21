input int idir
input int3 dir, stop
input real3 unit_vec

Kernel Qintrinsic(int idir, int3 dir, int3 stop, real3 unit_vec){

      int3 boundary = get_boundary(get_normal())   // grid indices of point on starting plane
      int3 sign=sgn(dir)
      stop *= sign

      l,m,n = boundary.x, boundary.y, boundary.z

      idx     = DEVICE_VTXBUF_IDX(l-dir.x,m-dir.y,n-dir.z)
      idx_for = DEVICE_VTXBUF_IDX(l,m,n)
      kapparho_mid = KAPPARHO[idx]
      kapparho_for = KAPPARHO[idx_for]

      srad_mid = SRAD[idx]
      srad_for = SRAD[idx_for]

      dlength_mid = dlength[idir][n-dir.z]
      dlength_for = dlength[idir][n]

      qrad_mid = 0.

      uu_mid = UU[idx]
      uu_for = UU[idx_for]
 
      while (l*sign.x < stop.x && m*sign.y < stop.y && n*sign.z < stop.z)
      {
        idx_back = idx
        idx      = idx_for
        idx_for  = DEVICE_VTXBUF_IDX(l+dir.x,m+dir.y,n+dir.z)

        kapparho_back = kapparho_mid
        kapparho_mid  = kapparho_for
        kapparho_for  = KAPPARHO[idx_for]

        srad_back = srad_mid
        srad_mid  = srad_for
        srad_for  = SRAD[idx_for]

        dlength_back = dlength_mid
        dlength_mid  = dlength_for
        dlength_for  = dlength[idir][n+dir.z]

        dtau_m=sqrt(kapparho_back*kapparho_mid)*0.5*(dlength_back+dlength_mid)
        dtau_p=sqrt(kapparho_mid *kapparho_for)*0.5*(dlength_mid +dlength_for)

        dtau_m = max(dtau_m,epsi)
        dtau_p = max(dtau_p,epsi)

        dSdtau_m = (srad_mid-srad_back)/dtau_m
        dSdtau_p = (srad_for-srad_mid)/dtau_p

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

        tau_mid  = tau_back + dtau_m
        TAU[idx] = tau_mid

        qrad_back = qrad_mid
        qrad_mid  = qrad_back*emdtau-Srad1st*emdtau1-Srad2nd*emdtau2

        if (ldoppler_rad) {

          dQdtau_m = (qrad_mid      - qrad_back)/dtau_m
          dQdtau_p = (QRAD[idx_for] - qrad_mid )/dtau_p   //???
          QRAD1st  = (dQdtau_p*dtau_m+dQdtau_m*dtau_p)/(dtau_m+dtau_p)

          uu_back = uu_mid
          uu_mid  = uu_for
          uu_for  = UU[idx_for]

          u_dot_n = dot(uu_mid,unit_vec)

          if (unit_vec.z)!=0.) {
            hemisign=sgn(unit_vec.z)
          } else {
            if (unit_vec.y)!=0.) {
              hemisign=sign(unit_vec.y)
            } else {
              hemisign=sign(unit_vec.x)
            }
          }

          u_dot_n_p= dot(uu_for,unit_vec)
          u_dot_n_m= dot(uu_back,unit_vec)

          dudtau_m=(u_dot_n-u_dot_n_m)/dtau_m
          dudtau_p=(u_dot_n_p-u_dot_n)/dtau_p

          u_dot_n_1st=(dudtau_p*dtau_m+dudtau_m*dtau_p)/(dtau_m+dtau_p)

          if (ldoppler_rad_includeQ) {

            qrad_mid +=    emdtau1* u_dot_n*(4.*srad_mid+hemisign*qrad_mid)*Q2fact
                       +4.*emdtau2*(u_dot_n*Srad1st+u_dot_n_1st*srad_mid)*Qfact
                       +   emdtau2*(u_dot_n*QRAD1st+u_dot_n_1st*qrad_mid)*Qderfact

          } else {

            qrad_mid +=    emdtau1*u_dot_n*4.*strad_mid
                       +4.*emdtau2*(u_dot_n*Srad1st+u_dot_n_1st*srad_mid)
          }
        }   //  if (ldoppler_rad)

        QRAD[idx] = qrad_mid

        tau_back = tau_mid
        qrad_back = qrad_mid

        l += dir.x
        m += dir.y
        n += dir.z
      }     //  while ...
}
