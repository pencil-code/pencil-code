Field3 RADFLUX
Field  QRADTOT
Field  RADPRESS[9]
// all correspond to f-array slots

// Provide weight*[idir] as weight*; unit_vec[idir] as unit_vec
Kernel Qfinalize(real weight, real weightn, real unit_vec[3]) {

// Standard kernel: one thread per grid point in subdomain

    idx = DEVICE_VTXBUF_IDX[l,m,n]
    QRADTOT[idx] += weight*QRAD*KAPPARHO

//  Calculate radiative flux.

    if (lradflux || lradpress) {tmp=(QRAD+SRAD)*KAPPARHO}

    if (lradflux) {
      RADFLUX[idx].x += weightn*unit_vec[0]*tmp
      RADFLUX[idx].y += weightn*unit_vec[1]*tmp
      RADFLUX[idx].z += weightn*unit_vec[2]*tmp
    }

//  Calculate radiative pressure.

    if (lradpress) {
      for j in range(0,3){
        for i in range(0,j){
          k=ij_table[i][j]-1     // order!!!
          RADPRESS[idx][k] += weightn*unit_vec[i]*unit_vec[j]*tmp/c_light
        }
      }
    }
}
