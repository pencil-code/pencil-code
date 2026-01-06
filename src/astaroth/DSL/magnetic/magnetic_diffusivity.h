// checked 17.6
etatot=0.
if (lresi_eta_const){ 
   rhs += eta * laplace(AA)
   etatot+=eta
}
if (lresi_hyper3){
   rhs += eta_hyper3 * del6(AA)
}
if (lresi_ydep) {
  diva = divergence(AA)
  m = vertexIdx.y
  rhs += eta_y[m]*laplace(AA)
  if (lspherical_coords) {
    rhs[2] += geta_y[m]*diva/AC_x[vertexIdx.x]
  } else {
    rhs[2] += geta_y[m]*diva
  }
  etatot += eta_y[m]
}

if (step_num == 0 && ldt && lcourant_dt) {
  dline_1, dxyz_2, dxyz_4, dxyz_6 = get_grid_mn()
  reduce_max(real(etatot*dxyz_2),AC_maxdiffeta)
}
