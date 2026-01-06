//checked 17.6.
real3 del2u
real3 graddivu
real3 sglnrho
Matrix sij

nutot=0.
if (lvisc_nu_const || lvisc_nu_profy_bound
#if LSHOCK
                   || lvisc_nu_shock
#endif
   ){
   del2u = laplace(UU)
   graddivu = gradient_of_divergence(UU)
   sij = traceless_rateof_strain(UU)
   sglnrho = sij * glnrho
}
if (lvisc_nu_const){
   rhs +=  nu * (del2u + (1./3.)*graddivu + 2.*sglnrho)
         + zeta * graddivu * rho1
   nutot += nu
}
if (lvisc_hyper3_nu_const){
   rhs += nu_hyper3 * ( del6(UU) + u_dot_grad(UU, gij5(UU), glnrho) )   //p%uij5glnrho
}

//  horizontal viscosity profile with two steps
if (lvisc_nu_profy_bound) {

        l = vertexIdx.x; m = vertexIdx.y
        real3 gradnu // = real3(0.,0.,0.)
        gradnu.x=0.; gradnu.z = 0.
        if (lspherical_coords || lcylindrical_coords){
          gradnu.y=gnu_y[m]/AC_x[l]
	}else{
	  if (lsphere_in_a_box) {
            //gradnu.y=gnu_y[m]/p%r_mn
	  }else{
            gradnu.y=gnu_y[m]
	  }
	}
//  Viscous force: nu(y)*(del2u+graddivu/3+2S.glnrho)+2S.gnu.
        rhs += nu_y[m]*(2.*sglnrho + del2u + 1./3.*graddivu) + 2.*sij*gradnu
	nutot += nu_y[m]
}
#if LSHOCK
if (lvisc_nu_shock){
   divu = divergence(UU)
   rhs += nu_shock*(SHOCK*( divu * glnrho + graddivu ) + divu*gradient(SHOCK))
   nutot += nu_shock*SHOCK
}
#endif

if (step_num == 0 && ldt && lcourant_dt) {
  dline_1, dxyz_2, dxyz_4, dxyz_6 = get_grid_mn()
  reduce_max(real(nutot*dxyz_2),AC_maxdiffnu)
}
