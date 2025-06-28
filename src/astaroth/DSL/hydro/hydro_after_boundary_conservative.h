Kernel hydro_after_boundary_conservative(real AC_t__mod__cdata){
  real cs201
  real cs2011
  real delx
  real rho
  real rho1
  real press
  real rho_gam21
  real rho_gam20
  real lorentz_gamma2
  real ss2
  real hydro_energy
  real hydro_energy1
  real rat
  real rat0
  real va2_pseudo
  real dely
  real delz
  int iter_relb
  int i
  int j
  if(AC_lconservative__mod__hydro)
  {
  	if (AC_lrelativistic_eos__mod__density) {
  	  cs201=1.+AC_cs20__mod__equationofstate
  	}
  	cs2011=1./cs201
  	if (AC_ldensity__mod__cparam) {
  	  if (AC_lmagnetic__mod__cparam) {
  	    if (AC_b_ext2__mod__magnetic!=0.) {
  	      hydro_energy=value(Field(AC_irho__mod__cdata-1))-0.5*AC_b_ext2__mod__magnetic
  	    }
  	    else {
  	      hydro_energy=value(Field(AC_irho__mod__cdata-1))
  	    }
  	  }
  	  else {
  	    hydro_energy=value(Field(AC_irho__mod__cdata-1))
  	  }
  	  if (AC_lhiggsless__mod__hydro) {
  	    if (AC_lhiggsless_old__mod__hydro) {
  	      for jhless in 1:AC_nhless__mod__hydro+1 {
  	        delx=2.*atan(tan(0.5*(AC_x__mod__cdata[vertexIdx.x]   -AC_xhless__mod__hydro[jhless-1])))
  	        dely=2.*atan(tan(0.5*(AC_y__mod__cdata[AC_m__mod__cdata-1]-AC_yhless__mod__hydro[jhless-1])))
  	        delz=2.*atan(tan(0.5*(AC_z__mod__cdata[AC_n__mod__cdata-1]-AC_zhless__mod__hydro[jhless-1])))
  	        if(sqrt((delx*delx)+(dely*dely)+(delz*delz))  < AC_vwall__mod__hydro*(max(AC_t__mod__cdata-AC_thless__mod__hydro[jhless-1],0.0))) {
  	          DF_HLESS=0.
  	        }
  	        hydro_energy=hydro_energy-value(Field(AC_ihless__mod__hydro-1))
  	      }
  	    }
  	    else {
  	      if(AC_t__mod__cdata < value(Field(AC_ihless__mod__hydro-1))) {
  	        hydro_energy=hydro_energy-AC_alpha_hless__mod__hydro/(1.+AC_alpha_hless__mod__hydro)
  	      }
  	    }
  	  }
  	  hydro_energy1=1./hydro_energy
  	}
  	else {
  	  hydro_energy=1.
  	  hydro_energy1=1.
  	}
  	if (AC_lrelativistic__mod__hydro || AC_llorentz_as_aux__mod__hydro) {
  	  ss=value(F_UVEC)
  	  ss2=(ss.x*ss.x)+(ss.y*ss.y)+(ss.z*ss.z)
  	  rat0=ss2*(hydro_energy1*hydro_energy1)
  	  if (AC_lmagnetic__mod__cparam) {
  	    va2_pseudo=AC_b_ext2__mod__magnetic*cs2011*hydro_energy1
  	    rat=rat0/((1.+va2_pseudo)*(1.+va2_pseudo))
  	  }
  	  else {
  	    rat=rat0
  	  }
  	  lorentz_gamma2=1./(1.-rat)
  	  if (AC_lrelativistic_eos__mod__density) {
  	    lorentz_gamma2=lorentz_gamma2*(0.5-rat*AC_cs20__mod__equationofstate*cs2011 +  sqrt(0.25-rat*AC_cs20__mod__equationofstate*(cs2011*cs2011)))
  	  }
  	  if (AC_lmagnetic__mod__cparam) {
  	    for iter_relb in 1:AC_niter_relb__mod__hydro+1 {
  	      if (AC_lrelativistic__mod__hydro) {
  	        rho1=(cs201*lorentz_gamma2-AC_cs20__mod__equationofstate)/hydro_energy
  	        rho_gam20=cs2011*rho1/lorentz_gamma2
  	        va2_pseudo=AC_b_ext2__mod__magnetic*rho_gam20
  	        rat=rat0/((1.+va2_pseudo)*(1.+va2_pseudo))
  	        lorentz_gamma2=(0.5-rat*AC_cs20__mod__equationofstate*cs2011+sqrt(0.25-rat*AC_cs20__mod__equationofstate*(cs2011*cs2011)))/(1.-rat)
  	      }
  	      else {
  	        if (AC_llorentz_as_aux__mod__hydro) {
  	          lorentz_gamma2=1./(1.-rat)
  	        }
  	      }
  	    }
  	    rho=hydro_energy/(cs201*lorentz_gamma2-AC_cs20__mod__equationofstate)
  	    rho_gam21=1./(cs201*rho*lorentz_gamma2+AC_b_ext2__mod__magnetic)
  	  }
  	  else {
  	    rho=hydro_energy/(cs201*lorentz_gamma2-AC_cs20__mod__equationofstate)
  	    rho_gam21=1./(cs201*rho*lorentz_gamma2)
  	  }
  	}
  	else {
  	  if (AC_lrelativistic_eos__mod__density) {
  	    rho=cs201*hydro_energy
  	  }
  	  else {
  	    rho=hydro_energy
  	  }
  	  rho_gam21=1./(cs201*rho)
  	}
  	if (AC_ilorentz__mod__cdata != 0) {
  	  DF_LORENTZ=lorentz_gamma2
  	}
  	if (AC_lhiggsless__mod__hydro) {
  	  if (AC_lhiggsless_old__mod__hydro) {
  	    press=rho*AC_cs20__mod__equationofstate-value(Field(AC_ihless__mod__hydro-1))
  	  }
  	  else {
  	    press=rho*AC_cs20__mod__equationofstate
  	    if(AC_t__mod__cdata < value(Field(AC_ihless__mod__hydro-1))) {
  	      press=press-AC_alpha_hless__mod__hydro/(1.+AC_alpha_hless__mod__hydro)
  	    }
  	  }
  	}
  	else {
  	  press=rho*AC_cs20__mod__equationofstate
  	}
  	DF_TIJ_0=rho_gam21*(F_UX*F_UX)+press
  	DF_TIJ_1=rho_gam21*(F_UY*F_UY)+press
  	DF_TIJ_2=rho_gam21*(F_UZ*F_UZ)+press
  	DF_TIJ_3=rho_gam21*F_UX*F_UY
  	DF_TIJ_4=rho_gam21*F_UY*F_UZ
  	DF_TIJ_5=rho_gam21*F_UZ*F_UX
  	if (AC_lvv_as_aux__mod__hydro  ||  AC_lvv_as_comaux__mod__hydro) {
  	  DF_VX=rho_gam21*F_UX
  	  DF_VY=rho_gam21*F_UY
  	  DF_VZ=rho_gam21*F_UZ
  	}
  	write(F_TIJ_0,DF_TIJ_0)
  	write(F_TIJ_1,DF_TIJ_1)
  	write(F_TIJ_2,DF_TIJ_2)
  	write(F_TIJ_3,DF_TIJ_3)
  	write(F_TIJ_4,DF_TIJ_4)
  	write(F_TIJ_5,DF_TIJ_5)
  }
}
