//TP: PC boundconds translated to DSL
//In mhd_modular/mhdsolver.ac see an example how to use them in BoundConds
//Unlike before now the syntax is 
//boundcond(BOUNDARY,params...)
//
//i.e. you do not have to give the field as input since the compiler can figure out which fields are written to and read from at least for these bcs (and it is not that crucial to get them 100% correct for the current use case since the interdependency makes the RHS calc still pretty much depend on all of the bcs and vice versa)
//
enum AC_TOP_BOT
{
        AC_top,
        AC_bot
}

//
//
const real density_scale_cgs=1.2435e21;
//TP: note in Fortran this is (-nghost:nghost) so in C this will then be  of length 3*2+1
//Also one has to index into it index+NGHOST_VAL
//
real AC_dx2_bound[7]
real AC_dy2_bound[7]
real AC_dz2_bound[7]

//TP: either make cparam_c.h DSL compatible or load these in gpu_astaroth.cc
//LOADABLE VARIABLES START
int AC_l1
int AC_l2

int AC_n1
int AC_n2

int AC_m1
int AC_m2
//TP: these two should come from cdata.f90 but does not for some reason
real AC_density_scale_factor
int AC_pretend_lnTT
//END


bc_steady_z(topbot,VtxBuffer j)
{
  int i;
  if(topbot == AC_bot) {
    if (j[vertexIdx.x][vertexIdx.y][AC_n1-1] <= 0.0) {
      for i in 1:NGHOST_VAL+1 {
        j[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=j[vertexIdx.x][vertexIdx.y][AC_n1-1];
      }
    }
    else {
      if (j[vertexIdx.x][vertexIdx.y][AC_n1-1] > j[vertexIdx.x][vertexIdx.y][AC_n1+1-1]) {
        j[vertexIdx.x][vertexIdx.y][AC_n1-1-1]=0.5*(j[vertexIdx.x][vertexIdx.y][AC_n1-1]    +j[vertexIdx.x][vertexIdx.y][1+AC_n1-1]);
      }
      else {
        j[vertexIdx.x][vertexIdx.y][AC_n1-1-1]=2.0* j[vertexIdx.x][vertexIdx.y][AC_n1-1]    -j[vertexIdx.x][vertexIdx.y][1+AC_n1-1];
      }
      for i in 2:NGHOST_VAL+1 {
        j[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=2.0* j[vertexIdx.x][vertexIdx.y][AC_n1-i+1-1]-j[vertexIdx.x][vertexIdx.y][AC_n1-i+2-1];
      }
    }
  }
  else if(topbot == AC_top) {
    if (j[vertexIdx.x][vertexIdx.y][AC_n2-1] >= 0.0) {
      for i in 1:NGHOST_VAL+1 {
        j[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=j[vertexIdx.x][vertexIdx.y][AC_n2-1];
      }
    }
    else {
      if (j[vertexIdx.x][vertexIdx.y][AC_n2-1] < j[vertexIdx.x][vertexIdx.y][AC_n2-1-1]) {
        j[vertexIdx.x][vertexIdx.y][1+AC_n2-1]=0.5*(j[vertexIdx.x][vertexIdx.y][AC_n2-1]    +j[vertexIdx.x][vertexIdx.y][AC_n2-1-1]);
      }
      else {
        j[vertexIdx.x][vertexIdx.y][1+AC_n2-1]=2.0* j[vertexIdx.x][vertexIdx.y][AC_n2-1]    -j[vertexIdx.x][vertexIdx.y][AC_n2-1-1];
      }
      for i in 2:NGHOST_VAL+1 {
        j[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=2.0* j[vertexIdx.x][vertexIdx.y][AC_n2+i-1-1]-j[vertexIdx.x][vertexIdx.y][AC_n2+i-2-1];
      }
    }
  }
  else {
  }
}



bc_ss_flux(topbot)
{
  real tmp_xy;
  real cs2_xy;
  real rho_xy;
  int i;
  if(topbot == AC_bot) {
    if (AC_pretend_lnTT) {
      tmp_xy=-FbotKbot/exp(SS[vertexIdx.x][vertexIdx.y][AC_n1-1]);
      for i in 1:NGHOST_VAL+1 {
        SS[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=SS[vertexIdx.x][vertexIdx.y][AC_n1+i-1]-AC_dz2_bound[-i+NGHOST_VAL]*tmp_xy;
      }
    }
    else {
      if (AC_ldensity_nolog) {
          rho_xy=LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-1];
      }
      else {
        rho_xy=exp(LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-1]);
      }
      cs2_xy = SS[vertexIdx.x][vertexIdx.y][AC_n1-1];
      if (AC_ldensity_nolog) {
        cs2_xy=AC_cs20*exp(AC_gamma_m1*(log(rho_xy)-AC_lnrho0)+AC_cv1*cs2_xy);
      }
      else {
        cs2_xy=AC_cs20*exp(AC_gamma_m1*(LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-1]-AC_lnrho0)+AC_cv1*cs2_xy);
      }
      if (AC_lheatc_chiconst) {
        tmp_xy=Fbot/(rho_xy*AC_chi*cs2_xy);
      }
      else if (AC_lheatc_kramers) {
        tmp_xy=Fbot*pow(rho_xy,(2*nkramers))*pow((AC_cp*AC_gamma_m1),(6.5*nkramers))  /(hcond0_kramers*pow(cs2_xy,(6.5*nkramers+1.)));
      }
      else {
        tmp_xy=FbotKbot/cs2_xy;
      }
      for i in 1:NGHOST_VAL+1 {
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][AC_n1+i-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-i-1];
        if (AC_ldensity_nolog) {
            rho_xy = rho_xy/LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-1];
        }
        SS[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=SS[vertexIdx.x][vertexIdx.y][AC_n1+i-1]+AC_cp*(AC_cp-AC_cv)*(rho_xy+AC_dz2_bound[-i+NGHOST_VAL]*tmp_xy);
      }
    }
  }
  else if(topbot == AC_top) {
    if (AC_pretend_lnTT) {
      tmp_xy=-FtopKtop/exp(SS[vertexIdx.x][vertexIdx.y][AC_n2-1]);
      for i in 1:NGHOST_VAL+1 {
        SS[vertexIdx.x][vertexIdx.y][AC_n2-i-1]=SS[vertexIdx.x][vertexIdx.y][AC_n2+i-1]-AC_dz2_bound[i+NGHOST_VAL]*tmp_xy;
      }
    }
    else {
      if (AC_ldensity_nolog) {
          rho_xy=LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-1];
      }
      else {
        rho_xy=exp(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-1]);
      }
      cs2_xy = SS[vertexIdx.x][vertexIdx.y][AC_n2-1];
      if (AC_ldensity_nolog) {
        cs2_xy=AC_cs20*exp(AC_gamma_m1*(log(rho_xy)-AC_lnrho0)+AC_cv1*cs2_xy);
      }
      else {
        cs2_xy=AC_cs20*exp(AC_gamma_m1*(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-1]-AC_lnrho0)+AC_cv1*cs2_xy);
      }
      if (AC_lheatc_chiconst) {
        tmp_xy=Ftop/(rho_xy*AC_chi*cs2_xy);
      }
      else if (AC_lheatc_kramers) {
        tmp_xy=Ftop*pow(rho_xy,(2*nkramers))*pow((AC_cp*AC_gamma_m1),(6.5*nkramers))  /(hcond0_kramers*pow(cs2_xy,(6.5*nkramers+1.)));
      }
      else {
        tmp_xy=FtopKtop/cs2_xy;
      }
      for i in 1:NGHOST_VAL+1 {
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2+i-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-i-1];
        if (AC_ldensity_nolog) {
            rho_xy = rho_xy/LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-1];
        }
        SS[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=SS[vertexIdx.x][vertexIdx.y][AC_n2-i-1]+AC_cp*(AC_cp-AC_cv)*(-rho_xy-AC_dz2_bound[i+NGHOST_VAL]*tmp_xy);
      }
    }
  }
  else {
  }
}



bc_ism(topbot,VtxBuffer j)
{
  int k;
  real density_scale1;
  real density_scale;
  if (AC_density_scale_factor==impossible) {
    density_scale=density_scale_cgs/AC_unit_length;
  }
  else {
    density_scale=AC_density_scale_factor;
  }
  density_scale1=1./density_scale;
  if(topbot == AC_bot) {
    for k in 1:NGHOST_VAL+1 {
      if (j==RHO  ||  j==LNRHO) {
        if (AC_ldensity_nolog) {
          j[vertexIdx.x][vertexIdx.y][k-1]=j[vertexIdx.x][vertexIdx.y][AC_n1-1]*exp(-(AC_z[AC_n1]-AC_z[k])*density_scale1);
        }
        else {
          j[vertexIdx.x][vertexIdx.y][k-1]=j[vertexIdx.x][vertexIdx.y][AC_n1-1] - (AC_z[AC_n1]-AC_z[k])*density_scale1;
        }
      }
      else if (j==SS) {
        if (AC_ldensity_nolog) {
          j[vertexIdx.x][vertexIdx.y][AC_n1-k-1]=j[vertexIdx.x][vertexIdx.y][AC_n1-1]+(AC_cp-AC_cv)*(log(RHO[vertexIdx.x][vertexIdx.y][AC_n1-1])-log(RHO[vertexIdx.x][vertexIdx.y][AC_n1-k-1])) +  AC_cv*log((AC_z[AC_n1]-AC_z[AC_n1-k])*density_scale+1.);
        }
        else {
          j[vertexIdx.x][vertexIdx.y][AC_n1-k-1]=j[vertexIdx.x][vertexIdx.y][AC_n1-1]+(AC_cp-AC_cv)*(LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-k-1])+  AC_cv*log((AC_z[AC_n1]-AC_z[AC_n1-k])*density_scale+1.);
        }
      }
      else {
      }
    }
  }
  else if(topbot == AC_top) {
    for k in 1:NGHOST_VAL+1 {
      if (j==RHO  ||  j==LNRHO) {
        if (AC_ldensity_nolog) {
          j[vertexIdx.x][vertexIdx.y][AC_n2+k-1]=j[vertexIdx.x][vertexIdx.y][AC_n2-1]*exp(-(AC_z[AC_n2+k]-AC_z[AC_n2])*density_scale1);
        }
        else {
          j[vertexIdx.x][vertexIdx.y][AC_n2+k-1]=j[vertexIdx.x][vertexIdx.y][AC_n2-1] - (AC_z[AC_n2+k]-AC_z[AC_n2])*density_scale1;
        }
      }
      else if (j==SS) {
        if (AC_ldensity_nolog) {
          j[vertexIdx.x][vertexIdx.y][AC_n2+k-1]=j[vertexIdx.x][vertexIdx.y][AC_n2-1]+(AC_cp-AC_cv)*(log(RHO[vertexIdx.x][vertexIdx.y][AC_n2-1])-log(RHO[vertexIdx.x][vertexIdx.y][AC_n2+k-1]))+  AC_cv*log((AC_z[AC_n2+k]-AC_z[AC_n2])*density_scale+1.);
        }
        else {
          j[vertexIdx.x][vertexIdx.y][AC_n2+k-1]=j[vertexIdx.x][vertexIdx.y][AC_n2-1]+(AC_cp-AC_cv)*(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n2+k-1])+  AC_cv*log((AC_z[AC_n2+k]-AC_z[AC_n2])*density_scale+1.);
        }
      }
      else {
      }
    }
  }
  else {
  }
}

bc_sym_x(sgn,topbot,VtxBuffer j,rel)
{
  int i;
  if(topbot == AC_bot) {
    if (rel) {
      for i in 1:NGHOST_VAL+1 {
        j[AC_l1-i-1][vertexIdx.y][vertexIdx.z]=2*j[AC_l1-1][vertexIdx.y][vertexIdx.z]+sgn*j[AC_l1+i-1][vertexIdx.y][vertexIdx.z];
      }
    }
    else {
      for i in 1:NGHOST_VAL+1 {
        j[AC_l1-i-1][vertexIdx.y][vertexIdx.z]=              sgn*j[AC_l1+i-1][vertexIdx.y][vertexIdx.z];
      }
      if (sgn<0) {
        j[AC_l1-1][vertexIdx.y][vertexIdx.z] = 0.;
      }
    }
  }
  else if(topbot == AC_top) {
    if (rel) {
      for i in 1:NGHOST_VAL+1 {
        j[AC_l2+i-1][vertexIdx.y][vertexIdx.z]=2*j[AC_l2-1][vertexIdx.y][vertexIdx.z]+sgn*j[AC_l2-i-1][vertexIdx.y][vertexIdx.z];
      }
    }
    else {
      for i in 1:NGHOST_VAL+1 {
        j[AC_l2+i-1][vertexIdx.y][vertexIdx.z]=              sgn*j[AC_l2-i-1][vertexIdx.y][vertexIdx.z];
      }
      if (sgn<0) {
        j[AC_l2-1][vertexIdx.y][vertexIdx.z] = 0.;
      }
    }
  }
  else {
  }
}


bc_sym_y(sgn,topbot,VtxBuffer j,rel)
{
  int i;
  if(topbot == AC_bot) {
    if (rel) {
      for i in 1:NGHOST_VAL+1 {
        j[vertexIdx.x][AC_m1-i-1][vertexIdx.z]=2*j[vertexIdx.x][AC_m1-1][vertexIdx.z]+sgn*j[vertexIdx.x][AC_m1+i-1][vertexIdx.z];
      }
    }
    else {
      for i in 1:NGHOST_VAL+1 {
        j[vertexIdx.x][AC_m1-i-1][vertexIdx.z]=              sgn*j[vertexIdx.x][AC_m1+i-1][vertexIdx.z];
      }
      if (sgn<0) {
        j[vertexIdx.x][AC_m1-1][vertexIdx.z] = 0.;
      }
    }
  }
  else if(topbot == AC_top) {
    if (rel) {
      for i in 1:NGHOST_VAL+1 {
        j[vertexIdx.x][AC_m2+i-1][vertexIdx.z]=2*j[vertexIdx.x][AC_m2-1][vertexIdx.z]+sgn*j[vertexIdx.x][AC_m2-i-1][vertexIdx.z];
      }
    }
    else {
      for i in 1:NGHOST_VAL+1 {
        j[vertexIdx.x][AC_m2+i-1][vertexIdx.z]=              sgn*j[vertexIdx.x][AC_m2-i-1][vertexIdx.z];
      }
      if (sgn<0) {
        j[vertexIdx.x][AC_m2-1][vertexIdx.z] = 0.;
      }
    }
  }
  else {
  }
}


bc_sym_z(sgn,topbot,VtxBuffer j,rel)
{
  int i;
  if(topbot == AC_bot) {
    if (rel) {
      for i in 1:NGHOST_VAL+1 {
        j[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=2*j[vertexIdx.x][vertexIdx.y][AC_n1-1]+sgn*j[vertexIdx.x][vertexIdx.y][AC_n1+i-1];
      }
    }
    else {
      for i in 1:NGHOST_VAL+1 {
        j[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=              sgn*j[vertexIdx.x][vertexIdx.y][AC_n1+i-1];
      }
      if (sgn<0) {
        j[vertexIdx.x][vertexIdx.y][AC_n1-1] = 0.;
      }
    }
  }
  else if(topbot == AC_top) {
    if (rel) {
      for i in 1:NGHOST_VAL+1 {
        j[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=j[vertexIdx.x][vertexIdx.y][AC_n2-1]+(j[vertexIdx.x][vertexIdx.y][AC_n2-1]+sgn*j[vertexIdx.x][vertexIdx.y][AC_n2-i-1]);
      }
    }
    else { 
      for i in 1:NGHOST_VAL+1 {
        j[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=              sgn*j[vertexIdx.x][vertexIdx.y][AC_n2-i-1];
      }
      if (sgn<0) {
        j[vertexIdx.x][vertexIdx.y][AC_n2-1] = 0.;
      }
    }
  }
  else {
  }
}   


bc_set_der_x(topbot,VtxBuffer j,val)
{
  int i;
  if(topbot == AC_bot) {
    for i in 1:NGHOST_VAL+1 {
      j[AC_l1-i-1][vertexIdx.y][vertexIdx.z] = j[AC_l1+i-1][vertexIdx.y][vertexIdx.z] - AC_dx2_bound[-i+NGHOST_VAL+1-1]*val;
    }
  }
  else if(topbot == AC_top) {
    for i in 1:NGHOST_VAL+1 {
      j[AC_l2+i-1][vertexIdx.y][vertexIdx.z] = j[AC_l2-i-1][vertexIdx.y][vertexIdx.z] + AC_dx2_bound[i+1+NGHOST_VAL-1]*val;
    }
  }
  else {
  }
}

bc_set_der_y(topbot,VtxBuffer j,val)
{
  int i;
  if(topbot == AC_bot) {
    for i in 1:NGHOST_VAL+1 {
      j[vertexIdx.x][AC_m1-i-1][vertexIdx.z] = j[vertexIdx.x][AC_m1+i-1][vertexIdx.z] - AC_dy2_bound[-i+NGHOST_VAL+1-1]*val;
    }
  }
  else if(topbot == AC_top) {
    for i in 1:NGHOST_VAL+1 {
      j[vertexIdx.x][AC_m2+i-1][vertexIdx.z] = j[vertexIdx.x][AC_m2-i-1][vertexIdx.z] + AC_dy2_bound[i+1+NGHOST_VAL-1]*val;
    }
  }
  else {
  }
}



bc_set_der_z(topbot,VtxBuffer j,val)
{
  int i;
  if(topbot == AC_bot) {
    for i in 1:NGHOST_VAL+1 {
      j[vertexIdx.x][vertexIdx.y][AC_n1-i-1] = j[vertexIdx.x][vertexIdx.y][AC_n1+i-1] - AC_dz2_bound[-i+NGHOST_VAL+1-1]*val;
    }
  }
  else if(topbot == AC_top) {
    for i in 1:NGHOST_VAL+1 {
      j[vertexIdx.x][vertexIdx.y][AC_n2+i-1] = j[vertexIdx.x][vertexIdx.y][AC_n2-i-1] + AC_dz2_bound[i+1+NGHOST_VAL-1]*val;
    }
  }
  else {
  }
}
