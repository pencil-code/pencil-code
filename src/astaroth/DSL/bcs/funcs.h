//TP: PC boundconds translated to DSL
//In mhd_modular/mhdsolver.ac see an example how to use them in BoundConds
//Unlike before now the syntax is 
//boundcond(BOUNDARY,params...)
//
//i.e. you do not have to give the field as input since the compiler can figure out which fields are written to and read from at least for these bcs (and it is not that crucial to get them 100% correct for the current use case since the interdependency makes the RHS calc still pretty much depend on all of the bcs and vice versa)
//
enum AC_TOP_BOT
{
        AC_bot,
        AC_top
}
#define BOT AC_bot
#define TOP AC_top

present(var)
{
	return true
}
bc_steady_z(boundary, topbot, VtxBuffer field)
{
  suppress_unused_warning(boundary)
  int i;
  if (topbot == AC_bot) {
    if (field[vertexIdx.x][vertexIdx.y][AC_n1-1] <= 0.0) {
      for i in 1:NGHOST+1 {
        field[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=field[vertexIdx.x][vertexIdx.y][AC_n1-1];
      }
    }
    else {
      if (field[vertexIdx.x][vertexIdx.y][AC_n1-1] > field[vertexIdx.x][vertexIdx.y][AC_n1+1-1]) {
        field[vertexIdx.x][vertexIdx.y][AC_n1-1-1]=0.5*(field[vertexIdx.x][vertexIdx.y][AC_n1-1]    +field[vertexIdx.x][vertexIdx.y][1+AC_n1-1]);
      }
      else {
        field[vertexIdx.x][vertexIdx.y][AC_n1-1-1]=2.0* field[vertexIdx.x][vertexIdx.y][AC_n1-1]    -field[vertexIdx.x][vertexIdx.y][1+AC_n1-1];
      }
      for i in 2:NGHOST+1 {
        field[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=2.0* field[vertexIdx.x][vertexIdx.y][AC_n1-i+1-1]-field[vertexIdx.x][vertexIdx.y][AC_n1-i+2-1];
      }
    }
  }
  else if (topbot == AC_top) {
    if (field[vertexIdx.x][vertexIdx.y][AC_n2-1] >= 0.0) {
      for i in 1:NGHOST+1 {
        field[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=field[vertexIdx.x][vertexIdx.y][AC_n2-1];
      }
    }
    else {
      if (field[vertexIdx.x][vertexIdx.y][AC_n2-1] < field[vertexIdx.x][vertexIdx.y][AC_n2-1-1]) {
        field[vertexIdx.x][vertexIdx.y][1+AC_n2-1]=0.5*(field[vertexIdx.x][vertexIdx.y][AC_n2-1]    +field[vertexIdx.x][vertexIdx.y][AC_n2-1-1]);
      }
      else {
        field[vertexIdx.x][vertexIdx.y][1+AC_n2-1]=2.0* field[vertexIdx.x][vertexIdx.y][AC_n2-1]    -field[vertexIdx.x][vertexIdx.y][AC_n2-1-1];
      }
      for i in 2:NGHOST+1 {
        field[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=2.0* field[vertexIdx.x][vertexIdx.y][AC_n2+i-1-1]-field[vertexIdx.x][vertexIdx.y][AC_n2+i-2-1];
      }
    }
  }
  else {
  }
}
#if LENERGY
bc_ss_flux(boundary, topbot)
{
  suppress_unused_warning(boundary)
  real tmp_xy;
  real cs2_xy;
  real rho_xy;
  int i;
  real lnrho
  lnrho = LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-1]
  if (topbot == AC_bot) {
    if (AC_pretend_lnTT) {
      tmp_xy=-FbotKbot/exp(SS[vertexIdx.x][vertexIdx.y][AC_n1-1]);
      for i in 1:NGHOST+1 {
        SS[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=SS[vertexIdx.x][vertexIdx.y][AC_n1+i-1]-AC_dz2_bound[-i+NGHOST]*tmp_xy;
      }
    } 
    else {
      if (AC_ldensity_nolog) {
          rho_xy=lnrho;
      }
      else {
        rho_xy=exp(lnrho);
      }
      cs2_xy = SS[vertexIdx.x][vertexIdx.y][AC_n1-1];
      if (AC_ldensity_nolog) {
        cs2_xy=AC_cs20*exp(AC_gamma_m1*(log(rho_xy)-AC_lnrho0)+AC_cv1*cs2_xy);
      }
      else {
        cs2_xy=AC_cs20*exp(AC_gamma_m1*(lnrho-AC_lnrho0)+AC_cv1*cs2_xy);
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
      for i in 1:NGHOST+1 {
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][AC_n1+i-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-i-1];
        if (AC_ldensity_nolog) {
            rho_xy = rho_xy/lnrho;
        }
        SS[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=SS[vertexIdx.x][vertexIdx.y][AC_n1+i-1]+AC_cp*(AC_cp-AC_cv)*(rho_xy+AC_dz2_bound[-i+NGHOST]*tmp_xy);
      }
    }
  }
  else if (topbot == AC_top) {
    if (AC_pretend_lnTT) {
      tmp_xy=-FtopKtop/exp(SS[vertexIdx.x][vertexIdx.y][AC_n2-1]);
      for i in 1:NGHOST+1 {
        SS[vertexIdx.x][vertexIdx.y][AC_n2-i-1]=SS[vertexIdx.x][vertexIdx.y][AC_n2+i-1]-AC_dz2_bound[i+NGHOST]*tmp_xy;
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
      for i in 1:NGHOST+1 {
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2+i-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-i-1];
        if (AC_ldensity_nolog) {
            rho_xy = rho_xy/LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-1];
        }
        SS[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=SS[vertexIdx.x][vertexIdx.y][AC_n2-i-1]+AC_cp*(AC_cp-AC_cv)*(-rho_xy-AC_dz2_bound[i+NGHOST]*tmp_xy);
      }
    }
  }
  else {
  }
}
#endif

bc_sym_x(AcBoundary boundary, AC_TOP_BOT topbot, VtxBuffer field,int sgn,bool rel)
{
  suppress_unused_warning(boundary)
  int i;
  if (topbot == AC_bot) {
    if (rel) {
      for i in 1:NGHOST+1 {
        field[AC_l1-i-1][vertexIdx.y][vertexIdx.z]=2*field[AC_l1-1][vertexIdx.y][vertexIdx.z]+sgn*field[AC_l1+i-1][vertexIdx.y][vertexIdx.z];
      }
    }
    else {
      for i in 1:NGHOST+1 {
        field[AC_l1-i-1][vertexIdx.y][vertexIdx.z]=              sgn*field[AC_l1+i-1][vertexIdx.y][vertexIdx.z];
      }
      if (sgn<0) {
        field[AC_l1-1][vertexIdx.y][vertexIdx.z] = 0.;
      }
    }
  }
  else if (topbot == AC_top) {
    if (rel) {
      for i in 1:NGHOST+1 {
        field[AC_l2+i-1][vertexIdx.y][vertexIdx.z]=2*field[AC_l2-1][vertexIdx.y][vertexIdx.z]+sgn*field[AC_l2-i-1][vertexIdx.y][vertexIdx.z];
      }
    }
    else {
      for i in 1:NGHOST+1 {
        field[AC_l2+i-1][vertexIdx.y][vertexIdx.z] = sgn*field[AC_l2-i-1][vertexIdx.y][vertexIdx.z];
      }
      if (sgn<0) {
        field[AC_l2-1][vertexIdx.y][vertexIdx.z] = 0.;
      }
    }
  }
  else {
  }
}

bc_sym_y(AcBoundary boundary, AC_TOP_BOT topbot, VtxBuffer field,int sgn,bool rel)
{
  suppress_unused_warning(boundary)
  int i;
  if (topbot == AC_bot) {
    if (rel) {
      for i in 1:NGHOST+1 {
        field[vertexIdx.x][AC_m1-i-1][vertexIdx.z]=2*field[vertexIdx.x][AC_m1-1][vertexIdx.z]+sgn*field[vertexIdx.x][AC_m1+i-1][vertexIdx.z];
      }
    }
    else {
      for i in 1:NGHOST+1 {
        field[vertexIdx.x][AC_m1-i-1][vertexIdx.z]=              sgn*field[vertexIdx.x][AC_m1+i-1][vertexIdx.z];
      }
      if (sgn<0) {
        field[vertexIdx.x][AC_m1-1][vertexIdx.z] = 0.;
      }
    }
  }
  else if (topbot == AC_top) {
    if (rel) {
      for i in 1:NGHOST+1 {
        field[vertexIdx.x][AC_m2+i-1][vertexIdx.z]=2*field[vertexIdx.x][AC_m2-1][vertexIdx.z]+sgn*field[vertexIdx.x][AC_m2-i-1][vertexIdx.z];
      }
    }
    else {
      for i in 1:NGHOST+1 {
        field[vertexIdx.x][AC_m2+i-1][vertexIdx.z]=              sgn*field[vertexIdx.x][AC_m2-i-1][vertexIdx.z];
      }
      if (sgn<0) {
        field[vertexIdx.x][AC_m2-1][vertexIdx.z] = 0.;
      }
    }
  }
  else {
  }
}

bc_sym_z(AcBoundary boundary, AC_TOP_BOT topbot, VtxBuffer field,int sgn,bool rel)
{
  suppress_unused_warning(boundary)
  int i;
  if (topbot == AC_bot) {
    if (rel) {
      for i in 1:NGHOST+1 {
        field[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=2*field[vertexIdx.x][vertexIdx.y][AC_n1-1]+sgn*field[vertexIdx.x][vertexIdx.y][AC_n1+i-1];
      }
    }
    else {
      for i in 1:NGHOST+1 {
        field[vertexIdx.x][vertexIdx.y][AC_n1-i-1] = sgn*field[vertexIdx.x][vertexIdx.y][AC_n1+i-1];
      }
      if (sgn<0) {
        field[vertexIdx.x][vertexIdx.y][AC_n1-1] = 0.;
      }
    }
  }
  else if (topbot == AC_top) {
    if (rel) {
      for i in 1:NGHOST+1 {
        field[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=field[vertexIdx.x][vertexIdx.y][AC_n2-1]+(field[vertexIdx.x][vertexIdx.y][AC_n2-1]+sgn*field[vertexIdx.x][vertexIdx.y][AC_n2-i-1]);
      }
    }
    else { 
      for i in 1:NGHOST+1 {
        field[vertexIdx.x][vertexIdx.y][AC_n2+i-1] = sgn*field[vertexIdx.x][vertexIdx.y][AC_n2-i-1];
      }
      if (sgn<0) {
        field[vertexIdx.x][vertexIdx.y][AC_n2-1] = 0.;
      }
    }
  }
  else {
  }
}   

bc_set_der_x(boundary, topbot, VtxBuffer field,val)
{
  suppress_unused_warning(boundary)
  int i;
  if (topbot == AC_bot) {
    for i in 1:NGHOST+1 {
      field[AC_l1-i-1][vertexIdx.y][vertexIdx.z] = field[AC_l1+i-1][vertexIdx.y][vertexIdx.z] - AC_dx2_bound[-i+NGHOST+1-1]*val;
    }
  }
  else if (topbot == AC_top) {
    for i in 1:NGHOST+1 {
      field[AC_l2+i-1][vertexIdx.y][vertexIdx.z] = field[AC_l2-i-1][vertexIdx.y][vertexIdx.z] + AC_dx2_bound[i+1+NGHOST-1]*val;
    }
  }
  else {
  }
}

bc_set_der_y(boundary, topbot, VtxBuffer field,val)
{
  suppress_unused_warning(boundary)
  int i;
  if (topbot == AC_bot) {
    for i in 1:NGHOST+1 {
      field[vertexIdx.x][AC_m1-i-1][vertexIdx.z] = field[vertexIdx.x][AC_m1+i-1][vertexIdx.z] - AC_dy2_bound[-i+NGHOST+1-1]*val;
    }
  }
  else if (topbot == AC_top) {
    for i in 1:NGHOST+1 {
      field[vertexIdx.x][AC_m2+i-1][vertexIdx.z] = field[vertexIdx.x][AC_m2-i-1][vertexIdx.z] + AC_dy2_bound[i+1+NGHOST-1]*val;
    }
  }
  else {
  }
}

bc_set_der_z(boundary, topbot, VtxBuffer field,val)
{
  suppress_unused_warning(boundary)
  int i;
  if (topbot == AC_bot) {
    for i in 1:NGHOST+1 {
      field[vertexIdx.x][vertexIdx.y][AC_n1-i-1] = field[vertexIdx.x][vertexIdx.y][AC_n1+i-1] - AC_dz2_bound[-i+NGHOST+1-1]*val;
    }
  }
  else if (topbot == AC_top) {
    for i in 1:NGHOST+1 {
      field[vertexIdx.x][vertexIdx.y][AC_n2+i-1] = field[vertexIdx.x][vertexIdx.y][AC_n2-i-1] + AC_dz2_bound[i+1+NGHOST-1]*val;
    }
  }
  else {
  }
}

#if LENERGY
bc_ss_temp_z(AcBoundary boundary, AC_TOP_BOT topbot,bool lone_sided)
{
suppress_unused_warning(boundary)
real tmp
real lnrho_xy
if(topbot == AC_bot) {;
if (lentropy  &&  ! AC_pretend_lntt) {;
tmp = 2*AC_cv*log(AC_cs2bot/AC_cs20);
if (AC_ldensity_nolog) {;
if (AC_lreference_state) {;
}
else {
lnrho_xy=log(LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-1]);
}
}
else {
lnrho_xy=LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-1];
}
SS[vertexIdx.x][vertexIdx.y][AC_n1-1] = 0.5*tmp - (AC_cp-AC_cv)*(lnrho_xy-AC_lnrho0);
loptest_return_value_1 = false;
if (present(lone_sided)) {;
loptest_return_value_1=lone_sided;
}
else if (false) {;
}
if (loptest_return_value_1) {;
print("not implemented set_ghosts_for_onesided_ders");
}
else {
if (AC_ldensity_nolog) {;
for i in 1:AC_NGHOST+1 {
SS[vertexIdx.x][vertexIdx.y][AC_n1-i-1] = -SS[vertexIdx.x][vertexIdx.y][AC_n1+i-1] + tmp  - 2*(AC_cp-AC_cv)*(lnrho_xy-AC_lnrho0);
}
}
else {
for i in 1:AC_NGHOST+1 {
SS[vertexIdx.x][vertexIdx.y][AC_n1-i-1] = -SS[vertexIdx.x][vertexIdx.y][AC_n1+i-1] + tmp  - (AC_cp-AC_cv)*(LNRHO[vertexIdx.x][vertexIdx.y][AC_n1+i-1]+LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-i-1]-2*AC_lnrho0);
}
}
}
}
else if (lentropy  &&  AC_pretend_lntt) {;
SS[vertexIdx.x][vertexIdx.y][AC_n1-1] = log(AC_cs2bot/AC_gamma_m1);
loptest_return_value_2 = false;
if (present(lone_sided)) {;
loptest_return_value_2=lone_sided;
}
else if (false) {;
}
if (loptest_return_value_2) {;
print("not implemented set_ghosts_for_onesided_ders");
}
else {
for i in 1:AC_NGHOST+1 {
SS[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=2*SS[vertexIdx.x][vertexIdx.y][AC_n1-1]-SS[vertexIdx.x][vertexIdx.y][AC_n1+i-1];
}
}
}
else if (ltemperature) {;
if (ltemperature_nolog) {;
TT[vertexIdx.x][vertexIdx.y][AC_n1-1]   = AC_cs2bot/AC_gamma_m1;
}
else {
LNTT[vertexIdx.x][vertexIdx.y][AC_n1-1] = log(AC_cs2bot/AC_gamma_m1);
}
loptest_return_value_3 = false;
if (present(lone_sided)) {;
loptest_return_value_3=lone_sided;
}
else if (false) {;
}
if (loptest_return_value_3) {;
print("not implemented set_ghosts_for_onesided_ders");
}
else {
for i in 1:AC_NGHOST+1 {
LNTT[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=2*LNTT[vertexIdx.x][vertexIdx.y][AC_n1-1]-LNTT[vertexIdx.x][vertexIdx.y][AC_n1+i-1];
}
}
}
}
else if(topbot == AC_top) {;
cs2top_loc=cs2top;
if (AC_lread_oldsnap) {;
cs2top_loc=AC_cs20*exp(AC_gamma*SS[AC_l2-1][AC_m2-1][AC_n2-1]/AC_cp+AC_gamma_m1*(LNRHO[AC_l2-1][AC_m2-1][AC_n2-1]-AC_lnrho0));
}
if (lentropy  &&  ! AC_pretend_lntt) {;
tmp = 2*AC_cv*log(cs2top_loc/AC_cs20);
if (AC_ldensity_nolog) {;
if (AC_lreference_state) {;
}
else {
lnrho_xy=log(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-1]);
}
}
else {
lnrho_xy=LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-1];
}
SS[vertexIdx.x][vertexIdx.y][AC_n2-1] = 0.5*tmp - (AC_cp-AC_cv)*(lnrho_xy-AC_lnrho0);
loptest_return_value_5 = false;
if (present(lone_sided)) {;
loptest_return_value_5=lone_sided;
}
else if (false) {;
}
if (loptest_return_value_5) {;
print("not implemented set_ghosts_for_onesided_ders");
}
else {
if (AC_ldensity_nolog) {;
for i in 1:AC_NGHOST+1 {
SS[vertexIdx.x][vertexIdx.y][AC_n2+i-1] = -SS[vertexIdx.x][vertexIdx.y][AC_n2-i-1] + tmp  - 2*(AC_cp-AC_cv)*(lnrho_xy-AC_lnrho0);
}
}
else {
for i in 1:AC_NGHOST+1 {
SS[vertexIdx.x][vertexIdx.y][AC_n2+i-1] = -SS[vertexIdx.x][vertexIdx.y][AC_n2-i-1] + tmp  - (AC_cp-AC_cv)*(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-i-1]+LNRHO[vertexIdx.x][vertexIdx.y][AC_n2+i-1]-2*AC_lnrho0);
}
}
}
}
else if (lentropy  &&  AC_pretend_lntt) {;
SS[vertexIdx.x][vertexIdx.y][AC_n2-1] = log(cs2top_loc/AC_gamma_m1);
loptest_return_value_6 = false;
if (present(lone_sided)) {;
loptest_return_value_6=lone_sided;
}
else if (false) {;
}
if (loptest_return_value_6) {;
print("not implemented set_ghosts_for_onesided_ders");
}
else {
for i in 1:AC_NGHOST+1 {
SS[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=2*SS[vertexIdx.x][vertexIdx.y][AC_n2-1]-SS[vertexIdx.x][vertexIdx.y][AC_n2-i-1];
}
}
}
else if (ltemperature) {;
if (ltemperature_nolog) {;
TT[vertexIdx.x][vertexIdx.y][AC_n2-1]   = cs2top_loc/AC_gamma_m1;
}
else {
LNTT[vertexIdx.x][vertexIdx.y][AC_n2-1] = log(cs2top_loc/AC_gamma_m1);
}
loptest_return_value_7 = false;
if (present(lone_sided)) {;
loptest_return_value_7=lone_sided;
}
else if (false) {;
}
if (loptest_return_value_7) {;
print("not implemented set_ghosts_for_onesided_ders");
}
else {
for i in 1:AC_NGHOST+1 {
LNTT[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=2*LNTT[vertexIdx.x][vertexIdx.y][AC_n2-1]-LNTT[vertexIdx.x][vertexIdx.y][AC_n2-i-1];
}
}
}
}
else {
}
}
#endif


