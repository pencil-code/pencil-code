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
//TP: old
#if Leos_idealgas_MODULE && Lentropy_MODULE
bc_ss_flux(boundary, topbot, bool lone_sided)
{
  suppress_unused_warning(boundary)
  suppress_unused_warning(lone_sided)
  real tmp_xy;
  real cs2_xy;
  real rho_xy;
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
        tmp_xy=AC_cp__mod__equationofstate*FbotKbot/cs2_xy;
      }
      for i in 1:NGHOST+1 {
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][AC_n1+i-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n1-i-1];
        if (AC_ldensity_nolog) {
            rho_xy = rho_xy/lnrho;
        }
        SS[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=SS[vertexIdx.x][vertexIdx.y][AC_n1+i-1]+(AC_cp-AC_cv)*(rho_xy+AC_dz2_bound[-i+NGHOST]*tmp_xy);
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
        tmp_xy=AC_cp*FtopKtop/cs2_xy;
      }
      for i in 1:NGHOST+1 {
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2+i-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-i-1];
        if (AC_ldensity_nolog) {
            rho_xy = rho_xy/LNRHO[vertexIdx.x][vertexIdx.y][AC_n2-1];
        }
        SS[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=SS[vertexIdx.x][vertexIdx.y][AC_n2-i-1]+(AC_cp-AC_cv)*(-rho_xy-AC_dz2_bound[i+NGHOST]*tmp_xy);
      }
    }
  }
  else {
  }
}
#endif

bc_copy_x(AcBoundary boundary, AC_TOP_BOT topbot,VtxBuffer j)
{
  suppress_unused_warning(boundary)
  if (topbot == AC_bot) {
    j[AC_l1-1-1][vertexIdx.y][vertexIdx.z]=j[AC_l1-1][vertexIdx.y][vertexIdx.z]
    j[AC_l1-2-1][vertexIdx.y][vertexIdx.z]=j[AC_l1-1][vertexIdx.y][vertexIdx.z]
    j[AC_l1-3-1][vertexIdx.y][vertexIdx.z]=j[AC_l1-1][vertexIdx.y][vertexIdx.z]
  }
  else if (topbot == AC_top) {
    j[1+AC_l2-1][vertexIdx.y][vertexIdx.z]=j[AC_l2-1][vertexIdx.y][vertexIdx.z]
    j[2+AC_l2-1][vertexIdx.y][vertexIdx.z]=j[AC_l2-1][vertexIdx.y][vertexIdx.z]
    j[3+AC_l2-1][vertexIdx.y][vertexIdx.z]=j[AC_l2-1][vertexIdx.y][vertexIdx.z]
  }
  else {
  }
}
bc_copy_y(AcBoundary boundary, AC_TOP_BOT topbot,Field j)
{
  suppress_unused_warning(boundary)
  if (topbot == AC_bot) {
    j[vertexIdx.x][m1-1-1][vertexIdx.z]=j[vertexIdx.x][m1-1][vertexIdx.z]
    j[vertexIdx.x][m1-2-1][vertexIdx.z]=j[vertexIdx.x][m1-1][vertexIdx.z]
    j[vertexIdx.x][m1-3-1][vertexIdx.z]=j[vertexIdx.x][m1-1][vertexIdx.z]
  }
  else if (topbot == AC_top) {
    j[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z]=j[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z]
    j[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z]=j[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z]
    j[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z]=j[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z]
  }
  else {
  }
}

bc_copy_z(AcBoundary boundary, AC_TOP_BOT topbot,Field j)
{
  suppress_unused_warning(boundary)
  if (topbot == AC_bot) {
    j[vertexIdx.x][vertexIdx.y][n1-1-1] = j[vertexIdx.x][vertexIdx.y][n1-1]
    j[vertexIdx.x][vertexIdx.y][n1-2-1] = j[vertexIdx.x][vertexIdx.y][n1-1]
    j[vertexIdx.x][vertexIdx.y][n1-3-1] = j[vertexIdx.x][vertexIdx.y][n1-1]
  }
  else if (topbot == AC_top) {
    j[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1] = j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]
    j[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1] = j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]
    j[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1] = j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]
  }
  else {
  }
}

bc_sym_x(AcBoundary boundary, AC_TOP_BOT topbot, VtxBuffer field,int sgn,bool rel)
{
  suppress_unused_warning(boundary)
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

#if Leos_idealgas_MODULE && Lentropy_MODULE
bc_ss_temp_z(AcBoundary boundary, AC_TOP_BOT topbot,bool lone_sided)
{
suppress_unused_warning(boundary)
real tmp
real lnrho_xy
if (topbot == AC_bot) {;
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
ac_set_ghosts_for_onesided_derivs(boundary,SS)
}
else {
if (AC_ldensity_nolog) {;
for i in 1:NGHOST+1 {
SS[vertexIdx.x][vertexIdx.y][AC_n1-i-1] = -SS[vertexIdx.x][vertexIdx.y][AC_n1+i-1] + tmp  - 2*(AC_cp-AC_cv)*(lnrho_xy-AC_lnrho0);
}
}
else {
for i in 1:NGHOST+1 {
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
ac_set_ghosts_for_onesided_derivs(boundary,SS)
}
else {
for i in 1:NGHOST+1 {
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
ac_set_ghosts_for_onesided_derivs(boundary,LNTT)
}
else {
for i in 1:NGHOST+1 {
LNTT[vertexIdx.x][vertexIdx.y][AC_n1-i-1]=2*LNTT[vertexIdx.x][vertexIdx.y][AC_n1-1]-LNTT[vertexIdx.x][vertexIdx.y][AC_n1+i-1];
}
}
}
}
else if (topbot == AC_top) {;
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
ac_set_ghosts_for_onesided_derivs(boundary,SS)
}
else {
if (AC_ldensity_nolog) {;
for i in 1:NGHOST+1 {
SS[vertexIdx.x][vertexIdx.y][AC_n2+i-1] = -SS[vertexIdx.x][vertexIdx.y][AC_n2-i-1] + tmp  - 2*(AC_cp-AC_cv)*(lnrho_xy-AC_lnrho0);
}
}
else {
for i in 1:NGHOST+1 {
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
ac_set_ghosts_for_onesided_derivs(boundary,SS)
}
else {
for i in 1:NGHOST+1 {
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
ac_set_ghosts_for_onesided_derivs(boundary,LNTT)
}
else {
for i in 1:NGHOST+1 {
LNTT[vertexIdx.x][vertexIdx.y][AC_n2+i-1]=2*LNTT[vertexIdx.x][vertexIdx.y][AC_n2-1]-LNTT[vertexIdx.x][vertexIdx.y][AC_n2-i-1];
}
}
}
}
else {
}
}
#endif

#if LVISCOSITY
#else
const bool AC_llambda_effect__mod__viscosity = false
const real AC_llambda_scale_with_nu__mod__viscosity = 0.0
const real AC_lambda_v0b__mod__viscosity = 0.0
const real AC_lambda_v1b__mod__viscosity = 0.0
const real AC_lambda_v0t__mod__viscosity = 0.0
const real AC_lambda_v1t__mod__viscosity = 0.0
const real AC_nu__mod__viscosity= 0.0
const real AC_lambda_h1__mod__viscosity = 0.0
gmem real AC_lh1_rprof__mod__viscosity[1]
#endif

bc_set_sfree_x(AcBoundary boundary, AC_TOP_BOT topbot,Field j)
{
  suppress_unused_warning(boundary)
  int k
  real fac
  real sth
  real lambda_exp
  if (topbot == AC_bot) {
    if (AC_llambda_effect__mod__viscosity) {
      if (AC_lspherical_coords__mod__cdata){
        if (j==AC_iuz__mod__cdata) {
          sth=AC_sinth__mod__cdata[vertexIdx.y]
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            lambda_exp=1.+(AC_lambda_v0b__mod__viscosity+AC_lambda_v1b__mod__viscosity*sth*sth)
          }
          else {
            lambda_exp=1.+(AC_lambda_v0b__mod__viscosity+AC_lambda_v1b__mod__viscosity*sth*sth)/AC_nu__mod__viscosity
          }
          fac=pow((1.-AC_dx2_bound__mod__cdata[-1+NGHOST+1-1]/AC_x__mod__cdata[1+l1-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[l1-1-1][vertexIdx.y][vertexIdx.z] = j[1+l1-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[l1-1-1][vertexIdx.y][vertexIdx.z] = (j[1+l1-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[1+l1-1]*sth)*fac  -AC_omega__mod__cdata*(AC_x__mod__cdata[1+l1-1]-AC_dx2_bound__mod__cdata[-1+NGHOST+1-1])*sth
          }
          fac=pow((1.-AC_dx2_bound__mod__cdata[-2+NGHOST+1-1]/AC_x__mod__cdata[2+l1-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[l1-2-1][vertexIdx.y][vertexIdx.z] = j[2+l1-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[l1-2-1][vertexIdx.y][vertexIdx.z] = (j[2+l1-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[2+l1-1]*sth)*fac  -AC_omega__mod__cdata*(AC_x__mod__cdata[2+l1-1]-AC_dx2_bound__mod__cdata[-2+NGHOST+1-1])*sth
          }
          fac=pow((1.-AC_dx2_bound__mod__cdata[-3+NGHOST+1-1]/AC_x__mod__cdata[3+l1-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[l1-3-1][vertexIdx.y][vertexIdx.z] = j[3+l1-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[l1-3-1][vertexIdx.y][vertexIdx.z] = (j[3+l1-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[3+l1-1]*sth)*fac  -AC_omega__mod__cdata*(AC_x__mod__cdata[3+l1-1]-AC_dx2_bound__mod__cdata[-3+NGHOST+1-1])*sth
          }
        }
      }
      else if (AC_lcylindrical_coords__mod__cdata) {
        if (j==AC_iuy__mod__cdata) {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            lambda_exp=1.+AC_lambda_v0b__mod__viscosity
          }
          else {
            lambda_exp=1.+AC_lambda_v0b__mod__viscosity/AC_nu__mod__viscosity
          }
          fac=pow((1.-AC_dx2_bound__mod__cdata[-1+NGHOST+1-1]/AC_x__mod__cdata[1+l1-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[l1-1-1][vertexIdx.y][vertexIdx.z] = j[1+l1-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[l1-1-1][vertexIdx.y][vertexIdx.z] = (j[1+l1-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[1+l1-1])*fac
          }
          fac=pow((1.-AC_dx2_bound__mod__cdata[-2+NGHOST+1-1]/AC_x__mod__cdata[2+l1-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[l1-2-1][vertexIdx.y][vertexIdx.z] = j[2+l1-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[l1-2-1][vertexIdx.y][vertexIdx.z] = (j[2+l1-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[2+l1-1])*fac
          }
          fac=pow((1.-AC_dx2_bound__mod__cdata[-3+NGHOST+1-1]/AC_x__mod__cdata[3+l1-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[l1-3-1][vertexIdx.y][vertexIdx.z] = j[3+l1-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[l1-3-1][vertexIdx.y][vertexIdx.z] = (j[3+l1-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[3+l1-1])*fac
          }
        }
      }
      else {
      }
    }
    else {
      j[l1-1-1][vertexIdx.y][vertexIdx.z] = j[1+l1-1][vertexIdx.y][vertexIdx.z]*(1.-AC_dx2_bound__mod__cdata[-1+NGHOST+1-1]/AC_x__mod__cdata[1+l1-1])
      j[l1-2-1][vertexIdx.y][vertexIdx.z] = j[2+l1-1][vertexIdx.y][vertexIdx.z]*(1.-AC_dx2_bound__mod__cdata[-2+NGHOST+1-1]/AC_x__mod__cdata[2+l1-1])
      j[l1-3-1][vertexIdx.y][vertexIdx.z] = j[3+l1-1][vertexIdx.y][vertexIdx.z]*(1.-AC_dx2_bound__mod__cdata[-3+NGHOST+1-1]/AC_x__mod__cdata[3+l1-1])
    }
  }
  else if (topbot == AC_top) {
    if (AC_llambda_effect__mod__viscosity) {
      if (AC_lspherical_coords__mod__cdata){
        if (j==AC_iuz__mod__cdata) {
          sth=AC_sinth__mod__cdata[vertexIdx.y]
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            lambda_exp=1.+(AC_lambda_v0t__mod__viscosity+AC_lambda_v1t__mod__viscosity*sth*sth)
          }
          else {
            lambda_exp=1.+(AC_lambda_v0t__mod__viscosity+AC_lambda_v1t__mod__viscosity*sth*sth)/AC_nu__mod__viscosity
          }
          fac=pow((1.+AC_dx2_bound__mod__cdata[2+NGHOST-1]/AC_x__mod__cdata[AC_l2__mod__cdata-1-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = j[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = (j[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[AC_l2__mod__cdata-1-1]*sth)*fac  -AC_omega__mod__cdata*(AC_x__mod__cdata[AC_l2__mod__cdata-1-1]+AC_dx2_bound__mod__cdata[2+NGHOST-1])*sth
          }
          fac=pow((1.+AC_dx2_bound__mod__cdata[3+NGHOST-1]/AC_x__mod__cdata[AC_l2__mod__cdata-2-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = j[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = (j[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[AC_l2__mod__cdata-2-1]*sth)*fac  -AC_omega__mod__cdata*(AC_x__mod__cdata[AC_l2__mod__cdata-2-1]+AC_dx2_bound__mod__cdata[3+NGHOST-1])*sth
          }
          fac=pow((1.+AC_dx2_bound__mod__cdata[4+NGHOST-1]/AC_x__mod__cdata[AC_l2__mod__cdata-3-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = j[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = (j[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[AC_l2__mod__cdata-3-1]*sth)*fac  -AC_omega__mod__cdata*(AC_x__mod__cdata[AC_l2__mod__cdata-3-1]+AC_dx2_bound__mod__cdata[4+NGHOST-1])*sth
          }
        }
      }
      else if (AC_lcylindrical_coords__mod__cdata) {
        if (j==AC_iuy__mod__cdata) {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            lambda_exp=1.+AC_lambda_v0t__mod__viscosity
          }
          else {
            lambda_exp=1.+AC_lambda_v0t__mod__viscosity/AC_nu__mod__viscosity
          }
          fac=pow((1.+AC_dx2_bound__mod__cdata[2+NGHOST-1]/AC_x__mod__cdata[AC_l2__mod__cdata-1-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = j[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = (j[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[AC_l2__mod__cdata-1-1])*fac
          }
          fac=pow((1.+AC_dx2_bound__mod__cdata[3+NGHOST-1]/AC_x__mod__cdata[AC_l2__mod__cdata-2-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = j[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = (j[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[AC_l2__mod__cdata-2-1])*fac
          }
          fac=pow((1.+AC_dx2_bound__mod__cdata[4+NGHOST-1]/AC_x__mod__cdata[AC_l2__mod__cdata-3-1]),lambda_exp)
          if (AC_omega__mod__cdata==0) {
            j[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = j[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]*fac
          }
          else {
            j[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = (j[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]+AC_omega__mod__cdata*AC_x__mod__cdata[AC_l2__mod__cdata-3-1])*fac
          }
        }
      }
      else {
      }
    }
    else {
      j[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]= j[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]*(1.+AC_dx2_bound__mod__cdata[2+NGHOST-1]/AC_x__mod__cdata[AC_l2__mod__cdata-1-1])
      j[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]= j[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]*(1.+AC_dx2_bound__mod__cdata[3+NGHOST-1]/AC_x__mod__cdata[AC_l2__mod__cdata-2-1])
      j[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]= j[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]*(1.+AC_dx2_bound__mod__cdata[4+NGHOST-1]/AC_x__mod__cdata[AC_l2__mod__cdata-3-1])
    }
  }
  else {
  }
}

bc_set_sfree_y(AcBoundary boundary,AC_TOP_BOT topbot,Field j)
{
  suppress_unused_warning(boundary)
  int k
  real cos2thm_k
  real cos2thmpk
  real somega
  real lh1
  if (AC_llambda_effect__mod__viscosity) {
    lh1=AC_lambda_h1__mod__viscosity*AC_lh1_rprof__mod__viscosity[vertexIdx.x]
  }
  if (topbot == AC_bot) {
    if (AC_llambda_effect__mod__viscosity && (j==AC_iuz__mod__cdata)) {
      if (AC_lambda_h1__mod__viscosity!=0.) {
        cos2thm_k=( AC_costh__mod__cdata[m1-1-1]* AC_costh__mod__cdata[m1-1-1])-(AC_sinth__mod__cdata[m1-1-1]*AC_sinth__mod__cdata[m1-1-1])
        cos2thmpk=( AC_costh__mod__cdata[1+m1-1]* AC_costh__mod__cdata[1+m1-1])-(AC_sinth__mod__cdata[1+m1-1]*AC_sinth__mod__cdata[1+m1-1])
        if (AC_omega__mod__cdata==0) {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            j[vertexIdx.x][m1-1-1][vertexIdx.z]= j[vertexIdx.x][1+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.))*AC_sin1th__mod__cdata[1+m1-1])  *(exp(-lh1*cos2thm_k/(4.))*AC_sinth__mod__cdata[m1-1-1])
          }
          else {
            j[vertexIdx.x][m1-1-1][vertexIdx.z]= j[vertexIdx.x][1+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[1+m1-1])  *(exp(-lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[m1-1-1])
          }
        }
        else {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[m1-1-1]*(  exp(cos2thmpk*lh1/(4.)) /exp((cos2thm_k)*lh1/(4.)) -1.)
            j[vertexIdx.x][m1-1-1][vertexIdx.z]= j[vertexIdx.x][1+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.))*AC_sin1th__mod__cdata[1+m1-1])  *(exp(-lh1*cos2thm_k/(4.))*AC_sinth__mod__cdata[m1-1-1])  +somega
          }
          else {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[m1-1-1]*(  exp(cos2thmpk*lh1/(4.*AC_nu__mod__viscosity)) /exp((cos2thm_k)*lh1/(4.*AC_nu__mod__viscosity)) -1.)
            j[vertexIdx.x][m1-1-1][vertexIdx.z]= j[vertexIdx.x][1+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[1+m1-1])  *(exp(-lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[m1-1-1])  +somega
          }
        }
        cos2thm_k=( AC_costh__mod__cdata[m1-2-1]* AC_costh__mod__cdata[m1-2-1])-(AC_sinth__mod__cdata[m1-2-1]*AC_sinth__mod__cdata[m1-2-1])
        cos2thmpk=( AC_costh__mod__cdata[2+m1-1]* AC_costh__mod__cdata[2+m1-1])-(AC_sinth__mod__cdata[2+m1-1]*AC_sinth__mod__cdata[2+m1-1])
        if (AC_omega__mod__cdata==0) {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            j[vertexIdx.x][m1-2-1][vertexIdx.z]= j[vertexIdx.x][2+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.))*AC_sin1th__mod__cdata[2+m1-1])  *(exp(-lh1*cos2thm_k/(4.))*AC_sinth__mod__cdata[m1-2-1])
          }
          else {
            j[vertexIdx.x][m1-2-1][vertexIdx.z]= j[vertexIdx.x][2+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[2+m1-1])  *(exp(-lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[m1-2-1])
          }
        }
        else {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[m1-2-1]*(  exp(cos2thmpk*lh1/(4.)) /exp((cos2thm_k)*lh1/(4.)) -1.)
            j[vertexIdx.x][m1-2-1][vertexIdx.z]= j[vertexIdx.x][2+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.))*AC_sin1th__mod__cdata[2+m1-1])  *(exp(-lh1*cos2thm_k/(4.))*AC_sinth__mod__cdata[m1-2-1])  +somega
          }
          else {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[m1-2-1]*(  exp(cos2thmpk*lh1/(4.*AC_nu__mod__viscosity)) /exp((cos2thm_k)*lh1/(4.*AC_nu__mod__viscosity)) -1.)
            j[vertexIdx.x][m1-2-1][vertexIdx.z]= j[vertexIdx.x][2+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[2+m1-1])  *(exp(-lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[m1-2-1])  +somega
          }
        }
        cos2thm_k=( AC_costh__mod__cdata[m1-3-1]* AC_costh__mod__cdata[m1-3-1])-(AC_sinth__mod__cdata[m1-3-1]*AC_sinth__mod__cdata[m1-3-1])
        cos2thmpk=( AC_costh__mod__cdata[3+m1-1]* AC_costh__mod__cdata[3+m1-1])-(AC_sinth__mod__cdata[3+m1-1]*AC_sinth__mod__cdata[3+m1-1])
        if (AC_omega__mod__cdata==0) {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            j[vertexIdx.x][m1-3-1][vertexIdx.z]= j[vertexIdx.x][3+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.))*AC_sin1th__mod__cdata[3+m1-1])  *(exp(-lh1*cos2thm_k/(4.))*AC_sinth__mod__cdata[m1-3-1])
          }
          else {
            j[vertexIdx.x][m1-3-1][vertexIdx.z]= j[vertexIdx.x][3+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[3+m1-1])  *(exp(-lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[m1-3-1])
          }
        }
        else {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[m1-3-1]*(  exp(cos2thmpk*lh1/(4.)) /exp((cos2thm_k)*lh1/(4.)) -1.)
            j[vertexIdx.x][m1-3-1][vertexIdx.z]= j[vertexIdx.x][3+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.))*AC_sin1th__mod__cdata[3+m1-1])  *(exp(-lh1*cos2thm_k/(4.))*AC_sinth__mod__cdata[m1-3-1])  +somega
          }
          else {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[m1-3-1]*(  exp(cos2thmpk*lh1/(4.*AC_nu__mod__viscosity)) /exp((cos2thm_k)*lh1/(4.*AC_nu__mod__viscosity)) -1.)
            j[vertexIdx.x][m1-3-1][vertexIdx.z]= j[vertexIdx.x][3+m1-1][vertexIdx.z]*  (exp(lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[3+m1-1])  *(exp(-lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[m1-3-1])  +somega
          }
        }
      }
    }
    else {
      j[vertexIdx.x][m1-1-1][vertexIdx.z]= j[vertexIdx.x][1+m1-1][vertexIdx.z]*(sin(AC_y__mod__cdata[1+m1-1]-AC_dy2_bound__mod__cdata[-1+NGHOST+1-1])*AC_sin1th__mod__cdata[1+m1-1])
      j[vertexIdx.x][m1-2-1][vertexIdx.z]= j[vertexIdx.x][2+m1-1][vertexIdx.z]*(sin(AC_y__mod__cdata[2+m1-1]-AC_dy2_bound__mod__cdata[-2+NGHOST+1-1])*AC_sin1th__mod__cdata[2+m1-1])
      j[vertexIdx.x][m1-3-1][vertexIdx.z]= j[vertexIdx.x][3+m1-1][vertexIdx.z]*(sin(AC_y__mod__cdata[3+m1-1]-AC_dy2_bound__mod__cdata[-3+NGHOST+1-1])*AC_sin1th__mod__cdata[3+m1-1])
    }
  }
  else if (topbot == AC_top) {
    if ((AC_llambda_effect__mod__viscosity) && (j==AC_iuz__mod__cdata)) {
      if (AC_lambda_h1__mod__viscosity!=0.) {
        cos2thm_k=( AC_costh__mod__cdata[AC_m2__mod__cdata-1-1]* AC_costh__mod__cdata[AC_m2__mod__cdata-1-1])-(AC_sinth__mod__cdata[AC_m2__mod__cdata-1-1]*AC_sinth__mod__cdata[AC_m2__mod__cdata-1-1])
        cos2thmpk=( AC_costh__mod__cdata[1+AC_m2__mod__cdata-1]* AC_costh__mod__cdata[1+AC_m2__mod__cdata-1])-(AC_sinth__mod__cdata[1+AC_m2__mod__cdata-1]*AC_sinth__mod__cdata[1+AC_m2__mod__cdata-1])
        if (AC_omega__mod__cdata==0){
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            j[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-1-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-1-1])  *(exp(-lh1*cos2thmpk/(4.))*AC_sinth__mod__cdata[1+AC_m2__mod__cdata-1])
          }
          else {
            j[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-1-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-1-1])  *(exp(-lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[1+AC_m2__mod__cdata-1])
          }
        }
        else {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[1+AC_m2__mod__cdata-1]*(  exp(cos2thm_k*lh1/(4.))     /exp(cos2thmpk*lh1/(4.))-1.)
            j[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-1-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-1-1])  *(exp(-lh1*cos2thmpk/(4.))*AC_sinth__mod__cdata[1+AC_m2__mod__cdata-1])  +somega
          }
          else {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[1+AC_m2__mod__cdata-1]*(  exp(cos2thm_k*lh1/(4.*AC_nu__mod__viscosity))     /exp(cos2thmpk*lh1/(4.*AC_nu__mod__viscosity))-1.)
            j[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-1-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-1-1])  *(exp(-lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[1+AC_m2__mod__cdata-1])  +somega
          }
        }
        cos2thm_k=( AC_costh__mod__cdata[AC_m2__mod__cdata-2-1]* AC_costh__mod__cdata[AC_m2__mod__cdata-2-1])-(AC_sinth__mod__cdata[AC_m2__mod__cdata-2-1]*AC_sinth__mod__cdata[AC_m2__mod__cdata-2-1])
        cos2thmpk=( AC_costh__mod__cdata[2+AC_m2__mod__cdata-1]* AC_costh__mod__cdata[2+AC_m2__mod__cdata-1])-(AC_sinth__mod__cdata[2+AC_m2__mod__cdata-1]*AC_sinth__mod__cdata[2+AC_m2__mod__cdata-1])
        if (AC_omega__mod__cdata==0){
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            j[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-2-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-2-1])  *(exp(-lh1*cos2thmpk/(4.))*AC_sinth__mod__cdata[2+AC_m2__mod__cdata-1])
          }
          else {
            j[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-2-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-2-1])  *(exp(-lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[2+AC_m2__mod__cdata-1])
          }
        }
        else {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[2+AC_m2__mod__cdata-1]*(  exp(cos2thm_k*lh1/(4.))     /exp(cos2thmpk*lh1/(4.))-1.)
            j[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-2-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-2-1])  *(exp(-lh1*cos2thmpk/(4.))*AC_sinth__mod__cdata[2+AC_m2__mod__cdata-1])  +somega
          }
          else {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[2+AC_m2__mod__cdata-1]*(  exp(cos2thm_k*lh1/(4.*AC_nu__mod__viscosity))     /exp(cos2thmpk*lh1/(4.*AC_nu__mod__viscosity))-1.)
            j[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-2-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-2-1])  *(exp(-lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[2+AC_m2__mod__cdata-1])  +somega
          }
        }
        cos2thm_k=( AC_costh__mod__cdata[AC_m2__mod__cdata-3-1]* AC_costh__mod__cdata[AC_m2__mod__cdata-3-1])-(AC_sinth__mod__cdata[AC_m2__mod__cdata-3-1]*AC_sinth__mod__cdata[AC_m2__mod__cdata-3-1])
        cos2thmpk=( AC_costh__mod__cdata[3+AC_m2__mod__cdata-1]* AC_costh__mod__cdata[3+AC_m2__mod__cdata-1])-(AC_sinth__mod__cdata[3+AC_m2__mod__cdata-1]*AC_sinth__mod__cdata[3+AC_m2__mod__cdata-1])
        if (AC_omega__mod__cdata==0){
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            j[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-3-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-3-1])  *(exp(-lh1*cos2thmpk/(4.))*AC_sinth__mod__cdata[3+AC_m2__mod__cdata-1])
          }
          else {
            j[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-3-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-3-1])  *(exp(-lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[3+AC_m2__mod__cdata-1])
          }
        }
        else {
          if (AC_llambda_scale_with_nu__mod__viscosity) {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[3+AC_m2__mod__cdata-1]*(  exp(cos2thm_k*lh1/(4.))     /exp(cos2thmpk*lh1/(4.))-1.)
            j[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-3-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-3-1])  *(exp(-lh1*cos2thmpk/(4.))*AC_sinth__mod__cdata[3+AC_m2__mod__cdata-1])  +somega
          }
          else {
            somega=AC_x__mod__cdata[vertexIdx.x]*AC_omega__mod__cdata*AC_sinth__mod__cdata[3+AC_m2__mod__cdata-1]*(  exp(cos2thm_k*lh1/(4.*AC_nu__mod__viscosity))     /exp(cos2thmpk*lh1/(4.*AC_nu__mod__viscosity))-1.)
            j[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-3-1][vertexIdx.z]*  (exp(lh1*cos2thm_k/(4.*AC_nu__mod__viscosity))*AC_sin1th__mod__cdata[AC_m2__mod__cdata-3-1])  *(exp(-lh1*cos2thmpk/(4.*AC_nu__mod__viscosity))*AC_sinth__mod__cdata[3+AC_m2__mod__cdata-1])  +somega
          }
        }
      }
    }
    else {
      j[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-1-1][vertexIdx.z]*(sin(AC_y__mod__cdata[AC_m2__mod__cdata-1-1]+AC_dy2_bound__mod__cdata[2+NGHOST-1])*AC_sin1th__mod__cdata[AC_m2__mod__cdata-1-1])
      j[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-2-1][vertexIdx.z]*(sin(AC_y__mod__cdata[AC_m2__mod__cdata-2-1]+AC_dy2_bound__mod__cdata[3+NGHOST-1])*AC_sin1th__mod__cdata[AC_m2__mod__cdata-2-1])
      j[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-3-1][vertexIdx.z]*(sin(AC_y__mod__cdata[AC_m2__mod__cdata-3-1]+AC_dy2_bound__mod__cdata[4+NGHOST-1])*AC_sin1th__mod__cdata[AC_m2__mod__cdata-3-1])
    }
  }
  else {
  }
}

#if Leos_idealgas_MODULE && Lentropy_MODULE
bc_ss_flux_x(AcBoundary boundary, AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  real work_yz
  real tmp_yz
  int stat
  real fac
  real cp_loc
  real get_cp_return_value_0_2
  real get_cv_return_value_1_2
  real get_cp_return_value_0_3
  real get_cv_return_value_1_3
  int i_5
  int ll_5
  int ia_5
  int ie_5
  bool heatflux_deriv_x_return_value_4_5
  int i_6
  int ll_6
  int ia_6
  int ie_6
  bool heatflux_deriv_x_return_value_4_6
  int i_7
  int ll_7
  int ia_7
  int ie_7
  bool heatflux_deriv_x_return_value_4_7
  int i_8
  int ll_8
  int ia_8
  int ie_8
  bool heatflux_deriv_x_return_value_4_8
  if (AC_lheatc_kramers__mod__energy) {
      get_cp_return_value_0_2=AC_cp__mod__equationofstate
      cp_loc=get_cp_return_value_0_2
  }
  if (AC_lheatc_chiconst__mod__energy) {
      get_cp_return_value_0_3=AC_cp__mod__equationofstate
      cp_loc=get_cp_return_value_0_3
  }
  fac=AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate
  if (topbot == AC_bot) {
    if (AC_pretend_lntt__mod__cdata) {
      SS[l1-1-1][vertexIdx.y][vertexIdx.z]=SS[1+l1-1][vertexIdx.y][vertexIdx.z]+AC_dx2_bound__mod__cdata[-1+NGHOST+1-1]*AC_FbotKbot__mod__energy/exp(SS[l1-1][vertexIdx.y][vertexIdx.z])
      SS[l1-2-1][vertexIdx.y][vertexIdx.z]=SS[2+l1-1][vertexIdx.y][vertexIdx.z]+AC_dx2_bound__mod__cdata[-2+NGHOST+1-1]*AC_FbotKbot__mod__energy/exp(SS[l1-1][vertexIdx.y][vertexIdx.z])
      SS[l1-3-1][vertexIdx.y][vertexIdx.z]=SS[3+l1-1][vertexIdx.y][vertexIdx.z]+AC_dx2_bound__mod__cdata[-3+NGHOST+1-1]*AC_FbotKbot__mod__energy/exp(SS[l1-1][vertexIdx.y][vertexIdx.z])
    }
    else {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lheatc_kramers__mod__energy || AC_lheatc_chiconst__mod__energy) {
          work_yz=RHO[l1-1][vertexIdx.y][vertexIdx.z]
        }
        if (AC_lreference_state__mod__cdata) {
          tmp_yz= AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(log(RHO[l1-1][vertexIdx.y][vertexIdx.z]+AC_reference_state__mod__density[1-1][iref_rho-1])-AC_lnrho0__mod__equationofstate)   +AC_gamma__mod__equationofstate*(SS[l1-1][vertexIdx.y][vertexIdx.z]+AC_reference_state__mod__density[1-1][iref_s-1]))
        }
        else {
          tmp_yz=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(log(RHO[l1-1][vertexIdx.y][vertexIdx.z])-AC_lnrho0__mod__equationofstate)+AC_gamma__mod__equationofstate*SS[l1-1][vertexIdx.y][vertexIdx.z])
        }
      }
      else {
        if (AC_lheatc_kramers__mod__energy || AC_lheatc_chiconst__mod__energy) {
          work_yz=exp(LNRHO[l1-1][vertexIdx.y][vertexIdx.z])
        }
        tmp_yz=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(LNRHO[l1-1][vertexIdx.y][vertexIdx.z]-AC_lnrho0__mod__equationofstate)+AC_gamma__mod__equationofstate*SS[l1-1][vertexIdx.y][vertexIdx.z])
      }
      if (AC_lheatc_kramers__mod__energy) {
        krho1kr_yz = AC_hcond0_kramers__mod__energy*pow(work_yz,(-2*AC_nkramers__mod__energy-1))*pow((tmp_yz/(cp_loc*AC_gamma_m1__mod__equationofstate)),(6.5*AC_nkramers__mod__energy))
        if (AC_chimin_kramers__mod__energy>0) {
          krho1kr_yz = max(krho1kr_yz, AC_chimin_kramers__mod__energy*cp_loc)
        }
        if (AC_chimax_kramers__mod__energy>0) {
          krho1kr_yz = min(krho1kr_yz, AC_chimax_kramers__mod__energy*cp_loc)
        }
        tmp_yz=AC_Fbot__mod__energy/(work_yz*krho1kr_yz*tmp_yz)
      }
      else if (AC_lheatc_chiconst__mod__energy) {
        tmp_yz=AC_Fbot__mod__energy/(work_yz*AC_chi__mod__energy*cp_loc*tmp_yz)
      }
      else {
        tmp_yz=AC_FbotKbot__mod__energy/tmp_yz
      }
      if (AC_lreference_state__mod__cdata) {
        work_yz= 1./(RHO[l1-1][vertexIdx.y][vertexIdx.z]+AC_reference_state__mod__density[1-1][iref_rho-1])
        tmp_yz = tmp_yz + AC_reference_state__mod__density[1-1][iref_gs-1]/fac + AC_reference_state__mod__density[1-1][iref_grho-1]*work_yz
        heatflux_deriv_x_return_value_4_5 = true
        if ( AC_lequidist__mod__cdata.x  ||  ! heatflux_deriv_x_return_value_4_5) {
          if (AC_bot==AC_bot) {
            ll_5=l1
            ia_5=1
            ie_5=NGHOST
          }
          else {
            ll_5=AC_l2__mod__cdata
            ia_5=-NGHOST
            ie_5=-1
          }
          for i_5 in ia_5:ie_5+1 {
            if (AC_ldensity_nolog__mod__cdata) {
              if (present(work_yz)) {
                SS[ll_5-i_5-1][vertexIdx.y][vertexIdx.z]=SS[ll_5+i_5-1][vertexIdx.y][vertexIdx.z]+fac*  ( (RHO[ll_5+i_5-1][vertexIdx.y][vertexIdx.z]-RHO[ll_5-i_5-1][vertexIdx.y][vertexIdx.z])*work_yz + AC_dx2_bound__mod__cdata[-i_5+NGHOST+1-1]*tmp_yz )
              }
              else {
                SS[ll_5-i_5-1][vertexIdx.y][vertexIdx.z]=SS[ll_5+i_5-1][vertexIdx.y][vertexIdx.z]+fac*  (log(RHO[ll_5+i_5-1][vertexIdx.y][vertexIdx.z]/RHO[ll_5-i_5-1][vertexIdx.y][vertexIdx.z]) + AC_dx2_bound__mod__cdata[-i_5+NGHOST+1-1]*tmp_yz)
              }
            }
            else {
              SS[ll_5-i_5-1][vertexIdx.y][vertexIdx.z]=SS[ll_5+i_5-1][vertexIdx.y][vertexIdx.z]+fac*  (LNRHO[ll_5+i_5-1][vertexIdx.y][vertexIdx.z]-LNRHO[ll_5-i_5-1][vertexIdx.y][vertexIdx.z] + AC_dx2_bound__mod__cdata[-i_5+NGHOST+1-1]*tmp_yz)
            }
          }
        }
      }
      else {
        heatflux_deriv_x_return_value_4_6 = true
        if ( AC_lequidist__mod__cdata.x  ||  ! heatflux_deriv_x_return_value_4_6) {
          if (AC_bot==AC_bot) {
            ll_6=l1
            ia_6=1
            ie_6=NGHOST
          }
          else {
            ll_6=AC_l2__mod__cdata
            ia_6=-NGHOST
            ie_6=-1
          }
          for i_6 in ia_6:ie_6+1 {
            if (AC_ldensity_nolog__mod__cdata) {
              if (false) {
              }
              else {
                SS[ll_6-i_6-1][vertexIdx.y][vertexIdx.z]=SS[ll_6+i_6-1][vertexIdx.y][vertexIdx.z]+fac*  (log(RHO[ll_6+i_6-1][vertexIdx.y][vertexIdx.z]/RHO[ll_6-i_6-1][vertexIdx.y][vertexIdx.z]) + AC_dx2_bound__mod__cdata[-i_6+NGHOST+1-1]*tmp_yz)
              }
            }
            else {
              SS[ll_6-i_6-1][vertexIdx.y][vertexIdx.z]=SS[ll_6+i_6-1][vertexIdx.y][vertexIdx.z]+fac*  (LNRHO[ll_6+i_6-1][vertexIdx.y][vertexIdx.z]-LNRHO[ll_6-i_6-1][vertexIdx.y][vertexIdx.z] + AC_dx2_bound__mod__cdata[-i_6+NGHOST+1-1]*tmp_yz)
            }
          }
        }
      }
    }
  }
  else if (topbot == AC_top) {
    if (AC_pretend_lntt__mod__cdata) {
      SS[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]-AC_dx2_bound__mod__cdata[2+NGHOST-1]*AC_FtopKtop__mod__energy/exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
      SS[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]-AC_dx2_bound__mod__cdata[3+NGHOST-1]*AC_FtopKtop__mod__energy/exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
      SS[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]-AC_dx2_bound__mod__cdata[4+NGHOST-1]*AC_FtopKtop__mod__energy/exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
    }
    else {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lheatc_kramers__mod__energy || AC_lheatc_chiconst__mod__energy) {
          work_yz=RHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]
        }
        if (AC_lreference_state__mod__cdata) {
          tmp_yz=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(log(RHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]+AC_reference_state__mod__density[nx-1][iref_rho-1])-AC_lnrho0__mod__equationofstate)  +AC_gamma__mod__equationofstate*(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]+AC_reference_state__mod__density[nx-1][iref_s-1]))
        }
        else {
          tmp_yz=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(log(RHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])-AC_lnrho0__mod__equationofstate)+AC_gamma__mod__equationofstate*SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
        }
      }
      else {
        if (AC_lheatc_kramers__mod__energy || AC_lheatc_chiconst__mod__energy) {
          work_yz=exp(LNRHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
        }
        tmp_yz=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(LNRHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-AC_lnrho0__mod__equationofstate)+AC_gamma__mod__equationofstate*SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
      }
      if (AC_lheatc_kramers__mod__energy) {
        krho1kr_yz = AC_hcond0_kramers__mod__energy*pow(work_yz,(-2*AC_nkramers__mod__energy-1))*pow((tmp_yz/(cp_loc*AC_gamma_m1__mod__equationofstate)),(6.5*AC_nkramers__mod__energy))
        if (AC_chimin_kramers__mod__energy>0) {
          krho1kr_yz = max(krho1kr_yz, AC_chimin_kramers__mod__energy*cp_loc)
        }
        if (AC_chimax_kramers__mod__energy>0) {
          krho1kr_yz = min(krho1kr_yz, AC_chimax_kramers__mod__energy*cp_loc)
        }
        tmp_yz=ftop/(work_yz*krho1kr_yz*tmp_yz)
      }
      else if (AC_lheatc_chiconst__mod__energy) {
        tmp_yz=ftop/(work_yz*AC_chi__mod__energy*cp_loc*tmp_yz)
      }
      else {
        tmp_yz=AC_FtopKtop__mod__energy/tmp_yz
      }
      if (AC_lreference_state__mod__cdata) {
        tmp_yz = tmp_yz + AC_reference_state__mod__density[nx-1][iref_gs-1]
      }
      if (AC_lreference_state__mod__cdata) {
        work_yz= 1./(RHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]+AC_reference_state__mod__density[nx-1][iref_rho-1])
        tmp_yz = tmp_yz + AC_reference_state__mod__density[nx-1][iref_gs-1]/fac + AC_reference_state__mod__density[nx-1][iref_grho-1]*work_yz
        heatflux_deriv_x_return_value_4_7 = true
        if ( AC_lequidist__mod__cdata.x  ||  ! heatflux_deriv_x_return_value_4_7) {
          if (AC_top==AC_bot) {
            ll_7=l1
            ia_7=1
            ie_7=NGHOST
          }
          else {
            ll_7=AC_l2__mod__cdata
            ia_7=-NGHOST
            ie_7=-1
          }
          for i_7 in ia_7:ie_7+1 {
            if (AC_ldensity_nolog__mod__cdata) {
              if (present(work_yz)) {
                SS[ll_7-i_7-1][vertexIdx.y][vertexIdx.z]=SS[ll_7+i_7-1][vertexIdx.y][vertexIdx.z]+fac*  ( (RHO[ll_7+i_7-1][vertexIdx.y][vertexIdx.z]-RHO[ll_7-i_7-1][vertexIdx.y][vertexIdx.z])*work_yz + AC_dx2_bound__mod__cdata[-i_7+NGHOST+1-1]*-tmp_yz )
              }
              else {
                SS[ll_7-i_7-1][vertexIdx.y][vertexIdx.z]=SS[ll_7+i_7-1][vertexIdx.y][vertexIdx.z]+fac*  (log(RHO[ll_7+i_7-1][vertexIdx.y][vertexIdx.z]/RHO[ll_7-i_7-1][vertexIdx.y][vertexIdx.z]) + AC_dx2_bound__mod__cdata[-i_7+NGHOST+1-1]*-tmp_yz)
              }
            }
            else {
              SS[ll_7-i_7-1][vertexIdx.y][vertexIdx.z]=SS[ll_7+i_7-1][vertexIdx.y][vertexIdx.z]+fac*  (LNRHO[ll_7+i_7-1][vertexIdx.y][vertexIdx.z]-LNRHO[ll_7-i_7-1][vertexIdx.y][vertexIdx.z] + AC_dx2_bound__mod__cdata[-i_7+NGHOST+1-1]*-tmp_yz)
            }
          }
        }
      }
      else {
        heatflux_deriv_x_return_value_4_8 = true
        if ( AC_lequidist__mod__cdata.x  ||  ! heatflux_deriv_x_return_value_4_8) {
          if (AC_top==AC_bot) {
            ll_8=l1
            ia_8=1
            ie_8=NGHOST
          }
          else {
            ll_8=AC_l2__mod__cdata
            ia_8=-NGHOST
            ie_8=-1
          }
          for i_8 in ia_8:ie_8+1 {
            if (AC_ldensity_nolog__mod__cdata) {
              if (false) {
              }
              else {
                SS[ll_8-i_8-1][vertexIdx.y][vertexIdx.z]=SS[ll_8+i_8-1][vertexIdx.y][vertexIdx.z]+fac*  (log(RHO[ll_8+i_8-1][vertexIdx.y][vertexIdx.z]/RHO[ll_8-i_8-1][vertexIdx.y][vertexIdx.z]) + AC_dx2_bound__mod__cdata[-i_8+NGHOST+1-1]*-tmp_yz)
              }
            }
            else {
              SS[ll_8-i_8-1][vertexIdx.y][vertexIdx.z]=SS[ll_8+i_8-1][vertexIdx.y][vertexIdx.z]+fac*  (LNRHO[ll_8+i_8-1][vertexIdx.y][vertexIdx.z]-LNRHO[ll_8-i_8-1][vertexIdx.y][vertexIdx.z] + AC_dx2_bound__mod__cdata[-i_8+NGHOST+1-1]*-tmp_yz)
            }
          }
        }
      }
    }
  }
}
#endif

#if Leos_idealgas_MODULE && Lentropy_MODULE
bc_ss_flux_turb_x(AcBoundary boundary, AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  real rho_yz
  real cs2_yz
  real dsdx_yz
  real tt_yz
  if (topbot == AC_bot) {
    if (AC_pretend_lntt__mod__cdata) {
      SS[l1-1-1][vertexIdx.y][vertexIdx.z]=SS[1+l1-1][vertexIdx.y][vertexIdx.z] +  AC_dx2_bound__mod__cdata[-1+NGHOST+1-1]*AC_sigmasbt__mod__equationofstate*(exp(SS[l1-1][vertexIdx.y][vertexIdx.z])*exp(SS[l1-1][vertexIdx.y][vertexIdx.z])*exp(SS[l1-1][vertexIdx.y][vertexIdx.z]))/AC_hcondxbot__mod__energy
      SS[l1-2-1][vertexIdx.y][vertexIdx.z]=SS[2+l1-1][vertexIdx.y][vertexIdx.z] +  AC_dx2_bound__mod__cdata[-2+NGHOST+1-1]*AC_sigmasbt__mod__equationofstate*(exp(SS[l1-1][vertexIdx.y][vertexIdx.z])*exp(SS[l1-1][vertexIdx.y][vertexIdx.z])*exp(SS[l1-1][vertexIdx.y][vertexIdx.z]))/AC_hcondxbot__mod__energy
      SS[l1-3-1][vertexIdx.y][vertexIdx.z]=SS[3+l1-1][vertexIdx.y][vertexIdx.z] +  AC_dx2_bound__mod__cdata[-3+NGHOST+1-1]*AC_sigmasbt__mod__equationofstate*(exp(SS[l1-1][vertexIdx.y][vertexIdx.z])*exp(SS[l1-1][vertexIdx.y][vertexIdx.z])*exp(SS[l1-1][vertexIdx.y][vertexIdx.z]))/AC_hcondxbot__mod__energy
    }
    else {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          rho_yz=LNRHO[l1-1][vertexIdx.y][vertexIdx.z]+AC_reference_state__mod__density[1-1][iref_rho-1]
        }
        else {
          rho_yz=LNRHO[l1-1][vertexIdx.y][vertexIdx.z]
        }
      }
      else {
        rho_yz=exp(LNRHO[l1-1][vertexIdx.y][vertexIdx.z])
      }
      if (AC_ldensity_nolog__mod__cdata) {
        cs2_yz=SS[l1-1][vertexIdx.y][vertexIdx.z]
        if (AC_lreference_state__mod__cdata) {
          cs2_yz = cs2_yz+AC_reference_state__mod__density[1-1][iref_s-1]
        }
        cs2_yz=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(log(rho_yz)-AC_lnrho0__mod__equationofstate)+AC_cv1__mod__equationofstate*cs2_yz)
      }
      else {
        cs2_yz=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(LNRHO[l1-1][vertexIdx.y][vertexIdx.z]-AC_lnrho0__mod__equationofstate)+AC_cv1__mod__equationofstate*SS[l1-1][vertexIdx.y][vertexIdx.z])
      }
      tt_yz=cs2_yz/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
      dlnrhodx_yz= AC_coeffs_1_x__mod__cdata[NGHOST+1-1][1-1]*(LNRHO[1+l1-1][vertexIdx.y][vertexIdx.z]-LNRHO[l1-1-1][vertexIdx.y][vertexIdx.z])  +AC_coeffs_1_x__mod__cdata[NGHOST+2-1][1-1]*(LNRHO[2+l1-1][vertexIdx.y][vertexIdx.z]-LNRHO[l1-2-1][vertexIdx.y][vertexIdx.z])  +AC_coeffs_1_x__mod__cdata[NGHOST+3-1][1-1]*(LNRHO[3+l1-1][vertexIdx.y][vertexIdx.z]-LNRHO[l1-3-1][vertexIdx.y][vertexIdx.z])
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          dlnrhodx_yz=dlnrhodx_yz + AC_reference_state__mod__density[1-1][iref_grho-1]
          dlnrhodx_yz=dlnrhodx_yz/(rho_yz + AC_reference_state__mod__density[1-1][iref_rho-1])
        }
        else {
          dlnrhodx_yz=dlnrhodx_yz/rho_yz
        }
      }
      if (AC_lheatc_kramers__mod__energy) {
        dsdx_yz=-AC_cv__mod__equationofstate*( (AC_sigmasbt__mod__equationofstate/AC_hcond0_kramers__mod__energy)*pow(tt_yz,(3-6.5*AC_nkramers__mod__energy))*pow(rho_yz,(2.*AC_nkramers__mod__energy))  +AC_gamma_m1__mod__equationofstate*dlnrhodx_yz)
      }
      else {
        dsdx_yz=-(AC_sigmasbt__mod__equationofstate*(tt_yz*tt_yz*tt_yz)+AC_hcondxbot__mod__energy*AC_gamma_m1__mod__equationofstate*dlnrhodx_yz)/  (AC_chit_prof1__mod__energy*AC_chi_t__mod__energy*rho_yz+AC_hcondxbot__mod__energy/AC_cv__mod__equationofstate)
      }
      if (AC_lreference_state__mod__cdata) {
        dsdx_yz = dsdx_yz - AC_reference_state__mod__density[1-1][iref_gs-1]
      }
      SS[l1-1-1][vertexIdx.y][vertexIdx.z]=SS[1+l1-1][vertexIdx.y][vertexIdx.z]-AC_dx2_bound__mod__cdata[-1+NGHOST+1-1]*dsdx_yz
      SS[l1-2-1][vertexIdx.y][vertexIdx.z]=SS[2+l1-1][vertexIdx.y][vertexIdx.z]-AC_dx2_bound__mod__cdata[-2+NGHOST+1-1]*dsdx_yz
      SS[l1-3-1][vertexIdx.y][vertexIdx.z]=SS[3+l1-1][vertexIdx.y][vertexIdx.z]-AC_dx2_bound__mod__cdata[-3+NGHOST+1-1]*dsdx_yz
    }
  }
  else if (topbot == AC_top) {
    if (AC_pretend_lntt__mod__cdata) {
      SS[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z] +  AC_dx2_bound__mod__cdata[2+NGHOST-1]*AC_sigmasbt__mod__equationofstate*(exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])*exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])*exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]))/AC_hcondxtop__mod__energy
      SS[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z] +  AC_dx2_bound__mod__cdata[3+NGHOST-1]*AC_sigmasbt__mod__equationofstate*(exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])*exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])*exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]))/AC_hcondxtop__mod__energy
      SS[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z] +  AC_dx2_bound__mod__cdata[4+NGHOST-1]*AC_sigmasbt__mod__equationofstate*(exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])*exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])*exp(SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]))/AC_hcondxtop__mod__energy
    }
    else {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          rho_yz=LNRHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]+AC_reference_state__mod__density[nx-1][iref_rho-1]
        }
        else {
          rho_yz=LNRHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]
        }
      }
      else {
        rho_yz=exp(LNRHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
      }
      if (AC_ldensity_nolog__mod__cdata) {
        cs2_yz=SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]
        if (AC_lreference_state__mod__cdata) {
          cs2_yz = cs2_yz+AC_reference_state__mod__density[nx-1][iref_s-1]
        }
        cs2_yz=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(log(rho_yz)-AC_lnrho0__mod__equationofstate)+AC_cv1__mod__equationofstate*cs2_yz)
      }
      else {
        cs2_yz=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(LNRHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-AC_lnrho0__mod__equationofstate)+AC_cv1__mod__equationofstate*SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
      }
      tt_yz=cs2_yz/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
      dlnrhodx_yz= AC_coeffs_1_x__mod__cdata[NGHOST+1][2-1]*(LNRHO[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-LNRHO[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z])  +AC_coeffs_1_x__mod__cdata[NGHOST+2][2-1]*(LNRHO[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-LNRHO[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z])  +AC_coeffs_1_x__mod__cdata[NGHOST+3][2-1]*(LNRHO[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-LNRHO[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z])
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          dlnrhodx_yz=dlnrhodx_yz + AC_reference_state__mod__density[nx-1][iref_grho-1]
          dlnrhodx_yz=dlnrhodx_yz/(rho_yz + AC_reference_state__mod__density[nx-1][iref_rho-1])
        }
        else {
          dlnrhodx_yz=dlnrhodx_yz/rho_yz
        }
      }
      if (AC_lheatc_kramers__mod__energy  ||  (AC_hcondxtop__mod__energy != 0.0)) {
        hcond_total = AC_hcondxtop__mod__energy
        if (AC_lheatc_kramers__mod__energy) {
          hcond_total = hcond_total +  AC_hcond0_kramers__mod__energy*pow(tt_yz,(6.5*AC_nkramers__mod__energy))*pow(rho_yz,(-2.*AC_nkramers__mod__energy))
        }
        dsdx_yz = -(AC_sigmasbt__mod__equationofstate*(tt_yz*tt_yz*tt_yz)+hcond_total*AC_gamma_m1__mod__equationofstate*dlnrhodx_yz) /  (AC_chit_prof2__mod__energy*AC_chi_t__mod__energy*rho_yz+hcond_total/AC_cv__mod__equationofstate)
        if (AC_lreference_state__mod__cdata) {
          dsdx_yz = dsdx_yz - AC_reference_state__mod__density[nx-1][iref_gs-1]
        }
        SS[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]+AC_dx2_bound__mod__cdata[1+NGHOST]*dsdx_yz
        SS[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]+AC_dx2_bound__mod__cdata[2+NGHOST]*dsdx_yz
        SS[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]+AC_dx2_bound__mod__cdata[3+NGHOST]*dsdx_yz
      }
    }
  }
  else {
  }
}
#endif
/**

bc_ss_temp_x(AcBoundary boundary, AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  real tmp
  if (topbot == AC_bot) {
    tmp = 2/AC_gamma__mod__equationofstate*log(AC_cs2bot__mod__equationofstate/AC_cs20__mod__equationofstate)
    SS[l1-1][vertexIdx.y][vertexIdx.z] = 0.5*tmp - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[l1-1][vertexIdx.y][vertexIdx.z]-AC_lnrho0__mod__equationofstate)
    SS[l1-1-1][vertexIdx.y][vertexIdx.z] = -SS[1+l1-1][vertexIdx.y][vertexIdx.z] + tmp  - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[1+l1-1][vertexIdx.y][vertexIdx.z]+LNRHO[l1-1-1][vertexIdx.y][vertexIdx.z]-2*AC_lnrho0__mod__equationofstate)
    SS[l1-2-1][vertexIdx.y][vertexIdx.z] = -SS[2+l1-1][vertexIdx.y][vertexIdx.z] + tmp  - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[2+l1-1][vertexIdx.y][vertexIdx.z]+LNRHO[l1-2-1][vertexIdx.y][vertexIdx.z]-2*AC_lnrho0__mod__equationofstate)
    SS[l1-3-1][vertexIdx.y][vertexIdx.z] = -SS[3+l1-1][vertexIdx.y][vertexIdx.z] + tmp  - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[3+l1-1][vertexIdx.y][vertexIdx.z]+LNRHO[l1-3-1][vertexIdx.y][vertexIdx.z]-2*AC_lnrho0__mod__equationofstate)
  }
  else if (topbot == AC_top) {
    tmp = 2/AC_gamma__mod__equationofstate*log(AC_cs2top__mod__equationofstate/AC_cs20__mod__equationofstate)
    SS[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = 0.5*tmp - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-AC_lnrho0__mod__equationofstate)
    SS[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = -SS[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z] + tmp  - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]+LNRHO[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-2*AC_lnrho0__mod__equationofstate)
    SS[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = -SS[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z] + tmp  - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]+LNRHO[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-2*AC_lnrho0__mod__equationofstate)
    SS[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = -SS[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z] + tmp  - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]+LNRHO[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-2*AC_lnrho0__mod__equationofstate)
  }
  else {
  }
}

bc_db_x(AcBoundary boundary, AC_TOP_BOT topbot,Field j)
{
  suppress_unused_warning(boundary)
  if (topbot == AC_bot) {
    fder=(-3*j[l1-1+1-1][vertexIdx.y][vertexIdx.z]+4*j[l1-1+2-1][vertexIdx.y][vertexIdx.z] -j[l1-1+3-1][vertexIdx.y][vertexIdx.z])/(2*AC_dx__mod__cdata)
    j[l1-1-1][vertexIdx.y][vertexIdx.z]=j[l1-1+2-1][vertexIdx.y][vertexIdx.z]-2*AC_dx__mod__cdata*fder
    fder=(-3*j[l1-2+1-1][vertexIdx.y][vertexIdx.z]+4*j[l1-2+2-1][vertexIdx.y][vertexIdx.z] -j[l1-2+3-1][vertexIdx.y][vertexIdx.z])/(2*AC_dx__mod__cdata)
    j[l1-2-1][vertexIdx.y][vertexIdx.z]=j[l1-2+2-1][vertexIdx.y][vertexIdx.z]-2*AC_dx__mod__cdata*fder
    fder=(-3*j[l1-3+1-1][vertexIdx.y][vertexIdx.z]+4*j[l1-3+2-1][vertexIdx.y][vertexIdx.z] -j[l1-3+3-1][vertexIdx.y][vertexIdx.z])/(2*AC_dx__mod__cdata)
    j[l1-3-1][vertexIdx.y][vertexIdx.z]=j[l1-3+2-1][vertexIdx.y][vertexIdx.z]-2*AC_dx__mod__cdata*fder
  }
  else if (topbot == AC_top) {
    fder=(3*j[AC_l2__mod__cdata+1-1-1][vertexIdx.y][vertexIdx.z]-4*j[AC_l2__mod__cdata+1-2-1][vertexIdx.y][vertexIdx.z] +j[AC_l2__mod__cdata+1-3-1][vertexIdx.y][vertexIdx.z])/(2*AC_dx__mod__cdata)
    j[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=j[AC_l2__mod__cdata+1-2-1][vertexIdx.y][vertexIdx.z]+2*AC_dx__mod__cdata*fder
    fder=(3*j[AC_l2__mod__cdata+2-1-1][vertexIdx.y][vertexIdx.z]-4*j[AC_l2__mod__cdata+2-2-1][vertexIdx.y][vertexIdx.z] +j[AC_l2__mod__cdata+2-3-1][vertexIdx.y][vertexIdx.z])/(2*AC_dx__mod__cdata)
    j[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=j[AC_l2__mod__cdata+2-2-1][vertexIdx.y][vertexIdx.z]+2*AC_dx__mod__cdata*fder
    fder=(3*j[AC_l2__mod__cdata+3-1-1][vertexIdx.y][vertexIdx.z]-4*j[AC_l2__mod__cdata+3-2-1][vertexIdx.y][vertexIdx.z] +j[AC_l2__mod__cdata+3-3-1][vertexIdx.y][vertexIdx.z])/(2*AC_dx__mod__cdata)
    j[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=j[AC_l2__mod__cdata+3-2-1][vertexIdx.y][vertexIdx.z]+2*AC_dx__mod__cdata*fder
  }
  else {
  }
}

bc_db_z(AcBoundary boundary, AC_TOP_BOT topbot,Field j)
{
  suppress_unused_warning(boundary)
  real fder
  if (topbot == AC_bot) {
    fder=(-3*j[vertexIdx.x][vertexIdx.y][n1-1+1-1]+4*j[vertexIdx.x][vertexIdx.y][n1-1+2-1] -j[vertexIdx.x][vertexIdx.y][n1-1+3-1])/(2*AC_dz__mod__cdata)
    j[vertexIdx.x][vertexIdx.y][n1-1-1]=j[vertexIdx.x][vertexIdx.y][n1-1+2-1]-2*AC_dz__mod__cdata*fder
    fder=(-3*j[vertexIdx.x][vertexIdx.y][n1-2+1-1]+4*j[vertexIdx.x][vertexIdx.y][n1-2+2-1] -j[vertexIdx.x][vertexIdx.y][n1-2+3-1])/(2*AC_dz__mod__cdata)
    j[vertexIdx.x][vertexIdx.y][n1-2-1]=j[vertexIdx.x][vertexIdx.y][n1-2+2-1]-2*AC_dz__mod__cdata*fder
    fder=(-3*j[vertexIdx.x][vertexIdx.y][n1-3+1-1]+4*j[vertexIdx.x][vertexIdx.y][n1-3+2-1] -j[vertexIdx.x][vertexIdx.y][n1-3+3-1])/(2*AC_dz__mod__cdata)
    j[vertexIdx.x][vertexIdx.y][n1-3-1]=j[vertexIdx.x][vertexIdx.y][n1-3+2-1]-2*AC_dz__mod__cdata*fder
  }
  else if (topbot == AC_top) {
    fder=(3*j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+1-1-1]-4*j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+1-2-1] +j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+1-3-1])/(2*AC_dz__mod__cdata)
    j[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1]=j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+1-2-1]+2*AC_dz__mod__cdata*fder
    fder=(3*j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+2-1-1]-4*j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+2-2-1] +j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+2-3-1])/(2*AC_dz__mod__cdata)
    j[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1]=j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+2-2-1]+2*AC_dz__mod__cdata*fder
    fder=(3*j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+3-1-1]-4*j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+3-2-1] +j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+3-3-1])/(2*AC_dz__mod__cdata)
    j[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1]=j[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+3-2-1]+2*AC_dz__mod__cdata*fder
  }
  else {
  }
}

bc_lnrho_temp_z(AcBoundary boundary, AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  real tmp
  real lnrho_xy
  if (topbot == AC_bot) {
    tmp = AC_cv__mod__equationofstate*log(AC_cs2bot__mod__equationofstate/AC_cs20__mod__equationofstate)
    if (AC_ldensity_nolog__mod__cdata) {
      if (AC_lreference_state__mod__cdata) {
        lnrho_xy= log(LNRHO[vertexIdx.x][vertexIdx.y][n1-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
      }
      else {
        lnrho_xy=log(LNRHO[vertexIdx.x][vertexIdx.y][n1-1])
      }
    }
    else {
      lnrho_xy=LNRHO[vertexIdx.x][vertexIdx.y][n1-1]
    }
    SS[vertexIdx.x][vertexIdx.y][n1-1] = tmp - (AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*(lnrho_xy-AC_lnrho0__mod__equationofstate)
    if (AC_lreference_state__mod__cdata) {
      SS[vertexIdx.x][vertexIdx.y][n1-1] = SS[vertexIdx.x][vertexIdx.y][n1-1] - AC_reference_state_padded__mod__density[vertexIdx.x][iref_s-1]
    }
    SS[vertexIdx.x][vertexIdx.y][n1-1-1] = 2*SS[vertexIdx.x][vertexIdx.y][n1-1]-SS[vertexIdx.x][vertexIdx.y][1+n1-1]
    SS[vertexIdx.x][vertexIdx.y][n1-2-1] = 2*SS[vertexIdx.x][vertexIdx.y][n1-1]-SS[vertexIdx.x][vertexIdx.y][2+n1-1]
    SS[vertexIdx.x][vertexIdx.y][n1-3-1] = 2*SS[vertexIdx.x][vertexIdx.y][n1-1]-SS[vertexIdx.x][vertexIdx.y][3+n1-1]
    tmp=-AC_gravz__mod__gravity/AC_cs2bot__mod__equationofstate
    LNRHO[vertexIdx.x][vertexIdx.y][n1-1-1]=LNRHO[vertexIdx.x][vertexIdx.y][1+n1-1] + AC_cp1__mod__equationofstate*(SS[vertexIdx.x][vertexIdx.y][1+n1-1]-SS[vertexIdx.x][vertexIdx.y][n1-1-1])+AC_dz2_bound__mod__cdata[-1+NGHOST+1-1]*tmp
    LNRHO[vertexIdx.x][vertexIdx.y][n1-2-1]=LNRHO[vertexIdx.x][vertexIdx.y][2+n1-1] + AC_cp1__mod__equationofstate*(SS[vertexIdx.x][vertexIdx.y][2+n1-1]-SS[vertexIdx.x][vertexIdx.y][n1-2-1])+AC_dz2_bound__mod__cdata[-2+NGHOST+1-1]*tmp
    LNRHO[vertexIdx.x][vertexIdx.y][n1-3-1]=LNRHO[vertexIdx.x][vertexIdx.y][3+n1-1] + AC_cp1__mod__equationofstate*(SS[vertexIdx.x][vertexIdx.y][3+n1-1]-SS[vertexIdx.x][vertexIdx.y][n1-3-1])+AC_dz2_bound__mod__cdata[-3+NGHOST+1-1]*tmp
  }
  else if (topbot == AC_top) {
    tmp = AC_cv__mod__equationofstate*log(AC_cs2top__mod__equationofstate/AC_cs20__mod__equationofstate)
    if (AC_ldensity_nolog__mod__cdata) {
      if (AC_lreference_state__mod__cdata) {
        lnrho_xy= log(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
      }
      else {
        lnrho_xy=log(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1])
      }
    }
    else {
      lnrho_xy=LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]
    }
    SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1] = tmp - (AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*(lnrho_xy-AC_lnrho0__mod__equationofstate)
    if (AC_lreference_state__mod__cdata) {
      SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1] = SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1] - AC_reference_state_padded__mod__density[vertexIdx.x][iref_s-1]
    }
    SS[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1] = 2*SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]-SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1]
    SS[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1] = 2*SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]-SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1]
    SS[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1] = 2*SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]-SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1]
    tmp=AC_gravz__mod__gravity/AC_cs2top__mod__equationofstate
    LNRHO[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1]=LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1] + AC_cp1__mod__equationofstate*(SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1]-SS[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1])+AC_dz2_bound__mod__cdata[2+NGHOST-1]*tmp
    LNRHO[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1]=LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1] + AC_cp1__mod__equationofstate*(SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1]-SS[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1])+AC_dz2_bound__mod__cdata[3+NGHOST-1]*tmp
    LNRHO[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1]=LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1] + AC_cp1__mod__equationofstate*(SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1]-SS[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1])+AC_dz2_bound__mod__cdata[4+NGHOST-1]*tmp
  }
  else {
  }
}
bc_ss_flux(AcBoundary boundary, AC_TOP_BOT topbot,bool lone_sided)
{
  suppress_unused_warning(boundary)
  real tmp_xy
  real cs2_xy
  real rho_xy
  bool loptest_return_value_1
  bool loptest_return_value_4
  if (topbot == AC_bot) {
    if (AC_pretend_lntt__mod__cdata) {
      tmp_xy=-AC_fbotkbot__mod__energy/exp(SS[vertexIdx.x][vertexIdx.y][n1-1])
      SS[vertexIdx.x][vertexIdx.y][n1-1-1]=SS[vertexIdx.x][vertexIdx.y][1+n1-1]-AC_dz2_bound__mod__cdata[-1+NGHOST+1-1]*tmp_xy
      SS[vertexIdx.x][vertexIdx.y][n1-2-1]=SS[vertexIdx.x][vertexIdx.y][2+n1-1]-AC_dz2_bound__mod__cdata[-2+NGHOST+1-1]*tmp_xy
      SS[vertexIdx.x][vertexIdx.y][n1-3-1]=SS[vertexIdx.x][vertexIdx.y][3+n1-1]-AC_dz2_bound__mod__cdata[-3+NGHOST+1-1]*tmp_xy
    }
    else {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          rho_xy= LNRHO[vertexIdx.x][vertexIdx.y][n1-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1]
        }
        else {
          rho_xy=LNRHO[vertexIdx.x][vertexIdx.y][n1-1]
        }
      }
      else {
        rho_xy=exp(LNRHO[vertexIdx.x][vertexIdx.y][n1-1])
      }
      cs2_xy = SS[vertexIdx.x][vertexIdx.y][n1-1]
      if (AC_lreference_state__mod__cdata) {
        cs2_xy = cs2_xy + AC_reference_state_padded__mod__density[vertexIdx.x][iref_s-1]
      }
      if (AC_ldensity_nolog__mod__cdata) {
        cs2_xy=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(log(rho_xy)-AC_lnrho0__mod__equationofstate)+AC_cv1__mod__equationofstate*cs2_xy)
      }
      else {
        cs2_xy=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(LNRHO[vertexIdx.x][vertexIdx.y][n1-1]-AC_lnrho0__mod__equationofstate)+AC_cv1__mod__equationofstate*cs2_xy)
      }
      if (AC_lheatc_chiconst__mod__energy) {
        tmp_xy=AC_fbot__mod__energy/(rho_xy*AC_chi__mod__energy*cs2_xy)
      }
      else if (AC_lheatc_kramers__mod__energy) {
        tmp_xy=AC_fbot__mod__energy*pow(rho_xy,(2*AC_nkramers__mod__energy))*pow((AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate),(6.5*AC_nkramers__mod__energy))  /(AC_hcond0_kramers__mod__energy*pow(cs2_xy,(6.5*AC_nkramers__mod__energy+1.)))
      }
      else {
        tmp_xy=AC_cp__mod__equationofstate*AC_fbotkbot__mod__energy/cs2_xy
      }
      if (present(lone_sided)) {
        loptest_return_value_1=lone_sided
      }
      else if (false) {
      }
      else {
        loptest_return_value_1=false
      }
      if (loptest_return_value_1) {
        print("not implemented getderlnrho_z")
        print("not implemented bval_from_arr")
        print("not implemented set_ghosts_for_onesided_ders")
      }
      else {
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][1+n1-1]-LNRHO[vertexIdx.x][vertexIdx.y][n1-1-1]
        if (AC_ldensity_nolog__mod__cdata) {
          if (AC_lreference_state__mod__cdata) {
            rho_xy = rho_xy/(LNRHO[vertexIdx.x][vertexIdx.y][n1-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
          }
          else {
            rho_xy = rho_xy/LNRHO[vertexIdx.x][vertexIdx.y][n1-1]
          }
        }
        SS[vertexIdx.x][vertexIdx.y][n1-1-1]=SS[vertexIdx.x][vertexIdx.y][1+n1-1]+(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*(rho_xy+AC_dz2_bound__mod__cdata[-1+NGHOST+1-1]*tmp_xy)
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][2+n1-1]-LNRHO[vertexIdx.x][vertexIdx.y][n1-2-1]
        if (AC_ldensity_nolog__mod__cdata) {
          if (AC_lreference_state__mod__cdata) {
            rho_xy = rho_xy/(LNRHO[vertexIdx.x][vertexIdx.y][n1-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
          }
          else {
            rho_xy = rho_xy/LNRHO[vertexIdx.x][vertexIdx.y][n1-1]
          }
        }
        SS[vertexIdx.x][vertexIdx.y][n1-2-1]=SS[vertexIdx.x][vertexIdx.y][2+n1-1]+(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*(rho_xy+AC_dz2_bound__mod__cdata[-2+NGHOST+1-1]*tmp_xy)
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][3+n1-1]-LNRHO[vertexIdx.x][vertexIdx.y][n1-3-1]
        if (AC_ldensity_nolog__mod__cdata) {
          if (AC_lreference_state__mod__cdata) {
            rho_xy = rho_xy/(LNRHO[vertexIdx.x][vertexIdx.y][n1-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
          }
          else {
            rho_xy = rho_xy/LNRHO[vertexIdx.x][vertexIdx.y][n1-1]
          }
        }
        SS[vertexIdx.x][vertexIdx.y][n1-3-1]=SS[vertexIdx.x][vertexIdx.y][3+n1-1]+(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*(rho_xy+AC_dz2_bound__mod__cdata[-3+NGHOST+1-1]*tmp_xy)
      }
    }
  }
  else if (topbot == AC_top) {
    if (AC_pretend_lntt__mod__cdata) {
      tmp_xy=-AC_ftopktop__mod__energy/exp(SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1])
      SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1]=SS[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1]-AC_dz2_bound__mod__cdata[2+NGHOST-1]*tmp_xy
      SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1]=SS[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1]-AC_dz2_bound__mod__cdata[3+NGHOST-1]*tmp_xy
      SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1]=SS[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1]-AC_dz2_bound__mod__cdata[4+NGHOST-1]*tmp_xy
    }
    else {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          rho_xy= LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1]
        }
        else {
          rho_xy=LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]
        }
      }
      else {
        rho_xy=exp(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1])
      }
      cs2_xy = SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]
      if (AC_lreference_state__mod__cdata) {
        cs2_xy = cs2_xy + AC_reference_state_padded__mod__density[vertexIdx.x][iref_s-1]
      }
      if (AC_ldensity_nolog__mod__cdata) {
        cs2_xy=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(log(rho_xy)-AC_lnrho0__mod__equationofstate)+AC_cv1__mod__equationofstate*cs2_xy)
      }
      else {
        cs2_xy=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]-AC_lnrho0__mod__equationofstate)+AC_cv1__mod__equationofstate*cs2_xy)
      }
      if (AC_lheatc_chiconst__mod__energy) {
        tmp_xy=AC_ftop__mod__energy/(rho_xy*AC_chi__mod__energy*cs2_xy)
      }
      else if (AC_lheatc_kramers__mod__energy) {
        tmp_xy=AC_ftop__mod__energy*pow(rho_xy,(2*AC_nkramers__mod__energy))*pow((AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate),(6.5*AC_nkramers__mod__energy))  /(AC_hcond0_kramers__mod__energy*pow(cs2_xy,(6.5*AC_nkramers__mod__energy+1.)))
      }
      else {
        tmp_xy=AC_cp__mod__equationofstate*AC_ftopktop__mod__energy/cs2_xy
      }
      if (present(lone_sided)) {
        loptest_return_value_4=lone_sided
      }
      else if (false) {
      }
      else {
        loptest_return_value_4=false
      }
      if (loptest_return_value_4) {
        print("not implemented getderlnrho_z")
        print("not implemented bval_from_arr")
        print("not implemented set_ghosts_for_onesided_ders")
      }
      else {
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1]
        if (AC_ldensity_nolog__mod__cdata) {
          if (AC_lreference_state__mod__cdata) {
            rho_xy = rho_xy/(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
          }
          else {
            rho_xy = rho_xy/LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]
          }
        }
        SS[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1]=SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1]+(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*(-rho_xy-AC_dz2_bound__mod__cdata[2+NGHOST-1]*tmp_xy)
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1]
        if (AC_ldensity_nolog__mod__cdata) {
          if (AC_lreference_state__mod__cdata) {
            rho_xy = rho_xy/(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
          }
          else {
            rho_xy = rho_xy/LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]
          }
        }
        SS[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1]=SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1]+(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*(-rho_xy-AC_dz2_bound__mod__cdata[3+NGHOST-1]*tmp_xy)
        rho_xy = LNRHO[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1]
        if (AC_ldensity_nolog__mod__cdata) {
          if (AC_lreference_state__mod__cdata) {
            rho_xy = rho_xy/(LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
          }
          else {
            rho_xy = rho_xy/LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]
          }
        }
        SS[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1]=SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1]+(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*(-rho_xy-AC_dz2_bound__mod__cdata[4+NGHOST-1]*tmp_xy)
      }
    }
  }
  else {
  }
}
**/
#if Leos_idealgas_MODULE && Lgravity_simple_MODULE
bc_lnrho_hds_z_iso(AcBoundary boundary, AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  real cs2
  real rho
  real ss
  real dlnrhodz
  real dssdz
  real cs2_point
  real potp
  real potm
  real getrho_s_return_value_0
  real getrho_s_return_value_0
  real lnrho__1
  real ss__1
  real lntt__1
  real ee__1
  real pp__1
  real cs2__1
  real tt__1
  real rho_1
  real eth_1
  real lnrho__2
  real ss__2
  real lntt__2
  real ee__2
  real pp__2
  real cs2__2
  real tt__2
  real rho_2
  real eth_2
  real lnrho__3
  real ss__3
  real lntt__3
  real ee__3
  real pp__3
  real cs2__3
  real tt__3
  real rho_3
  real eth_3
  real lnrho__4
  real ss__4
  real lntt__4
  real ee__4
  real pp__4
  real cs2__4
  real tt__4
  real rho_4
  real eth_4
  real potx_xpoint_5
  real poty_ypoint_5
  real potz_zpoint_5
  real prof_5
  real xdep_5
  real zdep_5
  real potx_xpoint_6
  real poty_ypoint_6
  real potz_zpoint_6
  real prof_6
  real xdep_6
  real zdep_6
  real getrho_s_return_value_7
  real getrho_s_return_value_7
  real lnrho__8
  real ss__8
  real lntt__8
  real ee__8
  real pp__8
  real cs2__8
  real tt__8
  real rho_8
  real eth_8
  real lnrho__9
  real ss__9
  real lntt__9
  real ee__9
  real pp__9
  real cs2__9
  real tt__9
  real rho_9
  real eth_9
  real lnrho__10
  real ss__10
  real lntt__10
  real ee__10
  real pp__10
  real cs2__10
  real tt__10
  real rho_10
  real eth_10
  real potx_xpoint_11
  real poty_ypoint_11
  real potz_zpoint_11
  real prof_11
  real xdep_11
  real zdep_11
  real potx_xpoint_12
  real poty_ypoint_12
  real potz_zpoint_12
  real prof_12
  real xdep_12
  real zdep_12
  if (topbot == AC_bot) {
    if (lentropy) {
      if (ldensity) {
        if (AC_ldensity_nolog__mod__cdata) {
          if (AC_lreference_state__mod__cdata) {
            getrho_s_return_value_0=LNRHO[l1-1][m1-1][n1-1]+AC_reference_state__mod__density[l1-l1+1-1][iref_rho-1]
          }
          else {
            getrho_s_return_value_0=LNRHO[l1-1][m1-1][n1-1]
          }
        }
        else {
          getrho_s_return_value_0=exp(LNRHO[l1-1][m1-1][n1-1])
        }
        rho=getrho_s_return_value_0
        ss=SS[l1-1][m1-1][n1-1]
        if (AC_lreference_state__mod__cdata) {
          ss=ss+AC_reference_state__mod__density[AC_xbot__mod__equationofstate-1][iref_s-1]
        }
        if (AC_irho_ss__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate || AC_irho_ss__mod__equationofstate == AC_irho_ss__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else if (AC_irho_ss__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate)  {
            lnrho__1 = rho
          }
          else {
            lnrho__1 = log(rho)
          }
          ss__1=ss
          lntt__1=AC_lntt0__mod__equationofstate+AC_cv1__mod__equationofstate*ss__1+AC_gamma_m1__mod__equationofstate*(lnrho__1-AC_lnrho0__mod__equationofstate)
          ee__1=AC_cv__mod__equationofstate*exp(lntt__1)
          pp__1=(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*exp(lntt__1+lnrho__1)
          cs2__1=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__1
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ilnrho_ee__mod__equationofstate || AC_irho_ss__mod__equationofstate == AC_irho_ee__mod__equationofstate) {
          ee__1=ss
          if (AC_irho_ss__mod__equationofstate==AC_ilnrho_ee__mod__equationofstate) {
            lnrho__1=rho
            pp__1=AC_gamma_m1__mod__equationofstate*ee__1*exp(rho)
          }
          else {
            lnrho__1=log(rho)
            pp__1=AC_gamma_m1__mod__equationofstate*ee__1*rho
          }
          lntt__1=log(AC_cv1__mod__equationofstate*ee__1)
          ss__1=AC_cv__mod__equationofstate*(lntt__1-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__1-AC_lnrho0__mod__equationofstate))
          cs2__1=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__1
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ilnrho_pp__mod__equationofstate || AC_irho_ss__mod__equationofstate == AC_irho_pp__mod__equationofstate) {
          pp__1=ss
          if (AC_irho_ss__mod__equationofstate==AC_ilnrho_pp__mod__equationofstate) {
            lnrho__1=rho
            ee__1=pp__1*exp(-lnrho__1)/AC_gamma_m1__mod__equationofstate
            ss__1=AC_cv__mod__equationofstate*(log(pp__1*exp(-lnrho__1)*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__1-AC_lnrho0__mod__equationofstate))
          }
          else {
            lnrho__1=log(rho)
            ee__1=pp__1/rho/AC_gamma_m1__mod__equationofstate
            ss__1=AC_cv__mod__equationofstate*(log(pp__1/rho*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__1-AC_lnrho0__mod__equationofstate))
          }
          lntt__1=log(AC_cv1__mod__equationofstate*ee__1)
          cs2__1=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__1
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ilnrho_lntt__mod__equationofstate) {
          lnrho__1=rho
          lntt__1=ss
          ss__1=AC_cv__mod__equationofstate*(lntt__1-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__1-AC_lnrho0__mod__equationofstate))
          ee__1=AC_cv__mod__equationofstate*exp(lntt__1)
          pp__1=ee__1*exp(lnrho__1)*AC_gamma_m1__mod__equationofstate
          cs2__1=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__1
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ilnrho_tt__mod__equationofstate) {
          lnrho__1=rho
          tt__1=ss
          ss__1=AC_cv__mod__equationofstate*(log(tt__1)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__1-AC_lnrho0__mod__equationofstate))
          ee__1=AC_cv__mod__equationofstate*tt__1
          pp__1=ee__1*exp(lnrho__1)*AC_gamma_m1__mod__equationofstate
          cs2__1=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__1
        }
        else if (AC_irho_ss__mod__equationofstate == AC_irho_tt__mod__equationofstate) {
          lnrho__1=log(rho)
          tt__1=ss
          ss__1=AC_cv__mod__equationofstate*(log(tt__1)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__1-AC_lnrho0__mod__equationofstate))
          ee__1=AC_cv__mod__equationofstate*tt__1
          pp__1=ee__1*rho*AC_gamma_m1__mod__equationofstate
          cs2__1=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__1
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ipp_cs2__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__1=log(rho)
              tt__1=exp(AC_lntt0__mod__equationofstate)
              pp__1=exp(lnrho__1)*AC_cs20__mod__equationofstate/AC_gamma__mod__equationofstate
            }
            else {
              if (AC_leos_isothermal__mod__equationofstate) {
                pp__1=rho
                lnrho__1=log(pp__1*AC_cs20__mod__equationofstate)
                tt__1=exp(AC_lntt0__mod__equationofstate)
              }
            }
          }
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ipp_ss__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__1=(rho)
              ss__1=ss
              cs2__1=exp(AC_gamma__mod__equationofstate*ss__1*AC_cp1__mod__equationofstate+AC_gamma_m1__mod__equationofstate*(lnrho__1-AC_lnrho0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__1=cs2__1/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
            }
            else {
              pp__1=rho
              ss__1=ss
              cs2__1=exp(ss__1*AC_cp1__mod__equationofstate+AC_gamma1__mod__equationofstate*AC_gamma_m1__mod__equationofstate*log(pp__1/AC_pp0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__1=cs2__1/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
              lnrho__1=log(AC_gamma__mod__equationofstate*pp__1/cs2__1)
            }
          }
        }
        else if (AC_irho_ss__mod__equationofstate == AC_irho_eth__mod__equationofstate || AC_irho_ss__mod__equationofstate == AC_ilnrho_eth__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else {
            if (AC_ldensity_nolog__mod__cdata) {
              rho_1 = rho
              if (false) {
                lnrho__1 = log(rho)
              }
            }
            else {
              rho_1 = exp(rho)
              if (false) {
                lnrho__1 = rho
              }
            }
            eth_1 = ss
          }
          if (false) {
            lntt__1 = log(AC_cv1__mod__equationofstate * eth_1/ rho_1)
          }
          if (false) {
            ee__1 = eth_1/ rho_1
          }
          if (false) {
            pp__1 = AC_gamma_m1__mod__equationofstate * eth_1
          }
          if (present(cs2_point)) {
            cs2__1 = AC_gamma__mod__equationofstate * AC_gamma_m1__mod__equationofstate * eth_1/ rho_1
          }
        }
        else {
        }
        if (present(cs2_point)) {
          cs2_point=cs2__1
        }
        dlnrhodz = AC_gamma__mod__equationofstate *AC_gravz__mod__gravity/cs2_point
        if (AC_ldensity_nolog__mod__cdata) {
          dlnrhodz=dlnrhodz*rho
        }
        dssdz = -AC_gamma_m1__mod__equationofstate*AC_gravz__mod__gravity/cs2_point
        LNRHO[vertexIdx.x][vertexIdx.y][n1-1-1] = LNRHO[vertexIdx.x][vertexIdx.y][1+n1-1] - AC_dz2_bound__mod__cdata[-1+NGHOST+1-1]*dlnrhodz
        SS[vertexIdx.x][vertexIdx.y][n1-1-1] = SS[vertexIdx.x][vertexIdx.y][1+n1-1] - AC_dz2_bound__mod__cdata[-1+NGHOST+1-1]*dssdz
        LNRHO[vertexIdx.x][vertexIdx.y][n1-2-1] = LNRHO[vertexIdx.x][vertexIdx.y][2+n1-1] - AC_dz2_bound__mod__cdata[-2+NGHOST+1-1]*dlnrhodz
        SS[vertexIdx.x][vertexIdx.y][n1-2-1] = SS[vertexIdx.x][vertexIdx.y][2+n1-1] - AC_dz2_bound__mod__cdata[-2+NGHOST+1-1]*dssdz
        LNRHO[vertexIdx.x][vertexIdx.y][n1-3-1] = LNRHO[vertexIdx.x][vertexIdx.y][3+n1-1] - AC_dz2_bound__mod__cdata[-3+NGHOST+1-1]*dlnrhodz
        SS[vertexIdx.x][vertexIdx.y][n1-3-1] = SS[vertexIdx.x][vertexIdx.y][3+n1-1] - AC_dz2_bound__mod__cdata[-3+NGHOST+1-1]*dssdz
      }
      else if (lanelastic) {
        if (AC_ipp_ss__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate || AC_ipp_ss__mod__equationofstate == AC_irho_ss__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else if (AC_ipp_ss__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate)  {
            lnrho__2 = log(RHO_B[l1-1][m1-1][n1-1])
          }
          else {
            lnrho__2 = log(log(RHO_B[l1-1][m1-1][n1-1]))
          }
          ss__2=SS_B[l1-1][m1-1][n1-1]
          lntt__2=AC_lntt0__mod__equationofstate+AC_cv1__mod__equationofstate*ss__2+AC_gamma_m1__mod__equationofstate*(lnrho__2-AC_lnrho0__mod__equationofstate)
          ee__2=AC_cv__mod__equationofstate*exp(lntt__2)
          pp__2=(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*exp(lntt__2+lnrho__2)
          cs2__2=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__2
        }
        else if (AC_ipp_ss__mod__equationofstate == AC_ilnrho_ee__mod__equationofstate || AC_ipp_ss__mod__equationofstate == AC_irho_ee__mod__equationofstate) {
          ee__2=SS_B[l1-1][m1-1][n1-1]
          if (AC_ipp_ss__mod__equationofstate==AC_ilnrho_ee__mod__equationofstate) {
            lnrho__2=log(RHO_B[l1-1][m1-1][n1-1])
            pp__2=AC_gamma_m1__mod__equationofstate*ee__2*exp(log(RHO_B[l1-1][m1-1][n1-1]))
          }
          else {
            lnrho__2=log(log(RHO_B[l1-1][m1-1][n1-1]))
            pp__2=AC_gamma_m1__mod__equationofstate*ee__2*log(RHO_B[l1-1][m1-1][n1-1])
          }
          lntt__2=log(AC_cv1__mod__equationofstate*ee__2)
          ss__2=AC_cv__mod__equationofstate*(lntt__2-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__2-AC_lnrho0__mod__equationofstate))
          cs2__2=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__2
        }
        else if (AC_ipp_ss__mod__equationofstate == AC_ilnrho_pp__mod__equationofstate || AC_ipp_ss__mod__equationofstate == AC_irho_pp__mod__equationofstate) {
          pp__2=SS_B[l1-1][m1-1][n1-1]
          if (AC_ipp_ss__mod__equationofstate==AC_ilnrho_pp__mod__equationofstate) {
            lnrho__2=log(RHO_B[l1-1][m1-1][n1-1])
            ee__2=pp__2*exp(-lnrho__2)/AC_gamma_m1__mod__equationofstate
            ss__2=AC_cv__mod__equationofstate*(log(pp__2*exp(-lnrho__2)*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__2-AC_lnrho0__mod__equationofstate))
          }
          else {
            lnrho__2=log(log(RHO_B[l1-1][m1-1][n1-1]))
            ee__2=pp__2/log(RHO_B[l1-1][m1-1][n1-1])/AC_gamma_m1__mod__equationofstate
            ss__2=AC_cv__mod__equationofstate*(log(pp__2/log(RHO_B[l1-1][m1-1][n1-1])*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__2-AC_lnrho0__mod__equationofstate))
          }
          lntt__2=log(AC_cv1__mod__equationofstate*ee__2)
          cs2__2=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__2
        }
        else if (AC_ipp_ss__mod__equationofstate == AC_ilnrho_lntt__mod__equationofstate) {
          lnrho__2=log(RHO_B[l1-1][m1-1][n1-1])
          lntt__2=SS_B[l1-1][m1-1][n1-1]
          ss__2=AC_cv__mod__equationofstate*(lntt__2-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__2-AC_lnrho0__mod__equationofstate))
          ee__2=AC_cv__mod__equationofstate*exp(lntt__2)
          pp__2=ee__2*exp(lnrho__2)*AC_gamma_m1__mod__equationofstate
          cs2__2=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__2
        }
        else if (AC_ipp_ss__mod__equationofstate == AC_ilnrho_tt__mod__equationofstate) {
          lnrho__2=log(RHO_B[l1-1][m1-1][n1-1])
          tt__2=SS_B[l1-1][m1-1][n1-1]
          ss__2=AC_cv__mod__equationofstate*(log(tt__2)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__2-AC_lnrho0__mod__equationofstate))
          ee__2=AC_cv__mod__equationofstate*tt__2
          pp__2=ee__2*exp(lnrho__2)*AC_gamma_m1__mod__equationofstate
          cs2__2=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__2
        }
        else if (AC_ipp_ss__mod__equationofstate == AC_irho_tt__mod__equationofstate) {
          lnrho__2=log(log(RHO_B[l1-1][m1-1][n1-1]))
          tt__2=SS_B[l1-1][m1-1][n1-1]
          ss__2=AC_cv__mod__equationofstate*(log(tt__2)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__2-AC_lnrho0__mod__equationofstate))
          ee__2=AC_cv__mod__equationofstate*tt__2
          pp__2=ee__2*log(RHO_B[l1-1][m1-1][n1-1])*AC_gamma_m1__mod__equationofstate
          cs2__2=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__2
        }
        else if (AC_ipp_ss__mod__equationofstate == AC_ipp_cs2__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__2=log(log(RHO_B[l1-1][m1-1][n1-1]))
              tt__2=exp(AC_lntt0__mod__equationofstate)
              pp__2=exp(lnrho__2)*AC_cs20__mod__equationofstate/AC_gamma__mod__equationofstate
            }
            else {
              if (AC_leos_isothermal__mod__equationofstate) {
                pp__2=log(RHO_B[l1-1][m1-1][n1-1])
                lnrho__2=log(pp__2*AC_cs20__mod__equationofstate)
                tt__2=exp(AC_lntt0__mod__equationofstate)
              }
            }
          }
        }
        else if (AC_ipp_ss__mod__equationofstate == AC_ipp_ss__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__2=(log(RHO_B[l1-1][m1-1][n1-1]))
              ss__2=SS_B[l1-1][m1-1][n1-1]
              cs2__2=exp(AC_gamma__mod__equationofstate*ss__2*AC_cp1__mod__equationofstate+AC_gamma_m1__mod__equationofstate*(lnrho__2-AC_lnrho0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__2=cs2__2/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
            }
            else {
              pp__2=log(RHO_B[l1-1][m1-1][n1-1])
              ss__2=SS_B[l1-1][m1-1][n1-1]
              cs2__2=exp(ss__2*AC_cp1__mod__equationofstate+AC_gamma1__mod__equationofstate*AC_gamma_m1__mod__equationofstate*log(pp__2/AC_pp0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__2=cs2__2/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
              lnrho__2=log(AC_gamma__mod__equationofstate*pp__2/cs2__2)
            }
          }
        }
        else if (AC_ipp_ss__mod__equationofstate == AC_irho_eth__mod__equationofstate || AC_ipp_ss__mod__equationofstate == AC_ilnrho_eth__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else {
            if (AC_ldensity_nolog__mod__cdata) {
              rho_2 = log(RHO_B[l1-1][m1-1][n1-1])
              if (false) {
                lnrho__2 = log(log(RHO_B[l1-1][m1-1][n1-1]))
              }
            }
            else {
              rho_2 = exp(log(RHO_B[l1-1][m1-1][n1-1]))
              if (false) {
                lnrho__2 = log(RHO_B[l1-1][m1-1][n1-1])
              }
            }
            eth_2 = SS_B[l1-1][m1-1][n1-1]
          }
          if (false) {
            lntt__2 = log(AC_cv1__mod__equationofstate * eth_2/ rho_2)
          }
          if (false) {
            ee__2 = eth_2/ rho_2
          }
          if (false) {
            pp__2 = AC_gamma_m1__mod__equationofstate * eth_2
          }
          if (present(cs2_point)) {
            cs2__2 = AC_gamma__mod__equationofstate * AC_gamma_m1__mod__equationofstate * eth_2/ rho_2
          }
        }
        else {
        }
        if (present(cs2_point)) {
          cs2_point=cs2__2
        }
        dlnrhodz = AC_gamma__mod__equationofstate *AC_gravz__mod__gravity/cs2_point
        dssdz    = AC_gamma_m1__mod__equationofstate*AC_gravz__mod__gravity/cs2_point
        RHO_B[vertexIdx.x][vertexIdx.y][n1-1-1] = RHO_B[vertexIdx.x][vertexIdx.y][1+n1-1] - AC_dz2_bound__mod__cdata[-1+NGHOST+1-1]*dlnrhodz*RHO_B[vertexIdx.x][vertexIdx.y][1+n1-1]
        SS_B[vertexIdx.x][vertexIdx.y][n1-1-1] = SS_B[vertexIdx.x][vertexIdx.y][1+n1-1] - AC_dz2_bound__mod__cdata[-1+NGHOST+1-1]*dssdz
        RHO_B[vertexIdx.x][vertexIdx.y][n1-2-1] = RHO_B[vertexIdx.x][vertexIdx.y][2+n1-1] - AC_dz2_bound__mod__cdata[-2+NGHOST+1-1]*dlnrhodz*RHO_B[vertexIdx.x][vertexIdx.y][1+n1-1]
        SS_B[vertexIdx.x][vertexIdx.y][n1-2-1] = SS_B[vertexIdx.x][vertexIdx.y][2+n1-1] - AC_dz2_bound__mod__cdata[-2+NGHOST+1-1]*dssdz
        RHO_B[vertexIdx.x][vertexIdx.y][n1-3-1] = RHO_B[vertexIdx.x][vertexIdx.y][3+n1-1] - AC_dz2_bound__mod__cdata[-3+NGHOST+1-1]*dlnrhodz*RHO_B[vertexIdx.x][vertexIdx.y][1+n1-1]
        SS_B[vertexIdx.x][vertexIdx.y][n1-3-1] = SS_B[vertexIdx.x][vertexIdx.y][3+n1-1] - AC_dz2_bound__mod__cdata[-3+NGHOST+1-1]*dssdz
      }
      else {
      }
    }
    else if (ltemperature) {
      if (AC_ltemperature_nolog__mod__cdata) {
        if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate || AC_ilnrho_tt__mod__equationofstate == AC_irho_ss__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate)  {
            lnrho__3 = LNRHO[l1-1][m1-1][n1-1]
          }
          else {
            lnrho__3 = log(LNRHO[l1-1][m1-1][n1-1])
          }
          ss__3=TT[l1-1][m1-1][n1-1]
          lntt__3=AC_lntt0__mod__equationofstate+AC_cv1__mod__equationofstate*ss__3+AC_gamma_m1__mod__equationofstate*(lnrho__3-AC_lnrho0__mod__equationofstate)
          ee__3=AC_cv__mod__equationofstate*exp(lntt__3)
          pp__3=(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*exp(lntt__3+lnrho__3)
          cs2__3=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__3
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_ee__mod__equationofstate || AC_ilnrho_tt__mod__equationofstate == AC_irho_ee__mod__equationofstate) {
          ee__3=TT[l1-1][m1-1][n1-1]
          if (AC_ilnrho_tt__mod__equationofstate==AC_ilnrho_ee__mod__equationofstate) {
            lnrho__3=LNRHO[l1-1][m1-1][n1-1]
            pp__3=AC_gamma_m1__mod__equationofstate*ee__3*exp(LNRHO[l1-1][m1-1][n1-1])
          }
          else {
            lnrho__3=log(LNRHO[l1-1][m1-1][n1-1])
            pp__3=AC_gamma_m1__mod__equationofstate*ee__3*LNRHO[l1-1][m1-1][n1-1]
          }
          lntt__3=log(AC_cv1__mod__equationofstate*ee__3)
          ss__3=AC_cv__mod__equationofstate*(lntt__3-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__3-AC_lnrho0__mod__equationofstate))
          cs2__3=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__3
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_pp__mod__equationofstate || AC_ilnrho_tt__mod__equationofstate == AC_irho_pp__mod__equationofstate) {
          pp__3=TT[l1-1][m1-1][n1-1]
          if (AC_ilnrho_tt__mod__equationofstate==AC_ilnrho_pp__mod__equationofstate) {
            lnrho__3=LNRHO[l1-1][m1-1][n1-1]
            ee__3=pp__3*exp(-lnrho__3)/AC_gamma_m1__mod__equationofstate
            ss__3=AC_cv__mod__equationofstate*(log(pp__3*exp(-lnrho__3)*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__3-AC_lnrho0__mod__equationofstate))
          }
          else {
            lnrho__3=log(LNRHO[l1-1][m1-1][n1-1])
            ee__3=pp__3/LNRHO[l1-1][m1-1][n1-1]/AC_gamma_m1__mod__equationofstate
            ss__3=AC_cv__mod__equationofstate*(log(pp__3/LNRHO[l1-1][m1-1][n1-1]*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__3-AC_lnrho0__mod__equationofstate))
          }
          lntt__3=log(AC_cv1__mod__equationofstate*ee__3)
          cs2__3=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__3
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_lntt__mod__equationofstate) {
          lnrho__3=LNRHO[l1-1][m1-1][n1-1]
          lntt__3=TT[l1-1][m1-1][n1-1]
          ss__3=AC_cv__mod__equationofstate*(lntt__3-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__3-AC_lnrho0__mod__equationofstate))
          ee__3=AC_cv__mod__equationofstate*exp(lntt__3)
          pp__3=ee__3*exp(lnrho__3)*AC_gamma_m1__mod__equationofstate
          cs2__3=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__3
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_tt__mod__equationofstate) {
          lnrho__3=LNRHO[l1-1][m1-1][n1-1]
          tt__3=TT[l1-1][m1-1][n1-1]
          ss__3=AC_cv__mod__equationofstate*(log(tt__3)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__3-AC_lnrho0__mod__equationofstate))
          ee__3=AC_cv__mod__equationofstate*tt__3
          pp__3=ee__3*exp(lnrho__3)*AC_gamma_m1__mod__equationofstate
          cs2__3=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__3
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_irho_tt__mod__equationofstate) {
          lnrho__3=log(LNRHO[l1-1][m1-1][n1-1])
          tt__3=TT[l1-1][m1-1][n1-1]
          ss__3=AC_cv__mod__equationofstate*(log(tt__3)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__3-AC_lnrho0__mod__equationofstate))
          ee__3=AC_cv__mod__equationofstate*tt__3
          pp__3=ee__3*LNRHO[l1-1][m1-1][n1-1]*AC_gamma_m1__mod__equationofstate
          cs2__3=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__3
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ipp_cs2__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__3=log(LNRHO[l1-1][m1-1][n1-1])
              tt__3=exp(AC_lntt0__mod__equationofstate)
              pp__3=exp(lnrho__3)*AC_cs20__mod__equationofstate/AC_gamma__mod__equationofstate
            }
            else {
              if (AC_leos_isothermal__mod__equationofstate) {
                pp__3=LNRHO[l1-1][m1-1][n1-1]
                lnrho__3=log(pp__3*AC_cs20__mod__equationofstate)
                tt__3=exp(AC_lntt0__mod__equationofstate)
              }
            }
          }
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ipp_ss__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__3=(LNRHO[l1-1][m1-1][n1-1])
              ss__3=TT[l1-1][m1-1][n1-1]
              cs2__3=exp(AC_gamma__mod__equationofstate*ss__3*AC_cp1__mod__equationofstate+AC_gamma_m1__mod__equationofstate*(lnrho__3-AC_lnrho0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__3=cs2__3/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
            }
            else {
              pp__3=LNRHO[l1-1][m1-1][n1-1]
              ss__3=TT[l1-1][m1-1][n1-1]
              cs2__3=exp(ss__3*AC_cp1__mod__equationofstate+AC_gamma1__mod__equationofstate*AC_gamma_m1__mod__equationofstate*log(pp__3/AC_pp0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__3=cs2__3/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
              lnrho__3=log(AC_gamma__mod__equationofstate*pp__3/cs2__3)
            }
          }
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_irho_eth__mod__equationofstate || AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_eth__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else {
            if (AC_ldensity_nolog__mod__cdata) {
              rho_3 = LNRHO[l1-1][m1-1][n1-1]
              if (false) {
                lnrho__3 = log(LNRHO[l1-1][m1-1][n1-1])
              }
            }
            else {
              rho_3 = exp(LNRHO[l1-1][m1-1][n1-1])
              if (false) {
                lnrho__3 = LNRHO[l1-1][m1-1][n1-1]
              }
            }
            eth_3 = TT[l1-1][m1-1][n1-1]
          }
          if (false) {
            lntt__3 = log(AC_cv1__mod__equationofstate * eth_3/ rho_3)
          }
          if (false) {
            ee__3 = eth_3/ rho_3
          }
          if (false) {
            pp__3 = AC_gamma_m1__mod__equationofstate * eth_3
          }
          if (present(cs2_point)) {
            cs2__3 = AC_gamma__mod__equationofstate * AC_gamma_m1__mod__equationofstate * eth_3/ rho_3
          }
        }
        else {
        }
        if (present(cs2_point)) {
          cs2_point=cs2__3
        }
      }
      else {
        if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate || AC_ilnrho_lntt__mod__equationofstate == AC_irho_ss__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate)  {
            lnrho__4 = LNRHO[l1-1][m1-1][n1-1]
          }
          else {
            lnrho__4 = log(LNRHO[l1-1][m1-1][n1-1])
          }
          ss__4=LNTT[l1-1][m1-1][n1-1]
          lntt__4=AC_lntt0__mod__equationofstate+AC_cv1__mod__equationofstate*ss__4+AC_gamma_m1__mod__equationofstate*(lnrho__4-AC_lnrho0__mod__equationofstate)
          ee__4=AC_cv__mod__equationofstate*exp(lntt__4)
          pp__4=(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*exp(lntt__4+lnrho__4)
          cs2__4=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__4
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_ee__mod__equationofstate || AC_ilnrho_lntt__mod__equationofstate == AC_irho_ee__mod__equationofstate) {
          ee__4=LNTT[l1-1][m1-1][n1-1]
          if (AC_ilnrho_lntt__mod__equationofstate==AC_ilnrho_ee__mod__equationofstate) {
            lnrho__4=LNRHO[l1-1][m1-1][n1-1]
            pp__4=AC_gamma_m1__mod__equationofstate*ee__4*exp(LNRHO[l1-1][m1-1][n1-1])
          }
          else {
            lnrho__4=log(LNRHO[l1-1][m1-1][n1-1])
            pp__4=AC_gamma_m1__mod__equationofstate*ee__4*LNRHO[l1-1][m1-1][n1-1]
          }
          lntt__4=log(AC_cv1__mod__equationofstate*ee__4)
          ss__4=AC_cv__mod__equationofstate*(lntt__4-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__4-AC_lnrho0__mod__equationofstate))
          cs2__4=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__4
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_pp__mod__equationofstate || AC_ilnrho_lntt__mod__equationofstate == AC_irho_pp__mod__equationofstate) {
          pp__4=LNTT[l1-1][m1-1][n1-1]
          if (AC_ilnrho_lntt__mod__equationofstate==AC_ilnrho_pp__mod__equationofstate) {
            lnrho__4=LNRHO[l1-1][m1-1][n1-1]
            ee__4=pp__4*exp(-lnrho__4)/AC_gamma_m1__mod__equationofstate
            ss__4=AC_cv__mod__equationofstate*(log(pp__4*exp(-lnrho__4)*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__4-AC_lnrho0__mod__equationofstate))
          }
          else {
            lnrho__4=log(LNRHO[l1-1][m1-1][n1-1])
            ee__4=pp__4/LNRHO[l1-1][m1-1][n1-1]/AC_gamma_m1__mod__equationofstate
            ss__4=AC_cv__mod__equationofstate*(log(pp__4/LNRHO[l1-1][m1-1][n1-1]*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__4-AC_lnrho0__mod__equationofstate))
          }
          lntt__4=log(AC_cv1__mod__equationofstate*ee__4)
          cs2__4=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__4
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_lntt__mod__equationofstate) {
          lnrho__4=LNRHO[l1-1][m1-1][n1-1]
          lntt__4=LNTT[l1-1][m1-1][n1-1]
          ss__4=AC_cv__mod__equationofstate*(lntt__4-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__4-AC_lnrho0__mod__equationofstate))
          ee__4=AC_cv__mod__equationofstate*exp(lntt__4)
          pp__4=ee__4*exp(lnrho__4)*AC_gamma_m1__mod__equationofstate
          cs2__4=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__4
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_tt__mod__equationofstate) {
          lnrho__4=LNRHO[l1-1][m1-1][n1-1]
          tt__4=LNTT[l1-1][m1-1][n1-1]
          ss__4=AC_cv__mod__equationofstate*(log(tt__4)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__4-AC_lnrho0__mod__equationofstate))
          ee__4=AC_cv__mod__equationofstate*tt__4
          pp__4=ee__4*exp(lnrho__4)*AC_gamma_m1__mod__equationofstate
          cs2__4=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__4
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_irho_tt__mod__equationofstate) {
          lnrho__4=log(LNRHO[l1-1][m1-1][n1-1])
          tt__4=LNTT[l1-1][m1-1][n1-1]
          ss__4=AC_cv__mod__equationofstate*(log(tt__4)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__4-AC_lnrho0__mod__equationofstate))
          ee__4=AC_cv__mod__equationofstate*tt__4
          pp__4=ee__4*LNRHO[l1-1][m1-1][n1-1]*AC_gamma_m1__mod__equationofstate
          cs2__4=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__4
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ipp_cs2__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__4=log(LNRHO[l1-1][m1-1][n1-1])
              tt__4=exp(AC_lntt0__mod__equationofstate)
              pp__4=exp(lnrho__4)*AC_cs20__mod__equationofstate/AC_gamma__mod__equationofstate
            }
            else {
              if (AC_leos_isothermal__mod__equationofstate) {
                pp__4=LNRHO[l1-1][m1-1][n1-1]
                lnrho__4=log(pp__4*AC_cs20__mod__equationofstate)
                tt__4=exp(AC_lntt0__mod__equationofstate)
              }
            }
          }
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ipp_ss__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__4=(LNRHO[l1-1][m1-1][n1-1])
              ss__4=LNTT[l1-1][m1-1][n1-1]
              cs2__4=exp(AC_gamma__mod__equationofstate*ss__4*AC_cp1__mod__equationofstate+AC_gamma_m1__mod__equationofstate*(lnrho__4-AC_lnrho0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__4=cs2__4/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
            }
            else {
              pp__4=LNRHO[l1-1][m1-1][n1-1]
              ss__4=LNTT[l1-1][m1-1][n1-1]
              cs2__4=exp(ss__4*AC_cp1__mod__equationofstate+AC_gamma1__mod__equationofstate*AC_gamma_m1__mod__equationofstate*log(pp__4/AC_pp0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__4=cs2__4/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
              lnrho__4=log(AC_gamma__mod__equationofstate*pp__4/cs2__4)
            }
          }
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_irho_eth__mod__equationofstate || AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_eth__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else {
            if (AC_ldensity_nolog__mod__cdata) {
              rho_4 = LNRHO[l1-1][m1-1][n1-1]
              if (false) {
                lnrho__4 = log(LNRHO[l1-1][m1-1][n1-1])
              }
            }
            else {
              rho_4 = exp(LNRHO[l1-1][m1-1][n1-1])
              if (false) {
                lnrho__4 = LNRHO[l1-1][m1-1][n1-1]
              }
            }
            eth_4 = LNTT[l1-1][m1-1][n1-1]
          }
          if (false) {
            lntt__4 = log(AC_cv1__mod__equationofstate * eth_4/ rho_4)
          }
          if (false) {
            ee__4 = eth_4/ rho_4
          }
          if (false) {
            pp__4 = AC_gamma_m1__mod__equationofstate * eth_4
          }
          if (present(cs2_point)) {
            cs2__4 = AC_gamma__mod__equationofstate * AC_gamma_m1__mod__equationofstate * eth_4/ rho_4
          }
        }
        else {
        }
        if (present(cs2_point)) {
          cs2_point=cs2__4
        }
      }
      dlnrhodz = AC_gamma__mod__equationofstate * AC_gravz__mod__gravity/cs2_point
      LNRHO[vertexIdx.x][vertexIdx.y][n1-1-1] = LNRHO[vertexIdx.x][vertexIdx.y][1+n1-1] - AC_dz2_bound__mod__cdata[-1+NGHOST+1-1]*dlnrhodz
      LNRHO[vertexIdx.x][vertexIdx.y][n1-2-1] = LNRHO[vertexIdx.x][vertexIdx.y][2+n1-1] - AC_dz2_bound__mod__cdata[-2+NGHOST+1-1]*dlnrhodz
      LNRHO[vertexIdx.x][vertexIdx.y][n1-3-1] = LNRHO[vertexIdx.x][vertexIdx.y][3+n1-1] - AC_dz2_bound__mod__cdata[-3+NGHOST+1-1]*dlnrhodz
    }
    else {
      potx_xpoint_5=0.0
      poty_ypoint_5=0.0
      potz_zpoint_5=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_5=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[n1-1-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[n1-1-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[n1-1-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[n1-1-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_5=-AC_gravz__mod__gravity*(AC_z__mod__cdata[n1-1-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_5=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[n1-1-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_5=0.5*((AC_z__mod__cdata[n1-1-1]*AC_z__mod__cdata[n1-1-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_5=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[n1-1-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[n1-1-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_5=0.5*((AC_z__mod__cdata[n1-1-1]*AC_z__mod__cdata[n1-1-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_5
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_5 = 1. +pow( (AC_z__mod__cdata[n1-1-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_5 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[n1-1-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[n1-1-1]))/pow(prof_5,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_5=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[n1-1-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potm = potx_xpoint_5 + poty_ypoint_5 + potz_zpoint_5
      potx_xpoint_6=0.0
      poty_ypoint_6=0.0
      potz_zpoint_6=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_6=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[1+n1-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[1+n1-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[1+n1-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[1+n1-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_6=-AC_gravz__mod__gravity*(AC_z__mod__cdata[1+n1-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_6=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[1+n1-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_6=0.5*((AC_z__mod__cdata[1+n1-1]*AC_z__mod__cdata[1+n1-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_6=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[1+n1-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[1+n1-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_6=0.5*((AC_z__mod__cdata[1+n1-1]*AC_z__mod__cdata[1+n1-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_6
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_6 = 1. +pow( (AC_z__mod__cdata[1+n1-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_6 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[1+n1-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[1+n1-1]))/pow(prof_6,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_6=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[1+n1-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potp = potx_xpoint_6 + poty_ypoint_6 + potz_zpoint_6
      cs2 = AC_cs2bot__mod__equationofstate
      if (AC_ldensity_nolog__mod__cdata) {
        RHO[vertexIdx.x][vertexIdx.y][n1-1-1]   = RHO[vertexIdx.x][vertexIdx.y][1+n1-1]*exp(-(potm-potp)/cs2)
      }
      else {
        LNRHO[vertexIdx.x][vertexIdx.y][n1-1-1] = LNRHO[vertexIdx.x][vertexIdx.y][1+n1-1] - (potm-potp)/cs2
      }
      potx_xpoint_5=0.0
      poty_ypoint_5=0.0
      potz_zpoint_5=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_5=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[n1-2-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[n1-2-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[n1-2-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[n1-2-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_5=-AC_gravz__mod__gravity*(AC_z__mod__cdata[n1-2-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_5=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[n1-2-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_5=0.5*((AC_z__mod__cdata[n1-2-1]*AC_z__mod__cdata[n1-2-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_5=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[n1-2-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[n1-2-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_5=0.5*((AC_z__mod__cdata[n1-2-1]*AC_z__mod__cdata[n1-2-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_5
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_5 = 1. +pow( (AC_z__mod__cdata[n1-2-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_5 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[n1-2-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[n1-2-1]))/pow(prof_5,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_5=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[n1-2-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potm = potx_xpoint_5 + poty_ypoint_5 + potz_zpoint_5
      potx_xpoint_6=0.0
      poty_ypoint_6=0.0
      potz_zpoint_6=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_6=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[2+n1-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[2+n1-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[2+n1-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[2+n1-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_6=-AC_gravz__mod__gravity*(AC_z__mod__cdata[2+n1-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_6=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[2+n1-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_6=0.5*((AC_z__mod__cdata[2+n1-1]*AC_z__mod__cdata[2+n1-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_6=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[2+n1-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[2+n1-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_6=0.5*((AC_z__mod__cdata[2+n1-1]*AC_z__mod__cdata[2+n1-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_6
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_6 = 1. +pow( (AC_z__mod__cdata[2+n1-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_6 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[2+n1-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[2+n1-1]))/pow(prof_6,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_6=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[2+n1-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potp = potx_xpoint_6 + poty_ypoint_6 + potz_zpoint_6
      cs2 = AC_cs2bot__mod__equationofstate
      if (AC_ldensity_nolog__mod__cdata) {
        RHO[vertexIdx.x][vertexIdx.y][n1-2-1]   = RHO[vertexIdx.x][vertexIdx.y][2+n1-1]*exp(-(potm-potp)/cs2)
      }
      else {
        LNRHO[vertexIdx.x][vertexIdx.y][n1-2-1] = LNRHO[vertexIdx.x][vertexIdx.y][2+n1-1] - (potm-potp)/cs2
      }
      potx_xpoint_5=0.0
      poty_ypoint_5=0.0
      potz_zpoint_5=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_5=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[n1-3-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[n1-3-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[n1-3-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[n1-3-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_5=-AC_gravz__mod__gravity*(AC_z__mod__cdata[n1-3-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_5=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[n1-3-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_5=0.5*((AC_z__mod__cdata[n1-3-1]*AC_z__mod__cdata[n1-3-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_5=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[n1-3-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[n1-3-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_5=0.5*((AC_z__mod__cdata[n1-3-1]*AC_z__mod__cdata[n1-3-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_5
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_5 = 1. +pow( (AC_z__mod__cdata[n1-3-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_5 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[n1-3-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[n1-3-1]))/pow(prof_5,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_5=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[n1-3-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potm = potx_xpoint_5 + poty_ypoint_5 + potz_zpoint_5
      potx_xpoint_6=0.0
      poty_ypoint_6=0.0
      potz_zpoint_6=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_6=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[3+n1-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[3+n1-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[3+n1-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[3+n1-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_6=-AC_gravz__mod__gravity*(AC_z__mod__cdata[3+n1-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_6=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[3+n1-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_6=0.5*((AC_z__mod__cdata[3+n1-1]*AC_z__mod__cdata[3+n1-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_6=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[3+n1-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[3+n1-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_6=0.5*((AC_z__mod__cdata[3+n1-1]*AC_z__mod__cdata[3+n1-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_6
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_6 = 1. +pow( (AC_z__mod__cdata[3+n1-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_6 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[3+n1-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[3+n1-1]))/pow(prof_6,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_6=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[3+n1-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potp = potx_xpoint_6 + poty_ypoint_6 + potz_zpoint_6
      cs2 = AC_cs2bot__mod__equationofstate
      if (AC_ldensity_nolog__mod__cdata) {
        RHO[vertexIdx.x][vertexIdx.y][n1-3-1]   = RHO[vertexIdx.x][vertexIdx.y][3+n1-1]*exp(-(potm-potp)/cs2)
      }
      else {
        LNRHO[vertexIdx.x][vertexIdx.y][n1-3-1] = LNRHO[vertexIdx.x][vertexIdx.y][3+n1-1] - (potm-potp)/cs2
      }
    }
  }

  else if (topbot == AC_top) {
    if (lentropy) {
      if (ldensity) {
        if (AC_ldensity_nolog__mod__cdata) {
          if (AC_lreference_state__mod__cdata) {
            getrho_s_return_value_7=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]+AC_reference_state__mod__density[nx-1][iref_rho-1]
          }
          else {
            getrho_s_return_value_7=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          }
        }
        else {
          getrho_s_return_value_7=exp(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
        }
        rho=getrho_s_return_value_7
        ss=SS[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
        if (AC_lreference_state__mod__cdata) {
          ss=ss+AC_reference_state__mod__density[AC_xtop__mod__equationofstate-1][iref_s-1]
        }
        if (AC_irho_ss__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate || AC_irho_ss__mod__equationofstate == AC_irho_ss__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else if (AC_irho_ss__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate)  {
            lnrho__8 = rho
          }
          else {
            lnrho__8 = log(rho)
          }
          ss__8=ss
          lntt__8=AC_lntt0__mod__equationofstate+AC_cv1__mod__equationofstate*ss__8+AC_gamma_m1__mod__equationofstate*(lnrho__8-AC_lnrho0__mod__equationofstate)
          ee__8=AC_cv__mod__equationofstate*exp(lntt__8)
          pp__8=(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*exp(lntt__8+lnrho__8)
          cs2__8=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__8
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ilnrho_ee__mod__equationofstate || AC_irho_ss__mod__equationofstate == AC_irho_ee__mod__equationofstate) {
          ee__8=ss
          if (AC_irho_ss__mod__equationofstate==AC_ilnrho_ee__mod__equationofstate) {
            lnrho__8=rho
            pp__8=AC_gamma_m1__mod__equationofstate*ee__8*exp(rho)
          }
          else {
            lnrho__8=log(rho)
            pp__8=AC_gamma_m1__mod__equationofstate*ee__8*rho
          }
          lntt__8=log(AC_cv1__mod__equationofstate*ee__8)
          ss__8=AC_cv__mod__equationofstate*(lntt__8-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__8-AC_lnrho0__mod__equationofstate))
          cs2__8=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__8
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ilnrho_pp__mod__equationofstate || AC_irho_ss__mod__equationofstate == AC_irho_pp__mod__equationofstate) {
          pp__8=ss
          if (AC_irho_ss__mod__equationofstate==AC_ilnrho_pp__mod__equationofstate) {
            lnrho__8=rho
            ee__8=pp__8*exp(-lnrho__8)/AC_gamma_m1__mod__equationofstate
            ss__8=AC_cv__mod__equationofstate*(log(pp__8*exp(-lnrho__8)*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__8-AC_lnrho0__mod__equationofstate))
          }
          else {
            lnrho__8=log(rho)
            ee__8=pp__8/rho/AC_gamma_m1__mod__equationofstate
            ss__8=AC_cv__mod__equationofstate*(log(pp__8/rho*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__8-AC_lnrho0__mod__equationofstate))
          }
          lntt__8=log(AC_cv1__mod__equationofstate*ee__8)
          cs2__8=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__8
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ilnrho_lntt__mod__equationofstate) {
          lnrho__8=rho
          lntt__8=ss
          ss__8=AC_cv__mod__equationofstate*(lntt__8-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__8-AC_lnrho0__mod__equationofstate))
          ee__8=AC_cv__mod__equationofstate*exp(lntt__8)
          pp__8=ee__8*exp(lnrho__8)*AC_gamma_m1__mod__equationofstate
          cs2__8=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__8
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ilnrho_tt__mod__equationofstate) {
          lnrho__8=rho
          tt__8=ss
          ss__8=AC_cv__mod__equationofstate*(log(tt__8)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__8-AC_lnrho0__mod__equationofstate))
          ee__8=AC_cv__mod__equationofstate*tt__8
          pp__8=ee__8*exp(lnrho__8)*AC_gamma_m1__mod__equationofstate
          cs2__8=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__8
        }
        else if (AC_irho_ss__mod__equationofstate == AC_irho_tt__mod__equationofstate) {
          lnrho__8=log(rho)
          tt__8=ss
          ss__8=AC_cv__mod__equationofstate*(log(tt__8)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__8-AC_lnrho0__mod__equationofstate))
          ee__8=AC_cv__mod__equationofstate*tt__8
          pp__8=ee__8*rho*AC_gamma_m1__mod__equationofstate
          cs2__8=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__8
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ipp_cs2__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__8=log(rho)
              tt__8=exp(AC_lntt0__mod__equationofstate)
              pp__8=exp(lnrho__8)*AC_cs20__mod__equationofstate/AC_gamma__mod__equationofstate
            }
            else {
              if (AC_leos_isothermal__mod__equationofstate) {
                pp__8=rho
                lnrho__8=log(pp__8*AC_cs20__mod__equationofstate)
                tt__8=exp(AC_lntt0__mod__equationofstate)
              }
            }
          }
        }
        else if (AC_irho_ss__mod__equationofstate == AC_ipp_ss__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__8=(rho)
              ss__8=ss
              cs2__8=exp(AC_gamma__mod__equationofstate*ss__8*AC_cp1__mod__equationofstate+AC_gamma_m1__mod__equationofstate*(lnrho__8-AC_lnrho0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__8=cs2__8/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
            }
            else {
              pp__8=rho
              ss__8=ss
              cs2__8=exp(ss__8*AC_cp1__mod__equationofstate+AC_gamma1__mod__equationofstate*AC_gamma_m1__mod__equationofstate*log(pp__8/AC_pp0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__8=cs2__8/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
              lnrho__8=log(AC_gamma__mod__equationofstate*pp__8/cs2__8)
            }
          }
        }
        else if (AC_irho_ss__mod__equationofstate == AC_irho_eth__mod__equationofstate || AC_irho_ss__mod__equationofstate == AC_ilnrho_eth__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else {
            if (AC_ldensity_nolog__mod__cdata) {
              rho_8 = rho
              if (false) {
                lnrho__8 = log(rho)
              }
            }
            else {
              rho_8 = exp(rho)
              if (false) {
                lnrho__8 = rho
              }
            }
            eth_8 = ss
          }
          if (false) {
            lntt__8 = log(AC_cv1__mod__equationofstate * eth_8/ rho_8)
          }
          if (false) {
            ee__8 = eth_8/ rho_8
          }
          if (false) {
            pp__8 = AC_gamma_m1__mod__equationofstate * eth_8
          }
          if (present(cs2_point)) {
            cs2__8 = AC_gamma__mod__equationofstate * AC_gamma_m1__mod__equationofstate * eth_8/ rho_8
          }
        }
        else {
        }
        if (present(cs2_point)) {
          cs2_point=cs2__8
        }
        dlnrhodz = AC_gamma__mod__equationofstate *AC_gravz__mod__gravity/cs2_point
        if (AC_ldensity_nolog__mod__cdata) {
          dlnrhodz=dlnrhodz*rho
        }
        dssdz    = -AC_gamma_m1__mod__equationofstate*AC_gravz__mod__gravity/cs2_point
        LNRHO[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1] = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1] + AC_dz2_bound__mod__cdata[2+NGHOST-1]*dlnrhodz
        SS[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1] = SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1] + AC_dz2_bound__mod__cdata[2+NGHOST-1]*dssdz
        LNRHO[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1] = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1] + AC_dz2_bound__mod__cdata[3+NGHOST-1]*dlnrhodz
        SS[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1] = SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1] + AC_dz2_bound__mod__cdata[3+NGHOST-1]*dssdz
        LNRHO[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1] = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1] + AC_dz2_bound__mod__cdata[4+NGHOST-1]*dlnrhodz
        SS[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1] = SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1] + AC_dz2_bound__mod__cdata[4+NGHOST-1]*dssdz
      }
      else {
      }
    }
    else if (ltemperature) {
      if (AC_ltemperature_nolog__mod__cdata) {
        if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate || AC_ilnrho_tt__mod__equationofstate == AC_irho_ss__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate)  {
            lnrho__9 = LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          }
          else {
            lnrho__9 = log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
          }
          ss__9=TT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          lntt__9=AC_lntt0__mod__equationofstate+AC_cv1__mod__equationofstate*ss__9+AC_gamma_m1__mod__equationofstate*(lnrho__9-AC_lnrho0__mod__equationofstate)
          ee__9=AC_cv__mod__equationofstate*exp(lntt__9)
          pp__9=(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*exp(lntt__9+lnrho__9)
          cs2__9=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__9
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_ee__mod__equationofstate || AC_ilnrho_tt__mod__equationofstate == AC_irho_ee__mod__equationofstate) {
          ee__9=TT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          if (AC_ilnrho_tt__mod__equationofstate==AC_ilnrho_ee__mod__equationofstate) {
            lnrho__9=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
            pp__9=AC_gamma_m1__mod__equationofstate*ee__9*exp(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
          }
          else {
            lnrho__9=log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
            pp__9=AC_gamma_m1__mod__equationofstate*ee__9*LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          }
          lntt__9=log(AC_cv1__mod__equationofstate*ee__9)
          ss__9=AC_cv__mod__equationofstate*(lntt__9-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__9-AC_lnrho0__mod__equationofstate))
          cs2__9=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__9
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_pp__mod__equationofstate || AC_ilnrho_tt__mod__equationofstate == AC_irho_pp__mod__equationofstate) {
          pp__9=TT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          if (AC_ilnrho_tt__mod__equationofstate==AC_ilnrho_pp__mod__equationofstate) {
            lnrho__9=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
            ee__9=pp__9*exp(-lnrho__9)/AC_gamma_m1__mod__equationofstate
            ss__9=AC_cv__mod__equationofstate*(log(pp__9*exp(-lnrho__9)*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__9-AC_lnrho0__mod__equationofstate))
          }
          else {
            lnrho__9=log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
            ee__9=pp__9/LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]/AC_gamma_m1__mod__equationofstate
            ss__9=AC_cv__mod__equationofstate*(log(pp__9/LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__9-AC_lnrho0__mod__equationofstate))
          }
          lntt__9=log(AC_cv1__mod__equationofstate*ee__9)
          cs2__9=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__9
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_lntt__mod__equationofstate) {
          lnrho__9=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          lntt__9=TT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          ss__9=AC_cv__mod__equationofstate*(lntt__9-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__9-AC_lnrho0__mod__equationofstate))
          ee__9=AC_cv__mod__equationofstate*exp(lntt__9)
          pp__9=ee__9*exp(lnrho__9)*AC_gamma_m1__mod__equationofstate
          cs2__9=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__9
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_tt__mod__equationofstate) {
          lnrho__9=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          tt__9=TT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          ss__9=AC_cv__mod__equationofstate*(log(tt__9)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__9-AC_lnrho0__mod__equationofstate))
          ee__9=AC_cv__mod__equationofstate*tt__9
          pp__9=ee__9*exp(lnrho__9)*AC_gamma_m1__mod__equationofstate
          cs2__9=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__9
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_irho_tt__mod__equationofstate) {
          lnrho__9=log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
          tt__9=TT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          ss__9=AC_cv__mod__equationofstate*(log(tt__9)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__9-AC_lnrho0__mod__equationofstate))
          ee__9=AC_cv__mod__equationofstate*tt__9
          pp__9=ee__9*LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]*AC_gamma_m1__mod__equationofstate
          cs2__9=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__9
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ipp_cs2__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__9=log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
              tt__9=exp(AC_lntt0__mod__equationofstate)
              pp__9=exp(lnrho__9)*AC_cs20__mod__equationofstate/AC_gamma__mod__equationofstate
            }
            else {
              if (AC_leos_isothermal__mod__equationofstate) {
                pp__9=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
                lnrho__9=log(pp__9*AC_cs20__mod__equationofstate)
                tt__9=exp(AC_lntt0__mod__equationofstate)
              }
            }
          }
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_ipp_ss__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__9=(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
              ss__9=TT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
              cs2__9=exp(AC_gamma__mod__equationofstate*ss__9*AC_cp1__mod__equationofstate+AC_gamma_m1__mod__equationofstate*(lnrho__9-AC_lnrho0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__9=cs2__9/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
            }
            else {
              pp__9=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
              ss__9=TT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
              cs2__9=exp(ss__9*AC_cp1__mod__equationofstate+AC_gamma1__mod__equationofstate*AC_gamma_m1__mod__equationofstate*log(pp__9/AC_pp0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__9=cs2__9/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
              lnrho__9=log(AC_gamma__mod__equationofstate*pp__9/cs2__9)
            }
          }
        }
        else if (AC_ilnrho_tt__mod__equationofstate == AC_irho_eth__mod__equationofstate || AC_ilnrho_tt__mod__equationofstate == AC_ilnrho_eth__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else {
            if (AC_ldensity_nolog__mod__cdata) {
              rho_9 = LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
              if (false) {
                lnrho__9 = log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
              }
            }
            else {
              rho_9 = exp(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
              if (false) {
                lnrho__9 = LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
              }
            }
            eth_9 = TT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          }
          if (false) {
            lntt__9 = log(AC_cv1__mod__equationofstate * eth_9/ rho_9)
          }
          if (false) {
            ee__9 = eth_9/ rho_9
          }
          if (false) {
            pp__9 = AC_gamma_m1__mod__equationofstate * eth_9
          }
          if (present(cs2_point)) {
            cs2__9 = AC_gamma__mod__equationofstate * AC_gamma_m1__mod__equationofstate * eth_9/ rho_9
          }
        }
        else {
        }
        if (present(cs2_point)) {
          cs2_point=cs2__9
        }
      }
      else {
        if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate || AC_ilnrho_lntt__mod__equationofstate == AC_irho_ss__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_ss__mod__equationofstate)  {
            lnrho__10 = LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          }
          else {
            lnrho__10 = log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
          }
          ss__10=LNTT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          lntt__10=AC_lntt0__mod__equationofstate+AC_cv1__mod__equationofstate*ss__10+AC_gamma_m1__mod__equationofstate*(lnrho__10-AC_lnrho0__mod__equationofstate)
          ee__10=AC_cv__mod__equationofstate*exp(lntt__10)
          pp__10=(AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*exp(lntt__10+lnrho__10)
          cs2__10=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__10
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_ee__mod__equationofstate || AC_ilnrho_lntt__mod__equationofstate == AC_irho_ee__mod__equationofstate) {
          ee__10=LNTT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          if (AC_ilnrho_lntt__mod__equationofstate==AC_ilnrho_ee__mod__equationofstate) {
            lnrho__10=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
            pp__10=AC_gamma_m1__mod__equationofstate*ee__10*exp(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
          }
          else {
            lnrho__10=log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
            pp__10=AC_gamma_m1__mod__equationofstate*ee__10*LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          }
          lntt__10=log(AC_cv1__mod__equationofstate*ee__10)
          ss__10=AC_cv__mod__equationofstate*(lntt__10-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__10-AC_lnrho0__mod__equationofstate))
          cs2__10=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__10
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_pp__mod__equationofstate || AC_ilnrho_lntt__mod__equationofstate == AC_irho_pp__mod__equationofstate) {
          pp__10=LNTT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          if (AC_ilnrho_lntt__mod__equationofstate==AC_ilnrho_pp__mod__equationofstate) {
            lnrho__10=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
            ee__10=pp__10*exp(-lnrho__10)/AC_gamma_m1__mod__equationofstate
            ss__10=AC_cv__mod__equationofstate*(log(pp__10*exp(-lnrho__10)*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__10-AC_lnrho0__mod__equationofstate))
          }
          else {
            lnrho__10=log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
            ee__10=pp__10/LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]/AC_gamma_m1__mod__equationofstate
            ss__10=AC_cv__mod__equationofstate*(log(pp__10/LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]*AC_gamma__mod__equationofstate/AC_cs20__mod__equationofstate)-AC_gamma_m1__mod__equationofstate*(lnrho__10-AC_lnrho0__mod__equationofstate))
          }
          lntt__10=log(AC_cv1__mod__equationofstate*ee__10)
          cs2__10=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__10
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_lntt__mod__equationofstate) {
          lnrho__10=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          lntt__10=LNTT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          ss__10=AC_cv__mod__equationofstate*(lntt__10-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__10-AC_lnrho0__mod__equationofstate))
          ee__10=AC_cv__mod__equationofstate*exp(lntt__10)
          pp__10=ee__10*exp(lnrho__10)*AC_gamma_m1__mod__equationofstate
          cs2__10=AC_gamma__mod__equationofstate*AC_gamma_m1__mod__equationofstate*ee__10
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_tt__mod__equationofstate) {
          lnrho__10=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          tt__10=LNTT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          ss__10=AC_cv__mod__equationofstate*(log(tt__10)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__10-AC_lnrho0__mod__equationofstate))
          ee__10=AC_cv__mod__equationofstate*tt__10
          pp__10=ee__10*exp(lnrho__10)*AC_gamma_m1__mod__equationofstate
          cs2__10=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__10
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_irho_tt__mod__equationofstate) {
          lnrho__10=log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
          tt__10=LNTT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          ss__10=AC_cv__mod__equationofstate*(log(tt__10)-AC_lntt0__mod__equationofstate-AC_gamma_m1__mod__equationofstate*(lnrho__10-AC_lnrho0__mod__equationofstate))
          ee__10=AC_cv__mod__equationofstate*tt__10
          pp__10=ee__10*LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]*AC_gamma_m1__mod__equationofstate
          cs2__10=AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate*tt__10
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ipp_cs2__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__10=log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
              tt__10=exp(AC_lntt0__mod__equationofstate)
              pp__10=exp(lnrho__10)*AC_cs20__mod__equationofstate/AC_gamma__mod__equationofstate
            }
            else {
              if (AC_leos_isothermal__mod__equationofstate) {
                pp__10=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
                lnrho__10=log(pp__10*AC_cs20__mod__equationofstate)
                tt__10=exp(AC_lntt0__mod__equationofstate)
              }
            }
          }
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_ipp_ss__mod__equationofstate) {
          if (lanelastic) {
            if (AC_lanelastic_lin__mod__equationofstate) {
              lnrho__10=(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
              ss__10=LNTT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
              cs2__10=exp(AC_gamma__mod__equationofstate*ss__10*AC_cp1__mod__equationofstate+AC_gamma_m1__mod__equationofstate*(lnrho__10-AC_lnrho0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__10=cs2__10/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
            }
            else {
              pp__10=LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
              ss__10=LNTT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
              cs2__10=exp(ss__10*AC_cp1__mod__equationofstate+AC_gamma1__mod__equationofstate*AC_gamma_m1__mod__equationofstate*log(pp__10/AC_pp0__mod__equationofstate))*AC_cs20__mod__equationofstate
              tt__10=cs2__10/(AC_gamma_m1__mod__equationofstate*AC_cp__mod__equationofstate)
              lnrho__10=log(AC_gamma__mod__equationofstate*pp__10/cs2__10)
            }
          }
        }
        else if (AC_ilnrho_lntt__mod__equationofstate == AC_irho_eth__mod__equationofstate || AC_ilnrho_lntt__mod__equationofstate == AC_ilnrho_eth__mod__equationofstate) {
          if (AC_lstratz__mod__cdata) {
          }
          else {
            if (AC_ldensity_nolog__mod__cdata) {
              rho_10 = LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
              if (false) {
                lnrho__10 = log(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
              }
            }
            else {
              rho_10 = exp(LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1])
              if (false) {
                lnrho__10 = LNRHO[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
              }
            }
            eth_10 = LNTT[AC_l2__mod__cdata-1][AC_m2__mod__cdata-1][AC_n2__mod__cdata-1]
          }
          if (false) {
            lntt__10 = log(AC_cv1__mod__equationofstate * eth_10/ rho_10)
          }
          if (false) {
            ee__10 = eth_10/ rho_10
          }
          if (false) {
            pp__10 = AC_gamma_m1__mod__equationofstate * eth_10
          }
          if (present(cs2_point)) {
            cs2__10 = AC_gamma__mod__equationofstate * AC_gamma_m1__mod__equationofstate * eth_10/ rho_10
          }
        }
        else {
        }
        if (present(cs2_point)) {
          cs2_point=cs2__10
        }
      }
      dlnrhodz =  AC_gamma__mod__equationofstate *AC_gravz__mod__gravity/cs2_point
      LNRHO[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1] = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1] + AC_dz2_bound__mod__cdata[2+NGHOST-1]*dlnrhodz
      LNRHO[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1] = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1] + AC_dz2_bound__mod__cdata[3+NGHOST-1]*dlnrhodz
      LNRHO[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1] = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1] + AC_dz2_bound__mod__cdata[4+NGHOST-1]*dlnrhodz
    }
    else {
      potx_xpoint_11=0.0
      poty_ypoint_11=0.0
      potz_zpoint_11=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_11=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[1+AC_n2__mod__cdata-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[1+AC_n2__mod__cdata-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[1+AC_n2__mod__cdata-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[1+AC_n2__mod__cdata-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_11=-AC_gravz__mod__gravity*(AC_z__mod__cdata[1+AC_n2__mod__cdata-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_11=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[1+AC_n2__mod__cdata-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_11=0.5*((AC_z__mod__cdata[1+AC_n2__mod__cdata-1]*AC_z__mod__cdata[1+AC_n2__mod__cdata-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_11=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[1+AC_n2__mod__cdata-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[1+AC_n2__mod__cdata-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_11=0.5*((AC_z__mod__cdata[1+AC_n2__mod__cdata-1]*AC_z__mod__cdata[1+AC_n2__mod__cdata-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_11
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_11 = 1. +pow( (AC_z__mod__cdata[1+AC_n2__mod__cdata-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_11 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[1+AC_n2__mod__cdata-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[1+AC_n2__mod__cdata-1]))/pow(prof_11,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_11=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[1+AC_n2__mod__cdata-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potp = potx_xpoint_11 + poty_ypoint_11 + potz_zpoint_11
      potx_xpoint_12=0.0
      poty_ypoint_12=0.0
      potz_zpoint_12=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_12=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-1-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-1-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-1-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[AC_n2__mod__cdata-1-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_12=-AC_gravz__mod__gravity*(AC_z__mod__cdata[AC_n2__mod__cdata-1-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_12=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[AC_n2__mod__cdata-1-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_12=0.5*((AC_z__mod__cdata[AC_n2__mod__cdata-1-1]*AC_z__mod__cdata[AC_n2__mod__cdata-1-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_12=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[AC_n2__mod__cdata-1-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[AC_n2__mod__cdata-1-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_12=0.5*((AC_z__mod__cdata[AC_n2__mod__cdata-1-1]*AC_z__mod__cdata[AC_n2__mod__cdata-1-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_12
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_12 = 1. +pow( (AC_z__mod__cdata[AC_n2__mod__cdata-1-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_12 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-1-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-1-1]))/pow(prof_12,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_12=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[AC_n2__mod__cdata-1-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potm = potx_xpoint_12 + poty_ypoint_12 + potz_zpoint_12
      cs2 = AC_cs2bot__mod__equationofstate
      if (AC_ldensity_nolog__mod__cdata) {
        RHO[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1]   = RHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1]*exp(-(potp-potm)/cs2)
      }
      else {
        LNRHO[vertexIdx.x][vertexIdx.y][1+AC_n2__mod__cdata-1] = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1-1] - (potp-potm)/cs2
      }
      potx_xpoint_11=0.0
      poty_ypoint_11=0.0
      potz_zpoint_11=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_11=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[2+AC_n2__mod__cdata-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[2+AC_n2__mod__cdata-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[2+AC_n2__mod__cdata-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[2+AC_n2__mod__cdata-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_11=-AC_gravz__mod__gravity*(AC_z__mod__cdata[2+AC_n2__mod__cdata-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_11=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[2+AC_n2__mod__cdata-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_11=0.5*((AC_z__mod__cdata[2+AC_n2__mod__cdata-1]*AC_z__mod__cdata[2+AC_n2__mod__cdata-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_11=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[2+AC_n2__mod__cdata-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[2+AC_n2__mod__cdata-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_11=0.5*((AC_z__mod__cdata[2+AC_n2__mod__cdata-1]*AC_z__mod__cdata[2+AC_n2__mod__cdata-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_11
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_11 = 1. +pow( (AC_z__mod__cdata[2+AC_n2__mod__cdata-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_11 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[2+AC_n2__mod__cdata-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[2+AC_n2__mod__cdata-1]))/pow(prof_11,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_11=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[2+AC_n2__mod__cdata-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potp = potx_xpoint_11 + poty_ypoint_11 + potz_zpoint_11
      potx_xpoint_12=0.0
      poty_ypoint_12=0.0
      potz_zpoint_12=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_12=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-2-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-2-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-2-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[AC_n2__mod__cdata-2-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_12=-AC_gravz__mod__gravity*(AC_z__mod__cdata[AC_n2__mod__cdata-2-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_12=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[AC_n2__mod__cdata-2-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_12=0.5*((AC_z__mod__cdata[AC_n2__mod__cdata-2-1]*AC_z__mod__cdata[AC_n2__mod__cdata-2-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_12=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[AC_n2__mod__cdata-2-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[AC_n2__mod__cdata-2-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_12=0.5*((AC_z__mod__cdata[AC_n2__mod__cdata-2-1]*AC_z__mod__cdata[AC_n2__mod__cdata-2-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_12
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_12 = 1. +pow( (AC_z__mod__cdata[AC_n2__mod__cdata-2-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_12 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-2-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-2-1]))/pow(prof_12,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_12=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[AC_n2__mod__cdata-2-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potm = potx_xpoint_12 + poty_ypoint_12 + potz_zpoint_12
      cs2 = AC_cs2bot__mod__equationofstate
      if (AC_ldensity_nolog__mod__cdata) {
        RHO[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1]   = RHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1]*exp(-(potp-potm)/cs2)
      }
      else {
        LNRHO[vertexIdx.x][vertexIdx.y][2+AC_n2__mod__cdata-1] = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-2-1] - (potp-potm)/cs2
      }
      potx_xpoint_11=0.0
      poty_ypoint_11=0.0
      potz_zpoint_11=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_11=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[3+AC_n2__mod__cdata-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[3+AC_n2__mod__cdata-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[3+AC_n2__mod__cdata-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[3+AC_n2__mod__cdata-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_11=-AC_gravz__mod__gravity*(AC_z__mod__cdata[3+AC_n2__mod__cdata-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_11=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[3+AC_n2__mod__cdata-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_11=0.5*((AC_z__mod__cdata[3+AC_n2__mod__cdata-1]*AC_z__mod__cdata[3+AC_n2__mod__cdata-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_11=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[3+AC_n2__mod__cdata-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[3+AC_n2__mod__cdata-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_11=0.5*((AC_z__mod__cdata[3+AC_n2__mod__cdata-1]*AC_z__mod__cdata[3+AC_n2__mod__cdata-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_11
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_11 = 1. +pow( (AC_z__mod__cdata[3+AC_n2__mod__cdata-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_11 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[3+AC_n2__mod__cdata-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[3+AC_n2__mod__cdata-1]))/pow(prof_11,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_11=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[3+AC_n2__mod__cdata-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potp = potx_xpoint_11 + poty_ypoint_11 + potz_zpoint_11
      potx_xpoint_12=0.0
      poty_ypoint_12=0.0
      potz_zpoint_12=0.0
      if (false) {
        if (AC_enum_gravx_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_const_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_linear_zdep_string) {
          zdep_12=(1.+AC_kappa_z1__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-3-1]+0.5*((AC_kappa_z1__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-3-1])*(AC_kappa_z1__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-3-1])))
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_kepler_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czbot1_string) {
        }
        else if (AC_enum_gravx_profile__mod__gravity == enum_czmid1_string) {
        }
        else {
        }
      }
      if (false) {
        if (AC_enum_gravy_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravy_profile__mod__gravity == enum_const_string) {
        }
        else {
        }
      }
      if (present(AC_z__mod__cdata[AC_n2__mod__cdata-3-1])) {
        if (AC_enum_gravz_profile__mod__gravity == enum_zero_string) {
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_const_string) {
          if (AC_zclip__mod__gravity==impossible) {
            potz_zpoint_12=-AC_gravz__mod__gravity*(AC_z__mod__cdata[AC_n2__mod__cdata-3-1]-AC_zinfty__mod__gravity)
          }
          else {
            potz_zpoint_12=-AC_gravz__mod__gravity*max(AC_z__mod__cdata[AC_n2__mod__cdata-3-1]-AC_zinfty__mod__gravity,AC_zclip__mod__gravity-AC_zinfty__mod__gravity)
          }
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_string) {
          potz_zpoint_12=0.5*((AC_z__mod__cdata[AC_n2__mod__cdata-3-1]*AC_z__mod__cdata[AC_n2__mod__cdata-3-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_spherical_string) {
          potz_zpoint_12=0.5*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*(AC_z1__mod__gravity*AC_z1__mod__gravity)*log(1.0+((AC_z__mod__cdata[AC_n2__mod__cdata-3-1]/AC_z1__mod__gravity)*(AC_z__mod__cdata[AC_n2__mod__cdata-3-1]/AC_z1__mod__gravity)))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_xdep_string) {
          potz_zpoint_12=0.5*((AC_z__mod__cdata[AC_n2__mod__cdata-3-1]*AC_z__mod__cdata[AC_n2__mod__cdata-3-1])-(AC_zinfty__mod__gravity*AC_zinfty__mod__gravity))*(AC_nu_epicycle__mod__gravity*AC_nu_epicycle__mod__gravity)*xdep_12
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_linear_smoothed_string) {
          prof_12 = 1. +pow( (AC_z__mod__cdata[AC_n2__mod__cdata-3-1]/AC_zref__mod__gravity),(2*AC_n_pot__mod__gravity))
          potz_zpoint_12 = 0.5*((AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-3-1])*(AC_nu_epicycle__mod__gravity*AC_z__mod__cdata[AC_n2__mod__cdata-3-1]))/pow(prof_12,(1./AC_n_pot__mod__gravity))
        }
        else if (AC_enum_gravz_profile__mod__gravity == enum_tanh_string) {
          potz_zpoint_12=AC_gravz__mod__gravity*AC_zref__mod__gravity*log(cosh(AC_z__mod__cdata[AC_n2__mod__cdata-3-1]/AC_zref__mod__gravity))
        }
        else {
        }
      }
      potm = potx_xpoint_12 + poty_ypoint_12 + potz_zpoint_12
      cs2 = AC_cs2bot__mod__equationofstate
      if (AC_ldensity_nolog__mod__cdata) {
        RHO[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1]   = RHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1]*exp(-(potp-potm)/cs2)
      }
      else {
        LNRHO[vertexIdx.x][vertexIdx.y][3+AC_n2__mod__cdata-1] = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-3-1] - (potp-potm)/cs2
      }
    }
  }
  else {
  }
}
#endif

//These don't actually do anything except set some flags in the CPU code so can be no-ops here
bc_freeze_var_x(AcBoundary boundary, AC_TOP_BOT topbot, Field f)
{
	suppress_unused_warning(boundary)
	suppress_unused_warning(topbot)
	suppress_unused_warning(f)
}

bc_freeze_var_y(AcBoundary boundary, AC_TOP_BOT topbot, Field f)
{
	suppress_unused_warning(boundary)
	suppress_unused_warning(topbot)
	suppress_unused_warning(f)
}

bc_freeze_var_z(AcBoundary boundary, AC_TOP_BOT topbot, Field f)
{
	suppress_unused_warning(boundary)
	suppress_unused_warning(topbot)
	suppress_unused_warning(f)
}

bcx_extrap_2_3(AcBoundary boundary, AC_TOP_BOT topbot,Field f)
{
  suppress_unused_warning(boundary)
  real yl1
  real ypi
  real ymi
  real xl1
  real xmi
  real xpi
  real yyi
  real xl2
  real yl2
  if (topbot == AC_bot) {
    yl1=log(f[AC_l1-1][vertexIdx.y][vertexIdx.z])
    ypi=log(f[1+AC_l1-1][vertexIdx.y][vertexIdx.z])
    xl1=log(AC_x__mod__cdata[AC_l1-1])
    xmi=log(AC_x__mod__cdata[AC_l1-1-1])
    xpi=log(AC_x__mod__cdata[1+AC_l1-1])
    yyi = yl1 - (ypi-yl1)*(xl1-xmi)/(xpi-xl1)
    f[AC_l1-1-1][vertexIdx.y][vertexIdx.z] = exp(yyi)
    yl1=log(f[AC_l1-1][vertexIdx.y][vertexIdx.z])
    ypi=log(f[2+AC_l1-1][vertexIdx.y][vertexIdx.z])
    xl1=log(AC_x__mod__cdata[AC_l1-1])
    xmi=log(AC_x__mod__cdata[AC_l1-2-1])
    xpi=log(AC_x__mod__cdata[2+AC_l1-1])
    yyi = yl1 - (ypi-yl1)*(xl1-xmi)/(xpi-xl1)
    f[AC_l1-2-1][vertexIdx.y][vertexIdx.z] = exp(yyi)
    yl1=log(f[AC_l1-1][vertexIdx.y][vertexIdx.z])
    ypi=log(f[3+AC_l1-1][vertexIdx.y][vertexIdx.z])
    xl1=log(AC_x__mod__cdata[AC_l1-1])
    xmi=log(AC_x__mod__cdata[AC_l1-3-1])
    xpi=log(AC_x__mod__cdata[3+AC_l1-1])
    yyi = yl1 - (ypi-yl1)*(xl1-xmi)/(xpi-xl1)
    f[AC_l1-3-1][vertexIdx.y][vertexIdx.z] = exp(yyi)
  }
  else if (topbot == AC_top) {
    yl2=log(f[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
    ymi=log(f[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z])
    xpi=log(AC_x__mod__cdata[1+AC_l2__mod__cdata-1])
    xl2=log(AC_x__mod__cdata[AC_l2__mod__cdata-1])
    xmi=log(AC_x__mod__cdata[AC_l2__mod__cdata-1-1])
    yyi = yl2 + (yl2-ymi)*(xpi-xl2)/(xl2-xmi)
    f[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = exp(yyi)
    yl2=log(f[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
    ymi=log(f[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z])
    xpi=log(AC_x__mod__cdata[2+AC_l2__mod__cdata-1])
    xl2=log(AC_x__mod__cdata[AC_l2__mod__cdata-1])
    xmi=log(AC_x__mod__cdata[AC_l2__mod__cdata-2-1])
    yyi = yl2 + (yl2-ymi)*(xpi-xl2)/(xl2-xmi)
    f[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = exp(yyi)
    yl2=log(f[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
    ymi=log(f[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z])
    xpi=log(AC_x__mod__cdata[3+AC_l2__mod__cdata-1])
    xl2=log(AC_x__mod__cdata[AC_l2__mod__cdata-1])
    xmi=log(AC_x__mod__cdata[AC_l2__mod__cdata-3-1])
    yyi = yl2 + (yl2-ymi)*(xpi-xl2)/(xl2-xmi)
    f[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = exp(yyi)
  }
  else {
  }
}
bc_outflow_y(AcBoundary boundary, AC_TOP_BOT topbot,Field f,bool lforce_ghost)
{
  suppress_unused_warning(boundary)
  int ix
  int iz
  bool lforce
  lforce = false
  if (present (lforce_ghost)) {
    lforce = lforce_ghost
  }
  if (topbot == AC_bot) {
    if (f[vertexIdx.x][AC_m1-1][vertexIdx.z]<0.0) {
      f[vertexIdx.x][AC_m1-1-1][vertexIdx.z]=+f[vertexIdx.x][1+AC_m1-1][vertexIdx.z]
      f[vertexIdx.x][AC_m1-2-1][vertexIdx.z]=+f[vertexIdx.x][2+AC_m1-1][vertexIdx.z]
      f[vertexIdx.x][AC_m1-3-1][vertexIdx.z]=+f[vertexIdx.x][3+AC_m1-1][vertexIdx.z]
    }
    else {
      f[vertexIdx.x][AC_m1-1-1][vertexIdx.z]=-f[vertexIdx.x][1+AC_m1-1][vertexIdx.z]
      f[vertexIdx.x][AC_m1-2-1][vertexIdx.z]=-f[vertexIdx.x][2+AC_m1-1][vertexIdx.z]
      f[vertexIdx.x][AC_m1-3-1][vertexIdx.z]=-f[vertexIdx.x][3+AC_m1-1][vertexIdx.z]
      f[vertexIdx.x][AC_m1-1][vertexIdx.z]=0.0
    }
    if (lforce) {
      if (f[vertexIdx.x][AC_m1-0-1][vertexIdx.z] > 0.0) {
        f[vertexIdx.x][AC_m1-0-1][vertexIdx.z] = 0.0
      }
      if (f[vertexIdx.x][AC_m1-1-1][vertexIdx.z] > 0.0) {
        f[vertexIdx.x][AC_m1-1-1][vertexIdx.z] = 0.0
      }
      if (f[vertexIdx.x][AC_m1-2-1][vertexIdx.z] > 0.0) {
        f[vertexIdx.x][AC_m1-2-1][vertexIdx.z] = 0.0
      }
      if (f[vertexIdx.x][AC_m1-3-1][vertexIdx.z] > 0.0) {
        f[vertexIdx.x][AC_m1-3-1][vertexIdx.z] = 0.0
      }
    }
  }
  else if (topbot == AC_top) {
    if (f[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z]>0.0) {
      f[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z]=+f[vertexIdx.x][AC_m2__mod__cdata-1-1][vertexIdx.z]
      f[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z]=+f[vertexIdx.x][AC_m2__mod__cdata-2-1][vertexIdx.z]
      f[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z]=+f[vertexIdx.x][AC_m2__mod__cdata-3-1][vertexIdx.z]
    }
    else {
      f[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z]=-f[vertexIdx.x][AC_m2__mod__cdata-1-1][vertexIdx.z]
      f[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z]=-f[vertexIdx.x][AC_m2__mod__cdata-2-1][vertexIdx.z]
      f[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z]=-f[vertexIdx.x][AC_m2__mod__cdata-3-1][vertexIdx.z]
      f[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z]=0.0
    }
    if (lforce) {
      if (f[vertexIdx.x][0+AC_m2__mod__cdata-1][vertexIdx.z] < 0.0) {
        f[vertexIdx.x][0+AC_m2__mod__cdata-1][vertexIdx.z] = 0.0
      }
      if (f[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z] < 0.0) {
        f[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z] = 0.0
      }
      if (f[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z] < 0.0) {
        f[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z] = 0.0
      }
      if (f[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z] < 0.0) {
        f[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z] = 0.0
      }
    }
  }
  else {
  }
}

#if Leos_idealgas_MODULE
bc_stratified_y(AcBoundary boundary, AC_TOP_BOT topbot,Field f)
{
  suppress_unused_warning(boundary)
  real rad=AC_x__mod__cdata[vertexIdx.x]
  if (topbot == AC_bot) {
    za=rad*AC_costh__mod__cdata[m1-1]
    h=AC_cs0__mod__equationofstate*rad/sqrt(AC_gamma__mod__equationofstate)
    zg=rad*AC_costh__mod__cdata[m1-1-1]
    if (AC_ldensity_nolog__mod__cdata) {
      lnrho = log(f[vertexIdx.x][m1-1][vertexIdx.z]) - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][m1-1-1][vertexIdx.z] = exp(lnrho)
    }
    else {
      lnrho = f[vertexIdx.x][m1-1][vertexIdx.z] - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][m1-1-1][vertexIdx.z] = lnrho
    }
    zg=rad*AC_costh__mod__cdata[m1-2-1]
    if (AC_ldensity_nolog__mod__cdata) {
      lnrho = log(f[vertexIdx.x][m1-1][vertexIdx.z]) - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][m1-2-1][vertexIdx.z] = exp(lnrho)
    }
    else {
      lnrho = f[vertexIdx.x][m1-1][vertexIdx.z] - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][m1-2-1][vertexIdx.z] = lnrho
    }
    zg=rad*AC_costh__mod__cdata[m1-3-1]
    if (AC_ldensity_nolog__mod__cdata) {
      lnrho = log(f[vertexIdx.x][m1-1][vertexIdx.z]) - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][m1-3-1][vertexIdx.z] = exp(lnrho)
    }
    else {
      lnrho = f[vertexIdx.x][m1-1][vertexIdx.z] - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][m1-3-1][vertexIdx.z] = lnrho
    }
  }
  else if (topbot == AC_top) {
    za=rad*AC_costh__mod__cdata[AC_m2__mod__cdata-1]
    h=AC_cs0__mod__equationofstate*rad/sqrt(AC_gamma__mod__equationofstate)
    zg=rad*AC_costh__mod__cdata[1+AC_m2__mod__cdata-1]
    if (AC_ldensity_nolog__mod__cdata) {
      lnrho = log(f[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z]) - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z] = exp(lnrho)
    }
    else {
      lnrho = f[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z] - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z] = lnrho
    }
    zg=rad*AC_costh__mod__cdata[2+AC_m2__mod__cdata-1]
    if (AC_ldensity_nolog__mod__cdata) {
      lnrho = log(f[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z]) - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z] = exp(lnrho)
    }
    else {
      lnrho = f[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z] - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z] = lnrho
    }
    zg=rad*AC_costh__mod__cdata[3+AC_m2__mod__cdata-1]
    if (AC_ldensity_nolog__mod__cdata) {
      lnrho = log(f[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z]) - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z] = exp(lnrho)
    }
    else {
      lnrho = f[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z] - ((zg*zg)-(za*za))/(2*(h*h))
      f[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z] = lnrho
    }
  }
  else {
  }
}
#endif

bc_outflow_x(AcBoundary boundary, AC_TOP_BOT topbot,Field f,bool lforce)
{
  suppress_unused_warning(boundary)
  if (topbot == AC_bot) {
    if (f[l1-1][vertexIdx.y][vertexIdx.z]<0.0) {
      f[l1-1-1][vertexIdx.y][vertexIdx.z]=+f[1+l1-1][vertexIdx.y][vertexIdx.z]
      f[l1-2-1][vertexIdx.y][vertexIdx.z]=+f[2+l1-1][vertexIdx.y][vertexIdx.z]
      f[l1-3-1][vertexIdx.y][vertexIdx.z]=+f[3+l1-1][vertexIdx.y][vertexIdx.z]
    }
    else {
      f[l1-1-1][vertexIdx.y][vertexIdx.z]=-f[1+l1-1][vertexIdx.y][vertexIdx.z]
      f[l1-2-1][vertexIdx.y][vertexIdx.z]=-f[2+l1-1][vertexIdx.y][vertexIdx.z]
      f[l1-3-1][vertexIdx.y][vertexIdx.z]=-f[3+l1-1][vertexIdx.y][vertexIdx.z]
      f[l1-1][vertexIdx.y][vertexIdx.z]=0.0
    }
    if (lforce) {
      if (f[l1-0-1][vertexIdx.y][vertexIdx.z] > 0.0) {
        f[l1-0-1][vertexIdx.y][vertexIdx.z] = 0.0
      }
      if (f[l1-1-1][vertexIdx.y][vertexIdx.z] > 0.0) {
        f[l1-1-1][vertexIdx.y][vertexIdx.z] = 0.0
      }
      if (f[l1-2-1][vertexIdx.y][vertexIdx.z] > 0.0) {
        f[l1-2-1][vertexIdx.y][vertexIdx.z] = 0.0
      }
      if (f[l1-3-1][vertexIdx.y][vertexIdx.z] > 0.0) {
        f[l1-3-1][vertexIdx.y][vertexIdx.z] = 0.0
      }
    }
  }
  else if (topbot == AC_top) {
    if (f[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]>0.0) {
      f[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=+f[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]
      f[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=+f[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]
      f[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=+f[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]
    }
    else {
      f[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=-f[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]
      f[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=-f[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]
      f[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=-f[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]
      f[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=0.0
    }
    if (lforce) {
      if (f[0+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] < 0.0) {
        f[0+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = 0.0
      }
      if (f[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] < 0.0) {
        f[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = 0.0
      }
      if (f[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] < 0.0) {
        f[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = 0.0
      }
      if (f[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] < 0.0) {
        f[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z] = 0.0
      }
    }
  }
  else {
  }
}


#if Lhydro_MODULE && Ldensity_MODULE && Leos_idealgas_MODULE
bc_ss_stemp_x(AcBoundary boundary,AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  int id_1
  int id_3
  real rho_yz
  if (topbot == AC_bot) {
    if (AC_lreference_state__mod__cdata) {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          rho_yz=LNRHO[l1-1][vertexIdx.y][vertexIdx.z]+reference_state__mod__density[AC_xbot__mod__equationofstate-1][iref_rho-1]
        }
        else {
          rho_yz=LNRHO[l1-1][vertexIdx.y][vertexIdx.z]
        }
      }
      else {
        rho_yz=exp(LNRHO[l1-1][vertexIdx.y][vertexIdx.z])
      }
    }
    for i in 1:NGHOST+1 {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          dlnrho = LNRHO[l1+i-1][vertexIdx.y][vertexIdx.z]-LNRHO[l1-i-1][vertexIdx.y][vertexIdx.z]
          if (AC_ldensity_nolog__mod__cdata) {
            if (AC_lreference_state__mod__cdata) {
              if (l1 <= (nx+1)/2) {
                id_1 = -i
              }
              else {
                id_1 = i
              }
              dlnrho = dlnrho + dx2_bound__mod__cdata[id_1+1+NGHOST-1]*AC_reference_state__mod__density[AC_l1-1][iref_grho-1]
            }
            dlnrho = dlnrho/rho_yz
          }
          SS[l1-i-1][vertexIdx.y][vertexIdx.z] =  SS[l1+i-1][vertexIdx.y][vertexIdx.z] + dx2_bound__mod__cdata[-i+NGHOST+1-1]*AC_reference_state__mod__density[AC_xbot__mod__equationofstate-1][iref_gs-1]  + (AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*dlnrho
        }
        else {
          SS[l1-i-1][vertexIdx.y][vertexIdx.z] =  SS[l1+i-1][vertexIdx.y][vertexIdx.z] + (cp__mod__equationofstate-AC_cv__mod__equationofstate)*(log(RHO[AC_l1+i-1][vertexIdx.y][vertexIdx.z]/RHO[l1-i-1][vertexIdx.y][vertexIdx.z]))
        }
      }
      else {
        SS[l1-i-1][vertexIdx.y][vertexIdx.z] =  SS[l1+i-1][vertexIdx.y][vertexIdx.z] + (cp__mod__equationofstate-AC_cv__mod__equationofstate)*(LNRHO[AC_l1+i-1][vertexIdx.y][vertexIdx.z]-LNRHO[l1-i-1][vertexIdx.y][vertexIdx.z])
      }
    }
  }
  else if (topbot == AC_top) {
    if (AC_lreference_state__mod__cdata) {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          rho_yz=LNRHO[l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]+AC_reference_state__mod__density[AC_xtop__mod__equationofstate-1][iref_rho-1]
        }
        else {
          rho_yz=LNRHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]
        }
      }
      else {
        rho_yz=exp(LNRHO[AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z])
      }
    }
    for i in 1:NGHOST+1 {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          dlnrho = LNRHO[AC_l2__mod__cdata+i-1][vertexIdx.y][vertexIdx.z]-LNRHO[AC_l2__mod__cdata-i-1][vertexIdx.y][vertexIdx.z]
          if (AC_ldensity_nolog__mod__cdata) {
            if (AC_lreference_state__mod__cdata) {
              if (l2__mod__cdata <= (AC_nx+1)/2) {
                id_3 = -i
              }
              else {
                id_3 = i
              }
              dlnrho = dlnrho + dx2_bound__mod__cdata[id_3+1+NGHOST-1]*AC_reference_state__mod__density[AC_l2__mod__cdata-1][iref_grho-1]
            }
            dlnrho = dlnrho/rho_yz
          }
          SS[l2__mod__cdata+i-1][vertexIdx.y][vertexIdx.z] =  SS[AC_l2__mod__cdata-i-1][vertexIdx.y][vertexIdx.z] - AC_dx2_bound__mod__cdata[i+1+NGHOST-1]*AC_reference_state__mod__density[AC_xtop__mod__equationofstate-1][iref_gs-1]  - (AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*dlnrho
        }
        else {
          SS[AC_l2__mod__cdata+i-1][vertexIdx.y][vertexIdx.z] = SS[AC_l2__mod__cdata-i-1][vertexIdx.y][vertexIdx.z] + (AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*log(RHO[AC_l2__mod__cdata-i-1][vertexIdx.y][vertexIdx.z]/RHO[AC_l2__mod__cdata+i-1][vertexIdx.y][vertexIdx.z])
        }
      }
      else {
        SS[AC_l2__mod__cdata+i-1][vertexIdx.y][vertexIdx.z] = SS[AC_l2__mod__cdata-i-1][vertexIdx.y][vertexIdx.z] + (AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*(LNRHO[AC_l2__mod__cdata-i-1][vertexIdx.y][vertexIdx.z]-LNRHO[AC_l2__mod__cdata+i-1][vertexIdx.y][vertexIdx.z])
      }
    }
  }
  else {
  }
}

bc_ss_stemp_y(AcBoundary boundary,AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  if (topbot == AC_bot) {
    for i in 1:NGHOST+1 {
      dlnrho = LNRHO[vertexIdx.x][m1+i-1][vertexIdx.z]-LNRHO[vertexIdx.x][m1-i-1][vertexIdx.z]
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          dlnrho = dlnrho/(LNRHO[vertexIdx.x][m1-1][vertexIdx.z]  +reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
        }
        else {
          dlnrho = dlnrho/LNRHO[vertexIdx.x][m1-1][vertexIdx.z]
        }
      }
      SS[vertexIdx.x][m1-i-1][vertexIdx.z] = SS[vertexIdx.x][m1+i-1][vertexIdx.z] + (AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*dlnrho
    }
  }
  else if (topbot == AC_top) {
    for i in 1:NGHOST+1 {
      dlnrho = LNRHO[vertexIdx.x][AC_m2__mod__cdata+i-1][vertexIdx.z]-LNRHO[vertexIdx.x][AC_m2__mod__cdata-i-1][vertexIdx.z]
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          dlnrho = dlnrho/(LNRHO[vertexIdx.x][m2__mod__cdata-1][vertexIdx.z]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
        }
        else {
          dlnrho = dlnrho/LNRHO[vertexIdx.x][AC_m2__mod__cdata-1][vertexIdx.z]
        }
      }
      SS[vertexIdx.x][AC_m2__mod__cdata+i-1][vertexIdx.z] = SS[vertexIdx.x][AC_m2__mod__cdata-i-1][vertexIdx.z] - (AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*dlnrho
    }
  }
  else {
  }
}

bc_ss_stemp_z(AcBoundary boundary,AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  real dlnrho
  if (topbot == AC_bot) {
    for i in 1:NGHOST+1 {
      dlnrho = LNRHO[vertexIdx.x][vertexIdx.y][n1+i-1]-LNRHO[vertexIdx.x][vertexIdx.y][n1-i-1]
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          dlnrho = dlnrho/(LNRHO[vertexIdx.x][vertexIdx.y][n1-1]  +reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
        }
        else {
          dlnrho = dlnrho/LNRHO[vertexIdx.x][vertexIdx.y][n1-1]
        }
      }
      SS[vertexIdx.x][vertexIdx.y][n1-i-1] = SS[vertexIdx.x][vertexIdx.y][n1+i-1] + (AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*dlnrho
    }
  }
  else if (topbot == AC_top) {
    for i in 1:NGHOST+1 {
      dlnrho = LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+i-1]-LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-i-1]
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          dlnrho = dlnrho/(LNRHO[vertexIdx.x][vertexIdx.y][n2__mod__cdata-1]  +AC_reference_state_padded__mod__density[vertexIdx.x][iref_rho-1])
        }
        else {
          dlnrho = dlnrho/LNRHO[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-1]
        }
      }
      SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata+i-1] = SS[vertexIdx.x][vertexIdx.y][AC_n2__mod__cdata-i-1] - (AC_cp__mod__equationofstate-AC_cv__mod__equationofstate)*dlnrho
    }
  }
  else {
  }
}
#endif

bc_set_nfr_y(AcBoundary boundary, AC_TOP_BOT topbot,Field j)
{
  int k
  suppress_unused_warning(boundary)
  if (topbot == AC_bot) {
    j[vertexIdx.x][m1-1-1][vertexIdx.z]= j[vertexIdx.x][1+m1-1][vertexIdx.z]*(AC_sinth__mod__cdata[1+m1-1]/sin(AC_y__mod__cdata[1+m1-1]-AC_dy2_bound__mod__cdata[-1+NGHOST+1-1]))
    j[vertexIdx.x][m1-2-1][vertexIdx.z]= j[vertexIdx.x][2+m1-1][vertexIdx.z]*(AC_sinth__mod__cdata[2+m1-1]/sin(AC_y__mod__cdata[2+m1-1]-AC_dy2_bound__mod__cdata[-2+NGHOST+1-1]))
    j[vertexIdx.x][m1-3-1][vertexIdx.z]= j[vertexIdx.x][3+m1-1][vertexIdx.z]*(AC_sinth__mod__cdata[3+m1-1]/sin(AC_y__mod__cdata[3+m1-1]-AC_dy2_bound__mod__cdata[-3+NGHOST+1-1]))
  }
  else if (topbot == AC_top) {
    j[vertexIdx.x][1+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-1-1][vertexIdx.z]*(AC_sinth__mod__cdata[AC_m2__mod__cdata-1-1]/sin(AC_y__mod__cdata[AC_m2__mod__cdata-1-1]+AC_dy2_bound__mod__cdata[2+NGHOST-1]))
    j[vertexIdx.x][2+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-2-1][vertexIdx.z]*(AC_sinth__mod__cdata[AC_m2__mod__cdata-2-1]/sin(AC_y__mod__cdata[AC_m2__mod__cdata-2-1]+AC_dy2_bound__mod__cdata[3+NGHOST-1]))
    j[vertexIdx.x][3+AC_m2__mod__cdata-1][vertexIdx.z]= j[vertexIdx.x][AC_m2__mod__cdata-3-1][vertexIdx.z]*(AC_sinth__mod__cdata[AC_m2__mod__cdata-3-1]/sin(AC_y__mod__cdata[AC_m2__mod__cdata-3-1]+AC_dy2_bound__mod__cdata[4+NGHOST-1]))
  }
  else {
  }
}

#if Ldensity_MODULE && Lentropy_MODULE && Leos_idealgas_MODULE
bc_ss_flux_condturb_x(AcBoundary boundary,AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  int id_1
  real dsdx_yz
  real tt_yz
  real rho_yz
  real dlnrhodx_yz
  real kxbot
  if (topbot == AC_bot) {
    if (AC_pretend_lntt__mod__cdata) {
    }
    else {
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          rho_yz=LNRHO[l1-1][vertexIdx.y][vertexIdx.z]+AC_reference_state__mod__density[AC_xbot__mod__equationofstate-1][iref_rho-1]
        }
        else {
          rho_yz=LNRHO[l1-1][vertexIdx.y][vertexIdx.z]
        }
      }
      else {
        rho_yz=exp(LNRHO[l1-1][vertexIdx.y][vertexIdx.z])
      }
      if (AC_ldensity_nolog__mod__cdata) {
        tt_yz=SS[l1-1][vertexIdx.y][vertexIdx.z]
        if (AC_lreference_state__mod__cdata) {
          tt_yz = tt_yz+AC_reference_state__mod__density[AC_xbot__mod__equationofstate-1][iref_s-1]
        }
        tt_yz=AC_cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(log(rho_yz)-AC_lnrho0__mod__equationofstate)+AC_cv1__mod__equationofstate*tt_yz)
      }
      else {
        tt_yz=cs20__mod__equationofstate*exp(AC_gamma_m1__mod__equationofstate*(LNRHO[AC_l1-1][vertexIdx.y][vertexIdx.z]-lnrho0__mod__equationofstate)+AC_cv1__mod__equationofstate*SS[AC_l1-1][vertexIdx.y][vertexIdx.z])
      }
      tt_yz=tt_yz/(AC_cp__mod__equationofstate*AC_gamma_m1__mod__equationofstate)
      dlnrhodx_yz= AC_coeffs_1_x__mod__cdata[1-1][1-1]*(LNRHO[1+AC_l1-1][vertexIdx.y][vertexIdx.z]-LNRHO[l1-1-1][vertexIdx.y][vertexIdx.z])  +AC_coeffs_1_x__mod__cdata[2-1][1-1]*(LNRHO[2+AC_l1-1][vertexIdx.y][vertexIdx.z]-LNRHO[l1-2-1][vertexIdx.y][vertexIdx.z])  +AC_coeffs_1_x__mod__cdata[3-1][1-1]*(LNRHO[3+AC_l1-1][vertexIdx.y][vertexIdx.z]-LNRHO[l1-3-1][vertexIdx.y][vertexIdx.z])
      if (AC_ldensity_nolog__mod__cdata) {
        if (AC_lreference_state__mod__cdata) {
          dlnrhodx_yz=dlnrhodx_yz + AC_reference_state__mod__density[AC_xbot__mod__equationofstate-1][iref_grho-1]
          dlnrhodx_yz=dlnrhodx_yz/(rho_yz + AC_reference_state__mod__density[AC_xbot__mod__equationofstate-1][iref_rho-1])
        }
        else {
          dlnrhodx_yz=dlnrhodx_yz/rho_yz
        }
      }
      if (AC_lheatc_kramers__mod__energy) {
        kxbot=AC_hcond0_kramers__mod__energy*pow(tt_yz,(6.5*AC_nkramers__mod__energy))/pow(rho_yz,(2.*AC_nkramers__mod__energy))
      }
      else {
        kxbot=AC_hcondxbot__mod__energy
      }
      dsdx_yz=(AC_Fbot__mod__energy/tt_yz)/(AC_chit_prof1__mod__energy*AC_chi_t__mod__energy*rho_yz + kxbot*AC_cv1__mod__equationofstate)
      if (AC_lreference_state__mod__cdata) {
        dsdx_yz = dsdx_yz + AC_reference_state__mod__density[AC_xbot__mod__equationofstate-1][iref_gs-1]
      }
      for i in 1:NGHOST+1 {
        dlnrhodx_yz = LNRHO[l1+i-1][vertexIdx.y][vertexIdx.z]-LNRHO[l1-i-1][vertexIdx.y][vertexIdx.z]
        if (AC_ldensity_nolog__mod__cdata) {
          if (AC_lreference_state__mod__cdata) {
            if (l1 <= (nx+1)/2) {
              id_1 = -i
            }
            else {
              id_1 = i
            }
            dlnrhodx_yz = dlnrhodx_yz + AC_dx2_bound__mod__cdata[id_1+1+NGHOST-1]*AC_reference_state__mod__density[AC_l1-1][iref_grho-1]
          }
          dlnrhodx_yz = dlnrhodx_yz/rho_yz
        }
        SS[l1-i-1][vertexIdx.y][vertexIdx.z]=SS[l1+i-1][vertexIdx.y][vertexIdx.z] + kxbot*AC_gamma_m1__mod__equationofstate/(kxbot*AC_cv1__mod__equationofstate+AC_chit_prof1__mod__energy*AC_chi_t__mod__energy*rho_yz)*  dlnrhodx_yz + AC_dx2_bound__mod__cdata[-i+NGHOST+1-1]*dsdx_yz
      }
    }
  }
  else if (topbot == AC_top) {
  }
}
#endif

bc_set_nfr_x(AcBoundary boundary, AC_TOP_BOT topbot,Field j)
{
  suppress_unused_warning(boundary)
  int k
  if (topbot == AC_bot) {
    j[l1-1-1][vertexIdx.y][vertexIdx.z]= j[1+l1-1][vertexIdx.y][vertexIdx.z]*(AC_x__mod__cdata[1+l1-1]/(AC_x__mod__cdata[1+l1-1]-AC_dx2_bound__mod__cdata[-1+NGHOST+1-1]))
    j[l1-2-1][vertexIdx.y][vertexIdx.z]= j[2+l1-1][vertexIdx.y][vertexIdx.z]*(AC_x__mod__cdata[2+l1-1]/(AC_x__mod__cdata[2+l1-1]-AC_dx2_bound__mod__cdata[-2+NGHOST+1-1]))
    j[l1-3-1][vertexIdx.y][vertexIdx.z]= j[3+l1-1][vertexIdx.y][vertexIdx.z]*(AC_x__mod__cdata[3+l1-1]/(AC_x__mod__cdata[3+l1-1]-AC_dx2_bound__mod__cdata[-3+NGHOST+1-1]))
  }
  else if (topbot == AC_top) {
    j[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]= j[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]*(AC_x__mod__cdata[AC_l2__mod__cdata-1-1]/(AC_x__mod__cdata[AC_l2__mod__cdata-1-1]+AC_dx2_bound__mod__cdata[2+NGHOST-1]))
    j[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]= j[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]*(AC_x__mod__cdata[AC_l2__mod__cdata-2-1]/(AC_x__mod__cdata[AC_l2__mod__cdata-2-1]+AC_dx2_bound__mod__cdata[3+NGHOST-1]))
    j[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]= j[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]*(AC_x__mod__cdata[AC_l2__mod__cdata-3-1]/(AC_x__mod__cdata[AC_l2__mod__cdata-3-1]+AC_dx2_bound__mod__cdata[4+NGHOST-1]))
  }
  else {
  }
} 

set_ghosts_for_onesided_ders(AcBoundary boundary,AC_TOP_BOT topbot, Field j, int idir)
{
	suppress_unused_warning(topbot)
	suppress_unused_warning(idir)
	ac_set_ghosts_for_onesided_derivs(boundary,j)
}
