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
//TP: old
#if LENERGY
bc_ss_flux(boundary, topbot, bool lone_sided)
{
  suppress_unused_warning(boundary)
  suppress_unused_warning(lone_sided)
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

bc_copy_x(AcBoundary boundary, AC_TOP_BOT topbot,VtxBuffer j)
{
  suppress_unused_warning(boundary)
  if(topbot == AC_bot) {
    j[AC_l1-1-1][vertexIdx.y][vertexIdx.z]=j[AC_l1-1][vertexIdx.y][vertexIdx.z]
    j[AC_l1-2-1][vertexIdx.y][vertexIdx.z]=j[AC_l1-1][vertexIdx.y][vertexIdx.z]
    j[AC_l1-3-1][vertexIdx.y][vertexIdx.z]=j[AC_l1-1][vertexIdx.y][vertexIdx.z]
  }
  else if(topbot == AC_top) {
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
  if(topbot == AC_bot) {
    j[vertexIdx.x][m1-1-1][vertexIdx.z]=j[vertexIdx.x][m1-1][vertexIdx.z]
    j[vertexIdx.x][m1-2-1][vertexIdx.z]=j[vertexIdx.x][m1-1][vertexIdx.z]
    j[vertexIdx.x][m1-3-1][vertexIdx.z]=j[vertexIdx.x][m1-1][vertexIdx.z]
  }
  else if(topbot == AC_top) {
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
  if(topbot == AC_bot) {
    j[vertexIdx.x][vertexIdx.y][n1-1-1] = j[vertexIdx.x][vertexIdx.y][n1-1]
    j[vertexIdx.x][vertexIdx.y][n1-2-1] = j[vertexIdx.x][vertexIdx.y][n1-1]
    j[vertexIdx.x][vertexIdx.y][n1-3-1] = j[vertexIdx.x][vertexIdx.y][n1-1]
  }
  else if(topbot == AC_top) {
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
print("not implemented set_ghosts_for_onesided_ders");
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
print("not implemented set_ghosts_for_onesided_ders");
}
else {
for i in 1:NGHOST+1 {
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
print("not implemented set_ghosts_for_onesided_ders");
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
print("not implemented set_ghosts_for_onesided_ders");
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
  if(topbot == AC_bot) {
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
  else if(topbot == AC_top) {
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
  if(topbot == AC_bot) {
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
  else if(topbot == AC_top) {
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

/**
bc_ss_flux_x(AcBoundary boundary, AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  real work_yz
  real tmp_yz
  int i
  int stat
  real fac
  real cp_loc
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
  if(topbot == AC_bot) {
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
        if( AC_lequidist__mod__cdata.x  ||  ! heatflux_deriv_x_return_value_4_5) {
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
        if( AC_lequidist__mod__cdata.x  ||  ! heatflux_deriv_x_return_value_4_6) {
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
  else if(topbot == AC_top) {
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
        if( AC_lequidist__mod__cdata.x  ||  ! heatflux_deriv_x_return_value_4_7) {
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
        if( AC_lequidist__mod__cdata.x  ||  ! heatflux_deriv_x_return_value_4_8) {
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
bc_ss_flux_turb_x(AcBoundary boundary, AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  int i
  real rho_yz
  real cs2_yz
  real dsdx_yz
  if(topbot == AC_bot) {
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
  else if(topbot == AC_top) {
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
      dlnrhodx_yz= AC_coeffs_1_x__mod__cdata[NGHOST+1-1][2-1]*(LNRHO[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-LNRHO[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z])  +AC_coeffs_1_x__mod__cdata[NGHOST+2-1][2-1]*(LNRHO[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-LNRHO[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z])  +AC_coeffs_1_x__mod__cdata[NGHOST+3-1][2-1]*(LNRHO[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]-LNRHO[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z])
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
        SS[1+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-1-1][vertexIdx.y][vertexIdx.z]+AC_dx2_bound__mod__cdata[2+NGHOST-1]*dsdx_yz
        SS[2+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-2-1][vertexIdx.y][vertexIdx.z]+AC_dx2_bound__mod__cdata[3+NGHOST-1]*dsdx_yz
        SS[3+AC_l2__mod__cdata-1][vertexIdx.y][vertexIdx.z]=SS[AC_l2__mod__cdata-3-1][vertexIdx.y][vertexIdx.z]+AC_dx2_bound__mod__cdata[4+NGHOST-1]*dsdx_yz
      }
    }
  }
  else {
  }
}

bc_ss_temp_x(AcBoundary boundary, AC_TOP_BOT topbot)
{
  suppress_unused_warning(boundary)
  real tmp
  if(topbot == AC_bot) {
    tmp = 2/AC_gamma__mod__equationofstate*log(AC_cs2bot__mod__equationofstate/AC_cs20__mod__equationofstate)
    SS[l1-1][vertexIdx.y][vertexIdx.z] = 0.5*tmp - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[l1-1][vertexIdx.y][vertexIdx.z]-AC_lnrho0__mod__equationofstate)
    SS[l1-1-1][vertexIdx.y][vertexIdx.z] = -SS[1+l1-1][vertexIdx.y][vertexIdx.z] + tmp  - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[1+l1-1][vertexIdx.y][vertexIdx.z]+LNRHO[l1-1-1][vertexIdx.y][vertexIdx.z]-2*AC_lnrho0__mod__equationofstate)
    SS[l1-2-1][vertexIdx.y][vertexIdx.z] = -SS[2+l1-1][vertexIdx.y][vertexIdx.z] + tmp  - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[2+l1-1][vertexIdx.y][vertexIdx.z]+LNRHO[l1-2-1][vertexIdx.y][vertexIdx.z]-2*AC_lnrho0__mod__equationofstate)
    SS[l1-3-1][vertexIdx.y][vertexIdx.z] = -SS[3+l1-1][vertexIdx.y][vertexIdx.z] + tmp  - AC_gamma_m1__mod__equationofstate/AC_gamma__mod__equationofstate*(LNRHO[3+l1-1][vertexIdx.y][vertexIdx.z]+LNRHO[l1-3-1][vertexIdx.y][vertexIdx.z]-2*AC_lnrho0__mod__equationofstate)
  }
  else if(topbot == AC_top) {
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
  if(topbot == AC_bot) {
    fder=(-3*j[l1-1+1-1][vertexIdx.y][vertexIdx.z]+4*j[l1-1+2-1][vertexIdx.y][vertexIdx.z] -j[l1-1+3-1][vertexIdx.y][vertexIdx.z])/(2*AC_dx__mod__cdata)
    j[l1-1-1][vertexIdx.y][vertexIdx.z]=j[l1-1+2-1][vertexIdx.y][vertexIdx.z]-2*AC_dx__mod__cdata*fder
    fder=(-3*j[l1-2+1-1][vertexIdx.y][vertexIdx.z]+4*j[l1-2+2-1][vertexIdx.y][vertexIdx.z] -j[l1-2+3-1][vertexIdx.y][vertexIdx.z])/(2*AC_dx__mod__cdata)
    j[l1-2-1][vertexIdx.y][vertexIdx.z]=j[l1-2+2-1][vertexIdx.y][vertexIdx.z]-2*AC_dx__mod__cdata*fder
    fder=(-3*j[l1-3+1-1][vertexIdx.y][vertexIdx.z]+4*j[l1-3+2-1][vertexIdx.y][vertexIdx.z] -j[l1-3+3-1][vertexIdx.y][vertexIdx.z])/(2*AC_dx__mod__cdata)
    j[l1-3-1][vertexIdx.y][vertexIdx.z]=j[l1-3+2-1][vertexIdx.y][vertexIdx.z]-2*AC_dx__mod__cdata*fder
  }
  else if(topbot == AC_top) {
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
  if(topbot == AC_bot) {
    fder=(-3*j[vertexIdx.x][vertexIdx.y][n1-1+1-1]+4*j[vertexIdx.x][vertexIdx.y][n1-1+2-1] -j[vertexIdx.x][vertexIdx.y][n1-1+3-1])/(2*AC_dz__mod__cdata)
    j[vertexIdx.x][vertexIdx.y][n1-1-1]=j[vertexIdx.x][vertexIdx.y][n1-1+2-1]-2*AC_dz__mod__cdata*fder
    fder=(-3*j[vertexIdx.x][vertexIdx.y][n1-2+1-1]+4*j[vertexIdx.x][vertexIdx.y][n1-2+2-1] -j[vertexIdx.x][vertexIdx.y][n1-2+3-1])/(2*AC_dz__mod__cdata)
    j[vertexIdx.x][vertexIdx.y][n1-2-1]=j[vertexIdx.x][vertexIdx.y][n1-2+2-1]-2*AC_dz__mod__cdata*fder
    fder=(-3*j[vertexIdx.x][vertexIdx.y][n1-3+1-1]+4*j[vertexIdx.x][vertexIdx.y][n1-3+2-1] -j[vertexIdx.x][vertexIdx.y][n1-3+3-1])/(2*AC_dz__mod__cdata)
    j[vertexIdx.x][vertexIdx.y][n1-3-1]=j[vertexIdx.x][vertexIdx.y][n1-3+2-1]-2*AC_dz__mod__cdata*fder
  }
  else if(topbot == AC_top) {
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
  if(topbot == AC_bot) {
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
  else if(topbot == AC_top) {
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
  if(topbot == AC_bot) {
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
  else if(topbot == AC_top) {
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
