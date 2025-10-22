#if Leos_ionization_MODULE
Kernel ioncalc(){
  real lnrho
  real ss
  real yh
  real lntt
  real dyhold_0
  real dyh_0
  real yhlow_0
  real yhhigh_0
  real ff_0
  real dff_0
  real lntt__0
  real dlntt__0
  real tt1__0
  real fractions1_0
  bool found_0
  int i_0
  const maxit_0 = 1000
  if (AC_ldensity_nolog__mod__cdata) {
    lnrho=log(value(Field(AC_ilnrho__mod__cdata-1)))
  }
  else {
    lnrho=value(Field(AC_ilnrho__mod__cdata-1))
  }
  ss=value(Field(AC_iss__mod__cdata-1))
  yh=value(Field(AC_iyh__mod__cdata-1))
  yhlow_0=AC_yhmin__mod__equationofstate
  yhhigh_0=AC_yhmax__mod__equationofstate
  dyh_0=yhhigh_0-yhlow_0
  dyhold_0=dyh_0
  found_0=false
  fractions1_0=1/(1+yh+AC_xhe__mod__equationofstate)
  lntt__0=(2.0/3.0)*((ss/AC_ss_ion__mod__equationofstate+(1-yh)*(log(1-yh+epsi)-AC_lnrho_h__mod__equationofstate)  +yh*(2*log(yh)-AC_lnrho_e__mod__equationofstate-AC_lnrho_h__mod__equationofstate)  +AC_xhe_term__mod__equationofstate)*fractions1_0+lnrho-2.5)
  tt1__0=exp(-lntt__0)
  ff_0=AC_lnrho_e__mod__equationofstate-lnrho+1.5*lntt__0-tt1__0+log(1-yh+epsi)-2*log(yh)
  dlntt__0=((2.0/3.0)*(-ff_0-tt1__0)-1)*fractions1_0
  dff_0=dlntt__0*(1.5+tt1__0)-1/(1-yh+epsi)-2/(yh+epsi)
  i_0 = 1
  while(i_0 <= maxit_0  &&  !found_0){
    i_0 = i_0 + 1
    if (!found_0) {
      if (      sign(1.,((yh-yhlow_0)*dff_0-ff_0))  == sign(1.,((yh-yhhigh_0)*dff_0-ff_0))   ||  abs(2*ff_0) > abs(dyhold_0*dff_0) ) {
        dyhold_0=dyh_0
        dyh_0=0.5*(yhhigh_0-yhlow_0)
        yh=yhhigh_0-dyh_0
      }
      else {
        dyhold_0=dyh_0
        dyh_0=ff_0/dff_0
        dyh_0=min(dyh_0,yh-AC_yhmin__mod__equationofstate)
        dyh_0=max(dyh_0,yh-AC_yhmax__mod__equationofstate)
        yh=yh-dyh_0
      }
    }
    if (abs(dyh_0)>max(AC_yhacc__mod__equationofstate,1e-31)*max(yh,1e-31)) {
      fractions1_0=1/(1+yh+AC_xhe__mod__equationofstate)
      lntt__0=(2.0/3.0)*((ss/AC_ss_ion__mod__equationofstate+(1-yh)*(log(1-yh+epsi)-AC_lnrho_h__mod__equationofstate)  +yh*(2*log(yh)-AC_lnrho_e__mod__equationofstate-AC_lnrho_h__mod__equationofstate)  +AC_xhe_term__mod__equationofstate)*fractions1_0+lnrho-2.5)
      tt1__0=exp(-lntt__0)
      ff_0=AC_lnrho_e__mod__equationofstate-lnrho+1.5*lntt__0-tt1__0+log(1-yh+epsi)-2*log(yh)
      dlntt__0=((2.0/3.0)*(-ff_0-tt1__0)-1)*fractions1_0
      dff_0=dlntt__0*(1.5+tt1__0)-1/(1-yh+epsi)-2/yh
      if (ff_0<0) {
        yhhigh_0=yh
      }
      else {
        yhlow_0=yh
      }
    }
    else {
      found_0=true
    }
  }
  lntt=(ss/AC_ss_ion__mod__equationofstate+(1-yh)*(log(1-yh+epsi)-AC_lnrho_h__mod__equationofstate)  +yh*(2*log(yh)-AC_lnrho_e__mod__equationofstate-AC_lnrho_h__mod__equationofstate)+AC_xhe_term__mod__equationofstate)/(1+yh+AC_xhe__mod__equationofstate)
  lntt=(2.0/3.0)*(lntt+lnrho-2.5)+AC_lntt_ion__mod__equationofstate
  write(F_YH,yh)
  write(F_LNTT,lntt)
}

#else

#if Leos_temperature_ionization_MODULE
Kernel ioncalc(){
  real yh
  real rho1
  real tt1
  real rhs
  real sqrtrhs
  real mu1
  real yh_term_cp
  real tt_term_cp

  real DF_YH
  if (AC_lconst_yh__mod__equationofstate) {
    DF_YH = AC_yh_const__mod__equationofstate
  }
  else {
    rho1 = exp(-value(Field(AC_ilnrho__mod__cdata-1)))
    tt1 = exp(-value(Field(AC_ilntt__mod__cdata-1)))
    rhs = AC_rho_e__mod__equationofstate*rho1*pow((tt1*AC_tt_ion__mod__equationofstate),(-1.5))*exp(-AC_tt_ion__mod__equationofstate*tt1)
    sqrtrhs = sqrt(rhs)
    yh = 2*sqrtrhs/(sqrtrhs+sqrt(4+rhs))
    DF_YH = yh
    if (AC_lcalc_cp_full__mod__equationofstate) {
      mu1 = AC_mu1_0__mod__equationofstate*(1 + yh + AC_xhe__mod__equationofstate)
      yh_term_cp = yh*(1-yh)/(2+AC_xhe__mod__equationofstate*(2-yh))
      tt_term_cp = 2.5 + AC_tt_ion__mod__equationofstate*tt1
      cp_full__mod__equationofstate = AC_rgas__mod__equationofstate*mu1*(2.5 + yh_term_cp*(tt_term_cp*tt_term_cp))
      write(AC_cp_full__mod__equationofstate,cp_full__mod__equationofstate)
    }
  }
  write(F_YH,DF_YH)
}
#else

Kernel ioncalc(){
}

#endif
#endif
