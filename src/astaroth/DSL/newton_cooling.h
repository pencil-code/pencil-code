#if LNEWTON_COOLING
Field TAU_BELOW,TAU_ABOVE
field_order(AC_itau__mod__special-1) Field F_TAU
//Reuse TAU buffer for DTAU
#define DTAU F_TAU
#define TAU F_TAU

Kernel calc_kappar_and_dtau(){
  real rho
  real tt
  real tt0
  real rho01
  real lntt0
  real kappa_cgs
  real ttdim
  real rhodim
  int i
  real k_0
  real a_0
  real b_0
  real logkk_0
  real logk_0
  const t1_0 = 132.0
  const t2_0 = 170.0
  const t3_0 = 375.0
  const t4_0 = 390.0
  const t5_0 = 580.0
  const t6_0 = 680.0
  const t7_0 = 960.0
  const t8_0 = 1570.0
  const t9_0 = 3730.0
  const t10_0 = 1e4
  const t11_0 = 1e5
  lntt0=log(AC_cs20__mod__equationofstate*AC_cp1__mod__special/AC_gamma_m1__mod__special)
  tt0=exp(lntt0)
  rho01=1./AC_rho0__mod__equationofstate
  if (AC_ldensity_nolog__mod__cdata) {
    rho=value(Field(AC_irho__mod__cdata-1))
  }
  else {
    rho=exp(value(Field(AC_ilnrho__mod__cdata-1)))
  }
  tt = tt0 *pow( (rho*rho01),AC_gamma_m1__mod__special) * exp(value(Field(AC_iss__mod__cdata-1))*AC_cv1__mod__special)
  ttdim=tt*AC_unit_temperature__mod__cdata
  rhodim=rho*AC_unit_density__mod__cdata
  if (ttdim <= t1_0) {
    k_0=2e-4
    a_0=0
    b_0= 2.1
    kappa_cgs = k_0 *pow( ttdim,b_0)
  }
  else if ((ttdim > t1_0)  &&  (ttdim <= t2_0)) {
    k_0=3.
    a_0=0
    b_0=-0.01
    kappa_cgs = k_0 *pow( ttdim,b_0)
  }
  else if ((ttdim > t2_0)  &&  (ttdim <= t3_0)) {
    k_0=0.01
    a_0=0
    b_0= 1.1
    kappa_cgs = k_0 *pow( ttdim,b_0)
  }
  else if ((ttdim > t3_0)  &&  (ttdim <= t4_0)) {
    k_0=5e4
    a_0=0
    b_0=-1.5
    kappa_cgs = k_0 *pow( ttdim,b_0)
  }
  else if ((ttdim > t4_0)  &&  (ttdim <= t5_0)) {
    k_0=0.1
    a_0=0
    b_0= 0.7
    kappa_cgs = k_0 *pow( ttdim,b_0)
  }
  else if ((ttdim > t5_0)  &&  (ttdim <= t6_0)) {
    k_0=2e15
    a_0=0
    b_0=-5.2
    kappa_cgs = k_0 *pow( ttdim,b_0)
  }
  else if ((ttdim > t6_0)  &&  (ttdim <= t7_0)) {
    k_0=0.02
    a_0=0
    b_0= 0.8
    kappa_cgs = k_0 *pow( ttdim,b_0)
  }
  else if ((ttdim > t7_0)  &&  (ttdim <= t8_0)) {
    logk_0=81.3010
    a_0=1.
    b_0=-24.
    logkk_0 = logk_0 + a_0*alog10(rhodim) + b_0*alog10(ttdim)
    kappa_cgs=pow(10,(logkk_0))
  }
  else if ((ttdim > t8_0)  &&  (ttdim <= t9_0)) {
    k_0=1e-8
    a_0=2./3
    b_0=3.
    kappa_cgs = k_0 *pow( rhodim,a_0) *pow( ttdim,b_0)
  }
  else if ((ttdim > t9_0)  &&  (ttdim <= t10_0)) {
    logk_0=-36.
    a_0=1./3
    b_0=10.
    logkk_0 = logk_0 + a_0*alog10(rhodim) + b_0*alog10(ttdim)
    kappa_cgs=pow(10,(logkk_0))
  }
  else if ((ttdim > t10_0)  &&  (ttdim <= t11_0)) {
    k_0=1.5e20
    a_0=1.
    b_0=-2.5
    kappa_cgs = k_0 *pow( rhodim,a_0) *pow( ttdim,b_0)
  }
  else if (ttdim > t11_0) {
    kappa_cgs=0.348
  }
  const real KAPPAR = kappa_cgs * (AC_unit_density__mod__cdata*AC_unit_length__mod__cdata)
  write(DTAU,KAPPAR*rho*AC_x__mod__cdata[vertexIdx.x]/AC_dy_1__mod__cdata[AC_m__mod__cdata-1])
}

Raytrace (+0,+1,+0) newton_cooling_up_ray
Raytrace (+0,-1,+0) newton_cooling_down_ray


Kernel integrate_tau_up()
{
	write(TAU_BELOW,incoming_newton_cooling_up_ray(TAU_BELOW)+TAU)
}

Kernel integrate_tau_down()
{
	write(TAU_ABOVE,incoming_newton_cooling_down_ray(TAU_ABOVE)+TAU)
}

Kernel calc_tau()
{
	write(F_TAU,min(TAU_BELOW,TAU_ABOVE))
}

#else
Kernel calc_kappar_and_dtau(){}
Kernel integrate_tau_up(){}
Kernel integrate_tau_down(){}
Kernel calc_tau(){}
#endif
