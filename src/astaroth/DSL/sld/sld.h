#if Lentropy_MODULE
energy_sld_sound_speed(){
  real cs2
  real prof_cs
  real fact_rho
  real fact_wsld
  real rhotop
  real lnrhotop
  real w_sldrat2
  real gamma
  real gamma_m1
  real cv
  real cv1
  real get_cp_return_value_0_2
  real get_cv_return_value_1_2
  real step_vector_return_value_3
  real step_vector_return_value_3
  real step_vector_return_value_3
  real DF_SLD_CHAR_ENTROPY
  if (AC_lslope_limit_diff__mod__cdata  &&  AC_llast__mod__cdata) {
    if (present(gamma)) {
      gamma=AC_gamma__mod__equationofstate
    }
    if (false) {
      get_cp_return_value_0_2=AC_cp__mod__equationofstate
    }
    if (present(cv)) {
      get_cv_return_value_1_2=AC_cv__mod__equationofstate
      cv=get_cv_return_value_1_2
    }
    gamma_m1=gamma-1.
    cv1=1./cv
    if (AC_ldensity_nolog__mod__cdata) {
      rhotop=exp(AC_lnrho0__mod__equationofstate)*(pow((AC_cs2top__mod__equationofstate/AC_cs20__mod__equationofstate),(1./gamma_m1)))
    }
    else {
      lnrhotop=AC_lnrho0__mod__equationofstate+(1./gamma_m1)*log(AC_cs2top__mod__equationofstate/AC_cs20__mod__equationofstate)
    }
    cs2=0.
    prof_cs=1.
    fact_rho=1.
    fact_wsld=1.
    if (AC_lsld_char_cslimit__mod__energy) {
      w_sldrat2=pow(AC_w_sldchar_ene2__mod__energy,2.)/(pow(AC_w_sldchar_ene__mod__energy,2.)+AC_tini__mod__cparam)
    }
    else {
      w_sldrat2=1.0
    }
    if (AC_lsld_char_wprofr__mod__energy) {
      fact_wsld=1 + (AC_w_sldchar_ene2__mod__energy/AC_w_sldchar_ene__mod__energy -1.)  *pow((AC_x__mod__cdata[vertexIdx.x]/AC_w_sldchar_ene_r0__mod__energy),AC_w_sldchar_ene_p__mod__energy)
    }
    if (AC_ldensity_nolog__mod__cdata) {
      cs2 = AC_cs20__mod__equationofstate*exp(gamma_m1*(log(value(Field(AC_irho__mod__cdata-1)))-AC_lnrho0__mod__equationofstate)+cv1*value(Field(AC_iss__mod__cdata-1)))
      if (AC_lsld_char_rholimit__mod__energy) {
        fact_rho=1.+pow((rhotop/value(Field(AC_irho__mod__cdata-1))),AC_cs2top__mod__equationofstate)/cs2
      }
    }
    else {
      cs2 = AC_cs20__mod__equationofstate*exp(gamma_m1*(value(Field(AC_ilnrho__mod__cdata-1))-AC_lnrho0__mod__equationofstate)+cv1*value(Field(AC_iss__mod__cdata-1)))
      if (AC_lsld_char_rholimit__mod__energy) {
        fact_rho=1.+ exp((lnrhotop-value(Field(AC_ilnrho__mod__cdata-1))))*AC_cs2top__mod__equationofstate/cs2
      }
    }
    if (AC_lsld_char_cslimit__mod__energy) {
      if (AC_w_sldchar_ene__mod__energy==0.0) {
        DF_SLD_CHAR_ENTROPY=AC_w_sldchar_ene2__mod__energy*sqrt(AC_cs2top__mod__equationofstate)
      }
      else {
        step_vector_return_value_3 = 0.5*(1+tanh((cs2-w_sldrat2*AC_cs2top__mod__equationofstate+AC_cs2top__mod__equationofstate/200.)/(AC_cs2top__mod__equationofstate/200.+AC_tini__mod__cparam)))
        prof_cs=step_vector_return_value_3
        DF_SLD_CHAR_ENTROPY =AC_w_sldchar_ene__mod__energy*prof_cs*sqrt(cs2)  + (1.-prof_cs)*AC_w_sldchar_ene2__mod__energy*sqrt(AC_cs2top__mod__equationofstate)
      }
    }
    else {
      DF_SLD_CHAR_ENTROPY =AC_w_sldchar_ene__mod__energy*fact_wsld*sqrt(cs2*fact_rho)
    }
  }
  return DF_SLD_CHAR_ENTROPY
}
#else
energy_sld_sound_speed(){
	return AC_w_sldchar_ene__mod__energy*AC_cs0__mod__equationofstate
}
#endif

Kernel sld_calc_char_speed(PC_SUB_STEP_NUMBER step_num)
{
#if Lmagnetic_MODULE && Leos_idealgas_MODULE && Lhydro_MODULE
	if(AC_lslope_limit_diff__mod__cdata && step_num == AC_num_substeps__mod__cdata-1)
	{
		res = calculate_characteristic_speed(AC_w_sldchar_hyd__mod__hydro, UU, 1.0, energy_sld_sound_speed(), AC_w_sldchar_mag__mod__magnetic, curl(AA), LNRHO, AC_mu0__mod__cdata)
		write(SLD_CHAR_SPEED,res)
	}
#endif
}
