Kernel opacity(){
  real ac_transformed_pencil_acc
  real ac_transformed_pencil_ssat
  real ac_transformed_pencil_ttc
  real ac_transformed_pencil_ywater
  real ac_transformed_pencil_lambda
  real ac_transformed_pencil_chem_conc[AC_nchemspec__mod__cparam]
  real ac_transformed_pencil_nucl_rmin
  real ac_transformed_pencil_nucl_rate
  real ac_transformed_pencil_conc_satm
  real ac_transformed_pencil_ff_cond
  real ac_transformed_pencil_latent_heat
  real ac_transformed_pencil_lnrho
  real ac_transformed_pencil_rho
  real ac_transformed_pencil_rho1
  real3 ac_transformed_pencil_glnrho
  real3 ac_transformed_pencil_grho
  real ac_transformed_pencil_uglnrho
  real ac_transformed_pencil_ugrho
  real ac_transformed_pencil_glnrho2
  real ac_transformed_pencil_del2lnrho
  real ac_transformed_pencil_del2rho
  real ac_transformed_pencil_del6lnrho
  real ac_transformed_pencil_del6rho
  Matrix ac_transformed_pencil_hlnrho
  real3 ac_transformed_pencil_sglnrho
  real3 ac_transformed_pencil_uij5glnrho
  real ac_transformed_pencil_transprho
  real ac_transformed_pencil_ekin
  real ac_transformed_pencil_uuadvec_glnrho
  real ac_transformed_pencil_uuadvec_grho
  real ac_transformed_pencil_rhos1
  real3 ac_transformed_pencil_glnrhos
  real ac_transformed_pencil_totenergy_rel
  real ac_transformed_pencil_divss
  real ac_transformed_pencil_rhod[AC_ndustspec__mod__cparam]
  real3 ac_transformed_pencil_udropav
  real ac_transformed_pencil_rhodsum
  real3 ac_transformed_pencil_glnrhodsum
  real3 ac_transformed_pencil_uud[AC_ndustspec__mod__cparam]
  real ac_transformed_pencil_divud[AC_ndustspec__mod__cparam]
  Matrix ac_transformed_pencil_sdij[AC_ndustspec__mod__cparam]
  real ac_transformed_pencil_ma2
  real ac_transformed_pencil_uglntt
  real ac_transformed_pencil_ugtt
  real ac_transformed_pencil_cvspec[AC_nchemspec__mod__cparam]
  real3 ac_transformed_pencil_fpres
  real ac_transformed_pencil_tcond
  real3 ac_transformed_pencil_sglntt
  real ac_transformed_pencil_advec_cs2
  real ac_transformed_pencil_ss
  real3 ac_transformed_pencil_gss
  real ac_transformed_pencil_ee
  real ac_transformed_pencil_pp
  real ac_transformed_pencil_lntt
  real ac_transformed_pencil_cs2
  real ac_transformed_pencil_nabla_ad
  real3 ac_transformed_pencil_glntt
  real ac_transformed_pencil_tt
  real ac_transformed_pencil_tt1
  real3 ac_transformed_pencil_gtt
  real ac_transformed_pencil_yh
  real ac_transformed_pencil_del2ss
  real ac_transformed_pencil_del2lntt
  real ac_transformed_pencil_del2tt
  real ac_transformed_pencil_cv
  real ac_transformed_pencil_cv1
  real ac_transformed_pencil_cp
  real ac_transformed_pencil_cp1
  real ac_transformed_pencil_gamma
  real ac_transformed_pencil_gamma_m1
  real ac_transformed_pencil_gamma1
  real ac_transformed_pencil_mu1
  Matrix ac_transformed_pencil_hlntt
  real3 ac_transformed_pencil_rho1gpp
  real ac_transformed_pencil_delta
  real3 ac_transformed_pencil_gradcp
  real ac_transformed_pencil_del6lntt
  real3 ac_transformed_pencil_glnmumol
  real ac_transformed_pencil_ppvap
  real ac_transformed_pencil_csvap2
  real ac_transformed_pencil_rho_anel
  real3 ac_transformed_pencil_fcont[AC_n_forcing_cont_max__mod__cparam]
  real3 ac_transformed_pencil_curlfcont[AC_n_forcing_cont_max__mod__cparam]
  real3 ac_transformed_pencil_gg
  real ac_transformed_pencil_epot
  real ac_transformed_pencil_x_mn
  real ac_transformed_pencil_y_mn
  real ac_transformed_pencil_z_mn
  real ac_transformed_pencil_r_mn
  real ac_transformed_pencil_r_mn1
  real ac_transformed_pencil_phix
  real ac_transformed_pencil_phiy
  real ac_transformed_pencil_pomx
  real ac_transformed_pencil_pomy
  real ac_transformed_pencil_rcyl_mn
  real ac_transformed_pencil_rcyl_mn1
  real ac_transformed_pencil_phi_mn
  real3 ac_transformed_pencil_evr
  real3 ac_transformed_pencil_rr
  real3 ac_transformed_pencil_evth
  real ac_transformed_pencil_divu
  real3 ac_transformed_pencil_oo
  real ac_transformed_pencil_o2
  real ac_transformed_pencil_ou
  real ac_transformed_pencil_oxu2
  real3 ac_transformed_pencil_oxu
  real ac_transformed_pencil_u2
  Matrix ac_transformed_pencil_uij
  real3 ac_transformed_pencil_uu
  real3 ac_transformed_pencil_curlo
  Matrix ac_transformed_pencil_sij
  real ac_transformed_pencil_sij2
  Matrix ac_transformed_pencil_uij5
  real3 ac_transformed_pencil_ugu
  real ac_transformed_pencil_ugu2
  Matrix ac_transformed_pencil_oij
  Matrix ac_transformed_pencil_d2uidxj
  Tensor ac_transformed_pencil_uijk
  real3 ac_transformed_pencil_ogu
  real ac_transformed_pencil_u3u21
  real ac_transformed_pencil_u1u32
  real ac_transformed_pencil_u2u13
  real3 ac_transformed_pencil_del2u
  real3 ac_transformed_pencil_del4u
  real3 ac_transformed_pencil_del6u
  real ac_transformed_pencil_u2u31
  real ac_transformed_pencil_u3u12
  real ac_transformed_pencil_u1u23
  real3 ac_transformed_pencil_graddivu
  real3 ac_transformed_pencil_del6u_bulk
  real3 ac_transformed_pencil_grad5divu
  real3 ac_transformed_pencil_rhougu
  real3 ac_transformed_pencil_der6u
  real3 ac_transformed_pencil_transpurho
  real ac_transformed_pencil_divu0
  Matrix ac_transformed_pencil_u0ij
  real3 ac_transformed_pencil_uu0
  real3 ac_transformed_pencil_uu_advec
  real3 ac_transformed_pencil_uuadvec_guu
  real3 ac_transformed_pencil_del6u_strict
  real3 ac_transformed_pencil_del4graddivu
  real3 ac_transformed_pencil_uu_sph
  Matrix ac_transformed_pencil_der6u_res
  real ac_transformed_pencil_lorentz
  real ac_transformed_pencil_hless
  real ac_transformed_pencil_advec_uu
  real ac_transformed_pencil_t00
  real3 ac_transformed_pencil_t0i
  real ac_transformed_pencil_tij[6]
  real3 ac_transformed_pencil_velx
  real ac_transformed_pencil_heat
  real ac_transformed_pencil_cool
  real ac_transformed_pencil_heatcool
  real3 ac_transformed_pencil_bb
  real3 ac_transformed_pencil_bbb
  Matrix ac_transformed_pencil_bij
  real3 ac_transformed_pencil_jxbr
  real ac_transformed_pencil_ss12
  real ac_transformed_pencil_b2
  real3 ac_transformed_pencil_uxb
  real3 ac_transformed_pencil_jj
  real3 ac_transformed_pencil_aa
  real ac_transformed_pencil_diva
  real3 ac_transformed_pencil_del2a
  Matrix ac_transformed_pencil_aij
  real3 ac_transformed_pencil_bunit
  real ac_transformed_pencil_va2
  real ac_transformed_pencil_j2
  real3 ac_transformed_pencil_el
  real ac_transformed_pencil_e2
  real3 ac_transformed_pencil_uun
  real ac_transformed_pencil_divun
  Matrix ac_transformed_pencil_snij
  real ac_transformed_pencil_rhop
  real3 ac_transformed_pencil_grhop
  real ac_transformed_pencil_peh
  real ac_transformed_pencil_tauascalar
  real ac_transformed_pencil_condensationrate
  real ac_transformed_pencil_watermixingratio
  real ac_transformed_pencil_part_heatcap
  real3 ac_transformed_pencil_gcc[AC_0]
  real ac_transformed_pencil_sgs_heat
  real ac_transformed_pencil_shock
  real3 ac_transformed_pencil_gshock
  real ac_transformed_pencil_shock_perp
  real3 ac_transformed_pencil_gshock_perp
  real3 ac_transformed_pencil_fvisc
  real ac_transformed_pencil_diffus_total
  real ac_transformed_pencil_diffus_total2
  real ac_transformed_pencil_diffus_total3
  real ac_transformed_pencil_visc_heat
  real ac_transformed_pencil_nu
  real3 ac_transformed_pencil_gradnu
  real ac_transformed_pencil_nu_smag
  real3 ac_transformed_pencil_gnu_smag
  real tmp
  real lnrho
  real lntt
  real yh
  real rho
  real tt
  real profile
  real kappa1
  real kappa2
  real kappae
  real kappa_rad
  real kappa_cond
  real kappa_tot
  real kappa0
  real kappa0_cgs
  real k1
  real k2
  bool lfirst
  int i
  real lnrho__0
  real lntt__0
  real yh__0
  real mu1__0
  real tt1_0
  real tmp_0
  real tmpy_0
  real tmpy1_0
  real lnrho__1
  real lntt__1
  real yh__1
  real mu1__1
  real tt1_1
  real tmp_1
  real tmpy_1
  real tmpy1_1
  real lnrho__2
  real lntt__2
  real yh__2
  real mu1__2
  real tt1_2
  real tmp_2
  real tmpy_2
  real tmpy1_2
  real lnrho__3
  real lntt__3
  real yh__3
  real mu1__3
  real tt1_3
  real tmp_3
  real tmpy_3
  real tmpy1_3
  real lnrho__4
  real lntt__4
  real yh__4
  real mu1__4
  real tt1_4
  real tmp_4
  real tmpy_4
  real tmpy1_4
  real lnrho__5
  real lntt__5
  real yh__5
  real mu1__5
  real tt1_5
  real tmp_5
  real tmpy_5
  real tmpy1_5
  real lnrho__6
  real lntt__6
  real yh__6
  real mu1__6
  real tt1_6
  real tmp_6
  real tmpy_6
  real tmpy1_6
  real cubic_step_mn_return_value_7
  real xi_7
  real relshift_7
  real cubic_step_mn_return_value_7
  real lnrho__8
  real lntt__8
  real yh__8
  real mu1__8
  real tt1_8
  real tmp_8
  real tmpy_8
  real tmpy1_8
  real lnrho__9
  real lntt__9
  real yh__9
  real mu1__9
  real tt1_9
  real tmp_9
  real tmpy_9
  real tmpy1_9
  real lnrho__10
  real lntt__10
  real yh__10
  real mu1__10
  real tt1_10
  real tmp_10
  real tmpy_10
  real tmpy1_10
  real lnrho__11
  real lntt__11
  real yh__11
  real mu1__11
  real tt1_11
  real tmp_11
  real tmpy_11
  real tmpy1_11
  real lnrho__12
  real lntt__12
  real yh__12
  real mu1__12
  real tt1_12
  real tmp_12
  real tmpy_12
  real tmpy1_12
  if(AC_enum_opacity_type__mod__radiation == AC_enum_hminus_string__mod__cparam) {
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__0=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__0=value(Field(AC_ilntt__mod__cdata-1))
      yh__0=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__0=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__0=value(Field(AC_ilntt__mod__cdata-1))
      yh__0=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (false) {
      tmp_0 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__0) - lnrho__0
    }
    if (present(tmp)) {
      tt1_0 = exp(-lntt__0)
      tmp_0 = 2*lnrho__0-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__0)+AC_tt_ion___mod__equationofstate*tt1_0
      tmpy_0 = yh__0+AC_ymetals__mod__equationofstate
      if ( tmpy_0==0. ) {
        tmp=0.
      }
      else {
        tmp = (1-yh__0)*AC_kappa0__mod__equationofstate*exp(min(tmp_0,log(AC_huge1__mod__cparam))+log(tmpy_0))
      }
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__0 = AC_mu1_0__mod__equationofstate*(1+yh__0+AC_xhe__mod__equationofstate)
        tmpy1_0 = 1./(1+yh__0)
        tmp = tmp*(1+4*AC_xhe__mod__equationofstate)*mu1__0*tmpy1_0
      }
    }
    DF_KAPPARHO=AC_kapparho_floor__mod__radiation+tmp*AC_scalefactor_kappa__mod__radiation[inu-1]
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_total_rosseland_mean_string__mod__cparam) {
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__1=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__1=value(Field(AC_ilntt__mod__cdata-1))
      yh__1=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__1=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__1=value(Field(AC_ilntt__mod__cdata-1))
      yh__1=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lnrho)) {
      lnrho=lnrho__1
    }
    if (false) {
      tmp_1 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__1) - lnrho__1
    }
    if (false) {
      tt1_1 = exp(-lntt__1)
      tmp_1 = 2*lnrho__1-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__1)+AC_tt_ion___mod__equationofstate*tt1_1
      tmpy_1 = yh__1+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__1 = AC_mu1_0__mod__equationofstate*(1+yh__1+AC_xhe__mod__equationofstate)
        tmpy1_1 = 1./(1+yh__1)
      }
    }
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__2=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__2=value(Field(AC_ilntt__mod__cdata-1))
      yh__2=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__2=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__2=value(Field(AC_ilntt__mod__cdata-1))
      yh__2=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lntt)) {
      lntt=lntt__2
    }
    if (false) {
      tmp_2 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__2) - lnrho__2
    }
    if (false) {
      tt1_2 = exp(-lntt__2)
      tmp_2 = 2*lnrho__2-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__2)+AC_tt_ion___mod__equationofstate*tt1_2
      tmpy_2 = yh__2+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__2 = AC_mu1_0__mod__equationofstate*(1+yh__2+AC_xhe__mod__equationofstate)
        tmpy1_2 = 1./(1+yh__2)
      }
    }
    kappa1=4.0e25*1.7381*0.0135*(AC_unit_density__mod__cdata*AC_unit_density__mod__cdata)*AC_unit_length__mod__cdata*  exp(lnrho)*pow((exp(lntt)*AC_unit_temperature__mod__cdata),(-3.5))
    kappa2=1.25d-29*0.0134*pow(AC_unit_density__mod__cdata,1.5)*AC_unit_length__mod__cdata*( AC_unit_temperature__mod__cdata* AC_unit_temperature__mod__cdata* AC_unit_temperature__mod__cdata* AC_unit_temperature__mod__cdata* AC_unit_temperature__mod__cdata* AC_unit_temperature__mod__cdata* AC_unit_temperature__mod__cdata* AC_unit_temperature__mod__cdata* AC_unit_temperature__mod__cdata)*exp(0.5*lnrho)*exp(9.0*lntt)
    kappae=0.2*1.7381*pow((1.+2.7e11*exp(lnrho-2*lntt)*AC_unit_density__mod__cdata/(AC_unit_temperature__mod__cdata*AC_unit_temperature__mod__cdata)),(-1.))
    kappa_cond=2.6e-7*AC_unit_length__mod__cdata*(AC_unit_temperature__mod__cdata*AC_unit_temperature__mod__cdata)*exp(2*lntt)*exp(-lnrho)
    kappa_rad=AC_kapparho_floor__mod__radiation+1./(1./(kappa1+kappae)+1./kappa2)
    kappa_tot=1./(1./kappa_rad+1./kappa_cond)
    if (AC_lcutoff_opticallythin__mod__radiation  &&  AC_z__mod__cdata[AC_n__mod__cdata-1] > AC_z_cutoff__mod__radiation) {
      kappa_tot=0.5*(1.-tanh((lntt-log(1.0e4))/log(2.0)))/(1./kappa_rad+1./kappa_cond)
    }
    kappa_tot = min(kappa_tot,AC_kappa_ceiling__mod__radiation)
    DF_KAPPARHO=exp(lnrho)*kappa_tot*AC_scalefactor_kappa__mod__radiation[inu-1]
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_kappa_es_string__mod__cparam) {
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__3=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__3=value(Field(AC_ilntt__mod__cdata-1))
      yh__3=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__3=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__3=value(Field(AC_ilntt__mod__cdata-1))
      yh__3=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lnrho)) {
      lnrho=lnrho__3
    }
    if (false) {
      tmp_3 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__3) - lnrho__3
    }
    if (false) {
      tt1_3 = exp(-lntt__3)
      tmp_3 = 2*lnrho__3-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__3)+AC_tt_ion___mod__equationofstate*tt1_3
      tmpy_3 = yh__3+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__3 = AC_mu1_0__mod__equationofstate*(1+yh__3+AC_xhe__mod__equationofstate)
        tmpy1_3 = 1./(1+yh__3)
      }
    }
    DF_KAPPARHO=AC_kapparho_floor__mod__radiation+AC_kappa_es__mod__cdata*exp(lnrho)
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_kappa_cst_string__mod__cparam) {
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__4=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__4=value(Field(AC_ilntt__mod__cdata-1))
      yh__4=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__4=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__4=value(Field(AC_ilntt__mod__cdata-1))
      yh__4=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lnrho)) {
      lnrho=lnrho__4
    }
    if (false) {
      tmp_4 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__4) - lnrho__4
    }
    if (false) {
      tt1_4 = exp(-lntt__4)
      tmp_4 = 2*lnrho__4-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__4)+AC_tt_ion___mod__equationofstate*tt1_4
      tmpy_4 = yh__4+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__4 = AC_mu1_0__mod__equationofstate*(1+yh__4+AC_xhe__mod__equationofstate)
        tmpy1_4 = 1./(1+yh__4)
      }
    }
    DF_KAPPARHO=AC_kapparho_floor__mod__radiation+AC_kappa_cst__mod__radiation[inu-1]*exp(lnrho)
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_kapparho_cst_string__mod__cparam) {
    DF_KAPPARHO=AC_kapparho_floor__mod__radiation+AC_kapparho_cst__mod__radiation
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_kappa_kconst_string__mod__cparam) {
    kappa0=16./3.*AC_sigmasb__mod__cdata/AC_kappa_kconst__mod__radiation
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__5=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__5=value(Field(AC_ilntt__mod__cdata-1))
      yh__5=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__5=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__5=value(Field(AC_ilntt__mod__cdata-1))
      yh__5=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lntt)) {
      lntt=lntt__5
    }
    if (false) {
      tmp_5 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__5) - lnrho__5
    }
    if (false) {
      tt1_5 = exp(-lntt__5)
      tmp_5 = 2*lnrho__5-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__5)+AC_tt_ion___mod__equationofstate*tt1_5
      tmpy_5 = yh__5+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__5 = AC_mu1_0__mod__equationofstate*(1+yh__5+AC_xhe__mod__equationofstate)
        tmpy1_5 = 1./(1+yh__5)
      }
    }
    tt=exp(lntt)
    DF_KAPPARHO=AC_kapparho_floor__mod__radiation+kappa0*(tt*tt*tt)
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_kappa_power_law_string__mod__cparam) {
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__6=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__6=value(Field(AC_ilntt__mod__cdata-1))
      yh__6=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__6=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__6=value(Field(AC_ilntt__mod__cdata-1))
      yh__6=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lnrho)) {
      lnrho=lnrho__6
    }
    if (present(lntt)) {
      lntt=lntt__6
    }
    if (false) {
      tmp_6 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__6) - lnrho__6
    }
    if (false) {
      tt1_6 = exp(-lntt__6)
      tmp_6 = 2*lnrho__6-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__6)+AC_tt_ion___mod__equationofstate*tt1_6
      tmpy_6 = yh__6+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__6 = AC_mu1_0__mod__equationofstate*(1+yh__6+AC_xhe__mod__equationofstate)
        tmpy1_6 = 1./(1+yh__6)
      }
    }
    rho=exp(lnrho)
    tt=exp(lntt)
    if (AC_knee_temp_opa__mod__radiation==0.0) {
      DF_KAPPARHO=AC_kapparho_floor__mod__radiation+rho*AC_kappa_cst__mod__radiation[inu-1]*pow(  (rho/AC_ref_rho_opa__mod__radiation),AC_expo_rho_opa__mod__radiation)*pow(  (tt/AC_ref_temp_opa__mod__radiation),AC_expo_temp_opa__mod__radiation)
    }
    else {
      if (false) {
      }
      else {
        relshift_7=0.0
      }
      xi_7 = (tt-AC_knee_temp_opa__mod__radiation)/(AC_width_temp_opa__mod__radiation+AC_tini__mod__cparam) - relshift_7
      xi_7 = max(xi_7,-1.0)
      xi_7 = min(xi_7,1.0)
      cubic_step_mn_return_value_7 = 0.5 + xi_7*(0.75-(xi_7*xi_7)*0.25)
      profile=1.0-cubic_step_mn_return_value_7
      DF_KAPPARHO=AC_kapparho_floor__mod__radiation+profile*rho*AC_kappa_cst__mod__radiation[inu-1]*pow(  (rho/AC_ref_rho_opa__mod__radiation),AC_expo_rho_opa__mod__radiation)*pow(  (tt/AC_ref_temp_opa__mod__radiation),AC_expo_temp_opa__mod__radiation) +  (1.0-profile)*rho*AC_kappa_cst__mod__radiation[inu-1]*pow(  (AC_knee_temp_opa__mod__radiation/AC_ref_temp_opa__mod__radiation),AC_expo_temp_opa__mod__radiation)*pow(  (tt/AC_knee_temp_opa__mod__radiation),AC_expo_temp_opa_buff__mod__radiation)
    }
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_kappa_double_power_law_string__mod__cparam) {
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__8=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__8=value(Field(AC_ilntt__mod__cdata-1))
      yh__8=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__8=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__8=value(Field(AC_ilntt__mod__cdata-1))
      yh__8=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lnrho)) {
      lnrho=lnrho__8
    }
    if (present(lntt)) {
      lntt=lntt__8
    }
    if (false) {
      tmp_8 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__8) - lnrho__8
    }
    if (false) {
      tt1_8 = exp(-lntt__8)
      tmp_8 = 2*lnrho__8-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__8)+AC_tt_ion___mod__equationofstate*tt1_8
      tmpy_8 = yh__8+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__8 = AC_mu1_0__mod__equationofstate*(1+yh__8+AC_xhe__mod__equationofstate)
        tmpy1_8 = 1./(1+yh__8)
      }
    }
    rho=exp(lnrho)
    tt=exp(lntt)
    kappa1=AC_kappa_cst__mod__radiation[inu-1]*pow((rho/AC_ref_rho_opa__mod__radiation),AC_expo_rho_opa__mod__radiation)*pow(  (tt/AC_ref_temp_opa__mod__radiation),AC_expo_temp_opa__mod__radiation)
    kappa2=AC_kappa20_cst__mod__radiation[inu-1]+AC_kappa_cst__mod__radiation[inu-1]*pow((rho/AC_ref_rho_opa__mod__radiation),AC_expo2_rho_opa__mod__radiation)*pow(  (tt/AC_ref_temp_opa__mod__radiation),AC_expo2_temp_opa__mod__radiation)
    DF_KAPPARHO= AC_kapparho_floor__mod__radiation+rho/(1./kappa1+1./kappa2)  *(1.+AC_ampl_bump__mod__radiation*exp(-0.5*(((tt-AC_tt_bump__mod__radiation)/AC_sigma_bump__mod__radiation)*((tt-AC_tt_bump__mod__radiation)/AC_sigma_bump__mod__radiation))))
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_tsquare_string__mod__cparam) {
    kappa0_cgs=2e-4
    kappa0=kappa0_cgs
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__9=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__9=value(Field(AC_ilntt__mod__cdata-1))
      yh__9=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__9=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__9=value(Field(AC_ilntt__mod__cdata-1))
      yh__9=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lnrho)) {
      lnrho=lnrho__9
    }
    if (present(lntt)) {
      lntt=lntt__9
    }
    if (false) {
      tmp_9 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__9) - lnrho__9
    }
    if (false) {
      tt1_9 = exp(-lntt__9)
      tmp_9 = 2*lnrho__9-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__9)+AC_tt_ion___mod__equationofstate*tt1_9
      tmpy_9 = yh__9+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__9 = AC_mu1_0__mod__equationofstate*(1+yh__9+AC_xhe__mod__equationofstate)
        tmpy1_9 = 1./(1+yh__9)
      }
    }
    DF_KAPPARHO=AC_kapparho_floor__mod__radiation+exp(lnrho)*kappa0*(((exp(lntt))*(exp(lntt))))
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_kramers_string__mod__cparam) {
    kappa0_cgs=6.6e22
    kappa0=kappa0_cgs
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__10=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__10=value(Field(AC_ilntt__mod__cdata-1))
      yh__10=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__10=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__10=value(Field(AC_ilntt__mod__cdata-1))
      yh__10=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lnrho)) {
      lnrho=lnrho__10
    }
    if (present(lntt)) {
      lntt=lntt__10
    }
    if (false) {
      tmp_10 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__10) - lnrho__10
    }
    if (false) {
      tt1_10 = exp(-lntt__10)
      tmp_10 = 2*lnrho__10-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__10)+AC_tt_ion___mod__equationofstate*tt1_10
      tmpy_10 = yh__10+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__10 = AC_mu1_0__mod__equationofstate*(1+yh__10+AC_xhe__mod__equationofstate)
        tmpy1_10 = 1./(1+yh__10)
      }
    }
    DF_KAPPARHO=AC_kapparho_floor__mod__radiation+kappa0*((exp(lnrho)*exp(lnrho)))*pow((exp(lntt)),(-3.5))
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_dustzinfrared_string__mod__cparam) {
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__11=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__11=value(Field(AC_ilntt__mod__cdata-1))
      yh__11=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__11=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__11=value(Field(AC_ilntt__mod__cdata-1))
      yh__11=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lnrho)) {
      lnrho=lnrho__11
    }
    if (present(lntt)) {
      lntt=lntt__11
    }
    if (false) {
      tmp_11 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__11) - lnrho__11
    }
    if (false) {
      tt1_11 = exp(-lntt__11)
      tmp_11 = 2*lnrho__11-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__11)+AC_tt_ion___mod__equationofstate*tt1_11
      tmpy_11 = yh__11+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__11 = AC_mu1_0__mod__equationofstate*(1+yh__11+AC_xhe__mod__equationofstate)
        tmpy1_11 = 1./(1+yh__11)
      }
    }
    if (exp(lntt) <= 150) {
      tmp=2e-4*(exp(lntt)*exp(lntt))
    }
    else if(exp(lntt)>=200) {
      tmp=exp(0.861353*lntt-4.56372)
    }
    else {
      tmp=exp(-5.22826*lntt+27.7010)
    }
    DF_KAPPARHO=AC_kapparho_floor__mod__radiation+exp(lnrho)*tmp
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_blob_string__mod__cparam) {
    if (lfirst) {
      DF_KAPPARHO=AC_kapparho_floor__mod__radiation+AC_kapparho_const__mod__radiation + AC_amplkapparho__mod__radiation  *exp(-((AC_x__mod__cdata[vertexIdx.x]/AC_radius_kapparho__mod__radiation)*(AC_x__mod__cdata[vertexIdx.x]/AC_radius_kapparho__mod__radiation)))  *exp(-((AC_y__mod__cdata[vertexIdx.y]/AC_radius_kapparho__mod__radiation)*(AC_y__mod__cdata[vertexIdx.y]/AC_radius_kapparho__mod__radiation)))  *exp(-((AC_z__mod__cdata[vertexIdx.z]/AC_radius_kapparho__mod__radiation)*(AC_z__mod__cdata[vertexIdx.z]/AC_radius_kapparho__mod__radiation)))
      lfirst=false
    }
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_cos_string__mod__cparam) {
    if (lfirst) {
      DF_KAPPARHO=AC_kapparho_floor__mod__radiation+AC_kapparho_const__mod__radiation + AC_amplkapparho__mod__radiation  *cos(AC_kx_kapparho__mod__radiation*AC_x__mod__cdata[vertexIdx.x])  *cos(AC_ky_kapparho__mod__radiation*AC_y__mod__cdata[vertexIdx.y])  *cos(AC_kz_kapparho__mod__radiation*AC_z__mod__cdata[vertexIdx.z])
      lfirst=false
    }
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_rad_ionization_string__mod__cparam) {
    if(AC_mx__mod__cparam == AC_nx__mod__cparam) {
      lnrho__12=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__12=value(Field(AC_ilntt__mod__cdata-1))
      yh__12=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(AC_mx__mod__cparam == AC_mx__mod__cparam) {
      lnrho__12=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__12=value(Field(AC_ilntt__mod__cdata-1))
      yh__12=value(Field(AC_iyh__mod__cdata-1))
    }
    else {
    }
    if (present(lnrho)) {
      lnrho=lnrho__12
    }
    if (present(yh)) {
      yh = yh__12
    }
    if (false) {
      tmp_12 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__12) - lnrho__12
    }
    if (false) {
      tt1_12 = exp(-lntt__12)
      tmp_12 = 2*lnrho__12-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__12)+AC_tt_ion___mod__equationofstate*tt1_12
      tmpy_12 = yh__12+AC_ymetals__mod__equationofstate
      if (AC_lhminus_opacity_correction__mod__equationofstate) {
        mu1__12 = AC_mu1_0__mod__equationofstate*(1+yh__12+AC_xhe__mod__equationofstate)
        tmpy1_12 = 1./(1+yh__12)
      }
    }
    DF_KAPPARHO=AC_kapparho_floor__mod__radiation+ AC_sigmah___mod__cdata*(1 - yh)*exp(2*lnrho+log(AC_m_h__mod__cdata))
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_b2_string__mod__cparam) {
    fatal_error_message(true,"calc_kapparho_b2")
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_b2zw2_string__mod__cparam) {
    fatal_error_message(true,"calc_kapparho_b2_w2")
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_read_file_string__mod__cparam) {
  }
  else if(AC_enum_opacity_type__mod__radiation == AC_enum_nothing_string__mod__cparam) {
    DF_KAPPARHO=0.0
  }
  else {
  }
  write(F_KAPPARHO,DF_KAPPARHO)
}
