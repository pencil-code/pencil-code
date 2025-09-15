#if LALPHADISK
field_order(AC_isigma__mod__alphadisk-1) Field F_SIGMA
field_order(AC_imdot__mod__alphadisk-1) Field F_MDOT
field_order(AC_itmid__mod__alphadisk-1) Field F_TMID
Kernel after_timestep_alphadisk(){
  real DF_MDOT = 0.0
  real DF_TMID = 0.0

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
  real ac_transformed_pencil_rho
  real ac_transformed_pencil_lnrho
  real ac_transformed_pencil_rho1
  real3 ac_transformed_pencil_glnrho
  real ac_transformed_pencil_del2rho
  real ac_transformed_pencil_del2lnrho
  Matrix ac_transformed_pencil_hlnrho
  real3 ac_transformed_pencil_grho
  real ac_transformed_pencil_glnrho2
  real ac_transformed_pencil_del6lnrho
  real3 ac_transformed_pencil_uij5glnrho
  real ac_transformed_pencil_uglnrho
  real ac_transformed_pencil_ugrho
  real3 ac_transformed_pencil_sglnrho
  real ac_transformed_pencil_ekin
  real ac_transformed_pencil_transprho
  real3 ac_transformed_pencil_glnrhos
  real ac_transformed_pencil_totenergy_rel
  real ac_transformed_pencil_rhod[AC_ndustspec__mod__cparam]
  real3 ac_transformed_pencil_udropav
  real ac_transformed_pencil_rhodsum
  real3 ac_transformed_pencil_glnrhodsum
  real3 ac_transformed_pencil_uud[AC_ndustspec__mod__cparam]
  real ac_transformed_pencil_divud[AC_ndustspec__mod__cparam]
  Matrix ac_transformed_pencil_sdij[AC_ndustspec__mod__cparam]
  real ac_transformed_pencil_ma2
  real3 ac_transformed_pencil_fpres
  real ac_transformed_pencil_tcond
  real3 ac_transformed_pencil_sglntt
  real ac_transformed_pencil_uglntt
  real ac_transformed_pencil_advec_cs2
  real ac_transformed_pencil_ss
  real3 ac_transformed_pencil_gss
  real ac_transformed_pencil_ee
  real ac_transformed_pencil_pp
  real ac_transformed_pencil_lntt
  real ac_transformed_pencil_cs2
  real ac_transformed_pencil_cv1
  real ac_transformed_pencil_cp1
  real ac_transformed_pencil_cp1tilde
  real3 ac_transformed_pencil_glntt
  real ac_transformed_pencil_tt
  real ac_transformed_pencil_tt1
  real ac_transformed_pencil_cp
  real ac_transformed_pencil_cv
  real3 ac_transformed_pencil_gtt
  real ac_transformed_pencil_mu1
  real3 ac_transformed_pencil_glnmu
  real3 ac_transformed_pencil_gmu1
  real ac_transformed_pencil_yh
  Matrix ac_transformed_pencil_hss
  Matrix ac_transformed_pencil_hlntt
  real ac_transformed_pencil_del2ss
  real ac_transformed_pencil_del6ss
  real ac_transformed_pencil_del2tt
  real ac_transformed_pencil_del2lntt
  real ac_transformed_pencil_del6tt
  real ac_transformed_pencil_del6lntt
  real3 ac_transformed_pencil_glnmumol
  real ac_transformed_pencil_csvap2
  real ac_transformed_pencil_rho_anel
  real3 ac_transformed_pencil_rho1gpp
  real3 ac_transformed_pencil_fcont[AC_n_forcing_cont_max__mod__cparam]
  real3 ac_transformed_pencil_gg
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
  real ac_transformed_pencil_lorentz_gamma2
  real ac_transformed_pencil_lorentz_gamma
  real ac_transformed_pencil_ss_rel2
  real3 ac_transformed_pencil_ss_rel
  Matrix ac_transformed_pencil_ss_rel_ij
  real ac_transformed_pencil_ss_rel_factor
  real ac_transformed_pencil_divss_rel
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
  real3 ac_transformed_pencil_mf_emf
  real ac_transformed_pencil_mf_emfdotb
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
  real ac_transformed_pencil_visc_heat
  real ac_transformed_pencil_nu
  real3 ac_transformed_pencil_gradnu
  real ac_transformed_pencil_nu_smag
  real lgsigma_0
  real lgsigma1_0
  real lgsigma2_0
  real lgmdot_0
  real sig_1
  real sdo_1
  real sup_1
  real lnsig_1
  real lnsdo_1
  real lnsup_1
  int isig_do_1
  int isig_up_1
  real lgsigma_2
  real lgsigma1_2
  real lgsigma2_2
  real lgmdot_2
  int i_2
  if(enum_temperature_model__mod__alphadisk == enum_hayashi_string) {
    if(enum_temperature_model__mod__alphadisk == enum_hayashi_string) {
      DF_MDOT = 3*pi*AC_nut_global__mod__alphadisk[vertexIdx.x-NGHOST_VAL]*value(Field(AC_isigma__mod__alphadisk-1))
    }
    else if(enum_temperature_model__mod__alphadisk == enum_radiative_string) {
      lgsigma1_0=(AC_c1__mod__alphadisk[vertexIdx.x-NGHOST]-AC_cprime__mod__alphadisk)/2.1
      lgsigma2_0=(AC_c3__mod__alphadisk[vertexIdx.x-NGHOST]-AC_c2__mod__alphadisk[vertexIdx.x-NGHOST])/0.9
      lgsigma_0=alog10(value(Field(AC_isigma__mod__alphadisk-1)))
      if (lgsigma_0<=lgsigma1_0) {
        lgmdot_0=AC_c1__mod__alphadisk[vertexIdx.x-NGHOST] + lgsigma_0
      }
      else if ((lgsigma_0>lgsigma1_0) && (lgsigma_0<lgsigma2_0)) {
        lgmdot_0=AC_c2__mod__alphadisk[vertexIdx.x-NGHOST] + 2.0*lgsigma_0
      }
      else if (lgsigma_0>=lgsigma2_0) {
        lgmdot_0=AC_c3__mod__alphadisk[vertexIdx.x-NGHOST] + 1.1*lgsigma_0
      }
      else {
      }
      DF_MDOT=pow(10,lgmdot_0)
    }
  }
  else if(enum_temperature_model__mod__alphadisk == enum_radiative_string) {
    sig_1=value(Field(AC_isigma__mod__alphadisk-1))
    if (sig_1>AC_maxsigma__mod__alphadisk) {
    }
    else if ((sig_1>=AC_sigma_middle__mod__alphadisk) && (sig_1<=AC_maxsigma__mod__alphadisk)) {
      isig_do_1 = floor((value(Field(AC_isigma__mod__alphadisk-1)) - AC_minsigma__mod__alphadisk)*AC_dsig1__mod__alphadisk) + 1
      isig_up_1 =  isig_do_1+1
      sdo_1 = AC_minsigma__mod__alphadisk + (isig_do_1-1)*AC_dsig__mod__alphadisk
      sup_1 = AC_minsigma__mod__alphadisk + (isig_up_1-1)*AC_dsig__mod__alphadisk
      DF_TMID = AC_dsig1__mod__alphadisk*(AC_tmid1_table__mod__alphadisk[isig_do_1-1][vertexIdx.x]*(sup_1-sig_1)+ AC_tmid1_table__mod__alphadisk[isig_up_1-1][vertexIdx.x]*(sig_1-sdo_1))
    }
    else if ((sig_1>=AC_sigma_floor__mod__alphadisk) && (sig_1<=AC_sigma_middle__mod__alphadisk)) {
      lnsig_1=log(sig_1)
      isig_do_1 = floor((lnsig_1 - AC_minlnsigma__mod__alphadisk)*AC_dlnsig1__mod__alphadisk) + 1
      isig_up_1 =  isig_do_1+1
      lnsdo_1=AC_minlnsigma__mod__alphadisk+(isig_do_1-1)*AC_dlnsig__mod__alphadisk
      lnsup_1=AC_minlnsigma__mod__alphadisk+(isig_up_1-1)*AC_dlnsig__mod__alphadisk
      DF_TMID = AC_dlnsig1__mod__alphadisk*(AC_tmid2_table__mod__alphadisk[isig_do_1-1][vertexIdx.x]*(lnsup_1-lnsig_1)+ AC_tmid2_table__mod__alphadisk[isig_up_1-1][vertexIdx.x]*(lnsig_1-lnsdo_1))
    }
    else {
    }
    if(enum_temperature_model__mod__alphadisk == enum_hayashi_string) {
      DF_MDOT = 3*pi*AC_nut_global__mod__alphadisk[vertexIdx.x-NGHOST_VAL]*value(Field(AC_isigma__mod__alphadisk-1))
    }
    else if(enum_temperature_model__mod__alphadisk == enum_radiative_string) {
      lgsigma1_2=(AC_c1__mod__alphadisk[vertexIdx.x-NGHOST]-AC_cprime__mod__alphadisk)/2.1
      lgsigma2_2=(AC_c3__mod__alphadisk[vertexIdx.x-NGHOST]-AC_c2__mod__alphadisk[vertexIdx.x-NGHOST])/0.9
      lgsigma_2=alog10(value(Field(AC_isigma__mod__alphadisk-1)))
      if (lgsigma_2<=lgsigma1_2) {
        lgmdot_2=AC_c1__mod__alphadisk[vertexIdx.x-NGHOST] + lgsigma_2
      }
      else if ((lgsigma_2>lgsigma1_2) && (lgsigma_2<lgsigma2_2)) {
        lgmdot_2=AC_c2__mod__alphadisk[vertexIdx.x-NGHOST] + 2.0*lgsigma_2
      }
      else if (lgsigma_2>=lgsigma2_2) {
        lgmdot_2=AC_c3__mod__alphadisk[vertexIdx.x-NGHOST] + 1.1*lgsigma_2
      }
      else {
      }
      DF_MDOT=pow(10,lgmdot_2)
    }
  }
  write(F_MDOT,DF_MDOT)
  write(F_TMID,DF_TMID)
}
#else
Kernel after_timestep_alphadisk(){}
#endif
