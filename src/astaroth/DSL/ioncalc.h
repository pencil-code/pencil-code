#if Leos_ionization_MODULE
Kernel ioncalc(){
  real lnrho
  real ss
  real yh
  real lntt
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
  real ac_transformed_pencil_ugss
  real ac_transformed_pencil_ma2
  real3 ac_transformed_pencil_fpres
  real ac_transformed_pencil_uglntt
  real3 ac_transformed_pencil_sglntt
  real ac_transformed_pencil_transprhos
  real ac_transformed_pencil_initss
  real ac_transformed_pencil_initlnrho
  real ac_transformed_pencil_uuadvec_gss
  real ac_transformed_pencil_advec_cs2
  real ac_transformed_pencil_cool_prof
  real ac_transformed_pencil_ss
  real3 ac_transformed_pencil_gss
  real ac_transformed_pencil_ee
  real ac_transformed_pencil_pp
  real ac_transformed_pencil_lntt
  real ac_transformed_pencil_cs2
  real ac_transformed_pencil_cp
  real ac_transformed_pencil_cp1
  real ac_transformed_pencil_cp1tilde
  real3 ac_transformed_pencil_glntt
  real ac_transformed_pencil_tt
  real ac_transformed_pencil_tt1
  real3 ac_transformed_pencil_gtt
  real ac_transformed_pencil_yh
  Matrix ac_transformed_pencil_hss
  Matrix ac_transformed_pencil_hlntt
  real ac_transformed_pencil_del2tt
  real ac_transformed_pencil_del6tt
  real ac_transformed_pencil_del6lntt
  real ac_transformed_pencil_del2ss
  real ac_transformed_pencil_del6ss
  real ac_transformed_pencil_del2lntt
  real ac_transformed_pencil_cv
  real ac_transformed_pencil_cv1
  real3 ac_transformed_pencil_glnmumol
  real ac_transformed_pencil_ppvap
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
  lntt__0=(2.0/3.0)*((ss/AC_ss_ion__mod__equationofstate+(1-yh)*(log(1-yh+AC_epsi__mod__cparam)-AC_lnrho_h__mod__equationofstate)  +yh*(2*log(yh)-AC_lnrho_e__mod__equationofstate-AC_lnrho_h__mod__equationofstate)  +AC_xhe_term__mod__equationofstate)*fractions1_0+lnrho-2.5)
  tt1__0=exp(-lntt__0)
  ff_0=AC_lnrho_e__mod__equationofstate-lnrho+1.5*lntt__0-tt1__0+log(1-yh+AC_epsi__mod__cparam)-2*log(yh)
  dlntt__0=((2.0/3.0)*(-ff_0-tt1__0)-1)*fractions1_0
  dff_0=dlntt__0*(1.5+tt1__0)-1/(1-yh+AC_epsi__mod__cparam)-2/(yh+AC_epsi__mod__cparam)
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
      lntt__0=(2.0/3.0)*((ss/AC_ss_ion__mod__equationofstate+(1-yh)*(log(1-yh+AC_epsi__mod__cparam)-AC_lnrho_h__mod__equationofstate)  +yh*(2*log(yh)-AC_lnrho_e__mod__equationofstate-AC_lnrho_h__mod__equationofstate)  +AC_xhe_term__mod__equationofstate)*fractions1_0+lnrho-2.5)
      tt1__0=exp(-lntt__0)
      ff_0=AC_lnrho_e__mod__equationofstate-lnrho+1.5*lntt__0-tt1__0+log(1-yh+AC_epsi__mod__cparam)-2*log(yh)
      dlntt__0=((2.0/3.0)*(-ff_0-tt1__0)-1)*fractions1_0
      dff_0=dlntt__0*(1.5+tt1__0)-1/(1-yh+AC_epsi__mod__cparam)-2/yh
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
  lntt=(ss/AC_ss_ion__mod__equationofstate+(1-yh)*(log(1-yh+AC_epsi__mod__cparam)-AC_lnrho_h__mod__equationofstate)  +yh*(2*log(yh)-AC_lnrho_e__mod__equationofstate-AC_lnrho_h__mod__equationofstate)+AC_xhe_term__mod__equationofstate)/(1+yh+AC_xhe__mod__equationofstate)
  lntt=(2.0/3.0)*(lntt+lnrho-2.5)+AC_lntt_ion__mod__equationofstate
  write(F_YH,yh)
  write(F_LNTT,lntt)
}
#else
Kernel ioncalc(){
}
#endif
