#if Leos_temperature_ionization_MODULE 
Kernel source_function(int inu){
  real srad__mod__radiation
  bool lfirst
  int ilntt_table
  real lntt
  int ierr
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
  Matrix aij_3
  real3 aa_3
  real3 bb_3
  real b2_3
  Matrix aij_4
  real3 aa_4
  real3 bb_4
  real b2_4
  Matrix uij_5
  real3 uu_5
  real3 oo_5
  real o2_5
  if(AC_enum_source_function_type__mod__radiation == enum_lte_string) {
    if (AC_lcutoff_opticallythin__mod__radiation) {
      if(mx == nx) {
        lnrho__0=value(Field(AC_ilnrho__mod__cdata-1))
        lntt__0=value(Field(AC_ilntt__mod__cdata-1))
        yh__0=value(Field(AC_iyh__mod__cdata-1))
      }
      else if(mx == mx) {
        lnrho__0=value(Field(AC_ilnrho__mod__cdata-1))
        lntt__0=value(Field(AC_ilntt__mod__cdata-1))
        yh__0=value(Field(AC_iyh__mod__cdata-1))
      }
      else {
      }
      if (present(lntt)) {
        lntt=lntt__0
      }
      if (false) {
        tmp_0 = 2.5 - 1.5*(AC_lntt_ion__mod__equationofstate-lntt__0) - lnrho__0
      }
      if (false) {
        tt1_0 = exp(-lntt__0)
        tmp_0 = 2*lnrho__0-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__0)+AC_tt_ion___mod__equationofstate*tt1_0
        tmpy_0 = yh__0+AC_ymetals__mod__equationofstate
        if (AC_lhminus_opacity_correction__mod__equationofstate) {
          mu1__0 = AC_mu1_0__mod__equationofstate*(1+yh__0+AC_xhe__mod__equationofstate)
          tmpy1_0 = 1./(1+yh__0)
        }
      }
      srad__mod__radiation=AC_arad__mod__radiation*exp(4*lntt)*AC_scalefactor_srad__mod__radiation[inu-1]
      if (AC_z__mod__cdata[AC_n__mod__cdata-1] > AC_z_cutoff__mod__radiation) {
        srad__mod__radiation=srad__mod__radiation*0.5*(1.-tanh((lntt-log(1.0e4))/log(2.0)))
      }
    }
    else {
      if(mx == nx) {
        lnrho__1=value(Field(AC_ilnrho__mod__cdata-1))
        lntt__1=value(Field(AC_ilntt__mod__cdata-1))
        yh__1=value(Field(AC_iyh__mod__cdata-1))
      }
      else if(mx == mx) {
        lnrho__1=value(Field(AC_ilnrho__mod__cdata-1))
        lntt__1=value(Field(AC_ilntt__mod__cdata-1))
        yh__1=value(Field(AC_iyh__mod__cdata-1))
      }
      else {
      }
      if (present(lntt)) {
        lntt=lntt__1
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
      srad__mod__radiation=AC_arad__mod__radiation*exp(4*lntt)*AC_scalefactor_srad__mod__radiation[inu-1]
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_twozcolored_string) {
    if(mx == nx) {
      lnrho__2=value(Field(AC_ilnrho__mod__cdata-1))
      lntt__2=value(Field(AC_ilntt__mod__cdata-1))
      yh__2=value(Field(AC_iyh__mod__cdata-1))
    }
    else if(mx == mx) {
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
    ilntt_table=max(min(int((lntt-AC_lntt_table0__mod__radiation)*AC_dlntt_table__mod__radiation),AC_nlntt_table__mod__radiation-1),1)
    srad__mod__radiation=exp(AC_lnss_table__mod__radiation[ilntt_table-1][inu-1]  +(AC_lnss_table__mod__radiation[1+ilntt_table-1][inu-1]-AC_lnss_table__mod__radiation[ilntt_table-1][inu-1])  *(lntt-AC_lntt_table__mod__radiation[ilntt_table-1])/(AC_lntt_table__mod__radiation[1+ilntt_table-1]-AC_lntt_table__mod__radiation[ilntt_table-1]))/AC_unit_flux__mod__cdata
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_blob_string) {
    if (lfirst) {
      srad__mod__radiation=AC_srad_const__mod__radiation+AC_amplsrad__mod__radiation*exp(-((AC_x__mod__cdata[vertexIdx.x]/AC_radius_srad__mod__radiation)*(AC_x__mod__cdata[vertexIdx.x]/AC_radius_srad__mod__radiation)))  *exp(-((AC_y__mod__cdata[vertexIdx.y]/AC_radius_srad__mod__radiation)*(AC_y__mod__cdata[vertexIdx.y]/AC_radius_srad__mod__radiation)))  *exp(-((AC_z__mod__cdata[vertexIdx.z]/AC_radius_srad__mod__radiation)*(AC_z__mod__cdata[vertexIdx.z]/AC_radius_srad__mod__radiation)))
      lfirst=false
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_cos_string) {
    if (lfirst) {
      srad__mod__radiation=AC_srad_const__mod__radiation+AC_amplsrad__mod__radiation*cos(AC_kx_srad__mod__radiation*AC_x__mod__cdata[vertexIdx.x])  *cos(AC_ky_srad__mod__radiation*AC_y__mod__cdata[vertexIdx.y])  *cos(AC_kz_srad__mod__radiation*AC_z__mod__cdata[vertexIdx.z])
      lfirst=false
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_b2_string) {
    if (AC_iaa__mod__cdata==0) {
    }
    else {
      aa_3=value(F_AVEC)
      aij_3 = gradient_tensor((Field3){Field(AC_iaa__mod__cdata-1), Field(AC_iaa__mod__cdata), Field(AC_iaa__mod__cdata+1)})
      bb_3=curl(aij_3,aa_3)
      b2_3 = dot(bb_3,bb_3)
      srad__mod__radiation=b2_3
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_b2zw2_string) {
    if (inu==1) {
      if (AC_iaa__mod__cdata==0) {
      }
      else {
        aa_4=value(F_AVEC)
        aij_4 = gradient_tensor((Field3){Field(AC_iaa__mod__cdata-1), Field(AC_iaa__mod__cdata), Field(AC_iaa__mod__cdata+1)})
        bb_4=curl(aij_4,aa_4)
        b2_4 = dot(bb_4,bb_4)
        srad__mod__radiation=b2_4
      }
    }
    else if (inu==2) {
      if (AC_iuu__mod__cdata==0) {
      }
      else {
        uu_5=value(F_UVEC)
        uij_5 = gradient_tensor((Field3){Field(AC_iuu__mod__cdata-1), Field(AC_iuu__mod__cdata), Field(AC_iuu__mod__cdata+1)})
        oo_5=curl(uij_5,uu_5)
        o2_5 = dot(oo_5,oo_5)
        srad__mod__radiation=o2_5
      }
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_read_file_string) {
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_nothing_string) {
    srad__mod__radiation=0.0
  }
  else {
  }
  write(SRAD,srad__mod__radiation)
}
#endif

#if Leos_ionization_MODULE
Kernel source_function(int inu){
  real srad__mod__radiation
  real ac_transformed_pencil_ma2
  real ac_transformed_pencil_sgs_heat
  real ac_transformed_pencil_t00
  real3 ac_transformed_pencil_t0i
  real ac_transformed_pencil_tt
  real ac_transformed_pencil_tt1
  real ac_transformed_pencil_tij[6]
  real ac_transformed_pencil_ywater
  real3 ac_transformed_pencil_aa
  real ac_transformed_pencil_acc
  real ac_transformed_pencil_advec_cs2
  real ac_transformed_pencil_advec_uu
  Matrix ac_transformed_pencil_aij
  real ac_transformed_pencil_b2
  real3 ac_transformed_pencil_bb
  real3 ac_transformed_pencil_bbb
  Matrix ac_transformed_pencil_bij
  real3 ac_transformed_pencil_bunit
  real ac_transformed_pencil_cc[1]
  real ac_transformed_pencil_cc1[1]
  real ac_transformed_pencil_chem_conc[1]
  real ac_transformed_pencil_conc_satm
  real ac_transformed_pencil_condensationrate
  real ac_transformed_pencil_cool
  real ac_transformed_pencil_cool_prof
  real ac_transformed_pencil_cp
  real ac_transformed_pencil_cp1
  real ac_transformed_pencil_cp1tilde
  real ac_transformed_pencil_cs2
  real ac_transformed_pencil_csvap2
  real3 ac_transformed_pencil_curlfcont[AC_n_forcing_cont_max__mod__cparam]
  real3 ac_transformed_pencil_curlo
  real ac_transformed_pencil_cv
  real ac_transformed_pencil_cv1
  Matrix ac_transformed_pencil_d2uidxj
  real3 ac_transformed_pencil_del2t
  real ac_transformed_pencil_del2tt
  real3 ac_transformed_pencil_del2a
  real ac_transformed_pencil_del2lntt
  real ac_transformed_pencil_del2lnrho
  real ac_transformed_pencil_del2rho
  real ac_transformed_pencil_del2ss
  real3 ac_transformed_pencil_del2u
  real3 ac_transformed_pencil_del4graddivu
  real3 ac_transformed_pencil_del4u
  real ac_transformed_pencil_del6tt
  real ac_transformed_pencil_del6lntt
  real ac_transformed_pencil_del6lnrho
  real ac_transformed_pencil_del6rho
  real ac_transformed_pencil_del6ss
  real3 ac_transformed_pencil_del6u
  real3 ac_transformed_pencil_del6u_bulk
  real3 ac_transformed_pencil_del6u_strict
  real3 ac_transformed_pencil_der6u
  Matrix ac_transformed_pencil_der6u_res
  real ac_transformed_pencil_diffus_total
  real ac_transformed_pencil_diffus_total2
  real ac_transformed_pencil_diffus_total3
  real ac_transformed_pencil_diva
  real ac_transformed_pencil_divss
  real ac_transformed_pencil_divu
  real ac_transformed_pencil_divu0
  real ac_transformed_pencil_divud[1]
  real ac_transformed_pencil_divun
  real ac_transformed_pencil_divutot
  real ac_transformed_pencil_e2
  real ac_transformed_pencil_ee
  real ac_transformed_pencil_ekin
  real3 ac_transformed_pencil_el
  real ac_transformed_pencil_epot
  real3 ac_transformed_pencil_evr
  real3 ac_transformed_pencil_evth
  real ac_transformed_pencil_ext_force[4]
  real3 ac_transformed_pencil_fcont[AC_n_forcing_cont_max__mod__cparam]
  real ac_transformed_pencil_ff_cond
  real3 ac_transformed_pencil_fpres
  real3 ac_transformed_pencil_fvisc
  real3 ac_transformed_pencil_gtt
  real ac_transformed_pencil_gcc[3][0]
  real3 ac_transformed_pencil_gg
  real3 ac_transformed_pencil_glntt
  real3 ac_transformed_pencil_glnmumol
  real3 ac_transformed_pencil_glnrho
  real ac_transformed_pencil_glnrho2
  real3 ac_transformed_pencil_glnrhodsum
  real3 ac_transformed_pencil_glnrhos
  real3 ac_transformed_pencil_gnu_smag
  real3 ac_transformed_pencil_grad5divu
  real3 ac_transformed_pencil_graddivu
  real3 ac_transformed_pencil_gradnu
  real3 ac_transformed_pencil_grho
  real3 ac_transformed_pencil_grhop
  real3 ac_transformed_pencil_gshock
  real3 ac_transformed_pencil_gshock_perp
  real3 ac_transformed_pencil_gss
  real ac_transformed_pencil_heat
  real ac_transformed_pencil_heatcool
  real ac_transformed_pencil_hless
  Matrix ac_transformed_pencil_hlntt
  Matrix ac_transformed_pencil_hlnrho
  Matrix ac_transformed_pencil_hss
  real ac_transformed_pencil_initlnrho
  real ac_transformed_pencil_initss
  real ac_transformed_pencil_j2
  real3 ac_transformed_pencil_jj
  real3 ac_transformed_pencil_jxbr
  real ac_transformed_pencil_lambda
  real ac_transformed_pencil_latent_heat
  real ac_transformed_pencil_lntt
  real ac_transformed_pencil_lnrho
  real ac_transformed_pencil_lorentz
  real ac_transformed_pencil_lorentz_gamma
  real ac_transformed_pencil_nu
  real ac_transformed_pencil_nu_smag
  real ac_transformed_pencil_nucl_rate
  real ac_transformed_pencil_nucl_rmin
  real ac_transformed_pencil_o2
  real3 ac_transformed_pencil_ogu
  Matrix ac_transformed_pencil_oij
  real3 ac_transformed_pencil_oo
  real ac_transformed_pencil_ou
  real3 ac_transformed_pencil_oxu
  real ac_transformed_pencil_oxu2
  real ac_transformed_pencil_part_heatcap
  real ac_transformed_pencil_peh
  real ac_transformed_pencil_phi_mn
  real ac_transformed_pencil_phix
  real ac_transformed_pencil_phiy
  real ac_transformed_pencil_pomx
  real ac_transformed_pencil_pomy
  real ac_transformed_pencil_pp
  real ac_transformed_pencil_ppvap
  real ac_transformed_pencil_r_mn
  real ac_transformed_pencil_r_mn1
  real ac_transformed_pencil_rcyl_mn
  real ac_transformed_pencil_rcyl_mn1
  real ac_transformed_pencil_rho
  real ac_transformed_pencil_rho1
  real3 ac_transformed_pencil_rho1gpp
  real ac_transformed_pencil_rho_anel
  real ac_transformed_pencil_rhod[1]
  real ac_transformed_pencil_rhodsum
  real ac_transformed_pencil_rhop
  real ac_transformed_pencil_rhos1
  real3 ac_transformed_pencil_rhougu
  real3 ac_transformed_pencil_rr
  Matrix ac_transformed_pencil_sdij[AC_ndustspec__mod__cparam]
  real3 ac_transformed_pencil_sglntt
  real3 ac_transformed_pencil_sglnrho
  real ac_transformed_pencil_shock
  real ac_transformed_pencil_shock_perp
  Matrix ac_transformed_pencil_sij
  real ac_transformed_pencil_sij2
  Matrix ac_transformed_pencil_snij
  real ac_transformed_pencil_ss
  real ac_transformed_pencil_ss12
  real ac_transformed_pencil_ssat
  real ac_transformed_pencil_tauascalar
  real ac_transformed_pencil_totenergy_rel
  real ac_transformed_pencil_transprho
  real ac_transformed_pencil_transprhos
  real3 ac_transformed_pencil_transpurho
  real ac_transformed_pencil_ttc
  Matrix ac_transformed_pencil_u0ij
  real ac_transformed_pencil_u1u23
  real ac_transformed_pencil_u1u32
  real ac_transformed_pencil_u2
  real ac_transformed_pencil_u2u13
  real ac_transformed_pencil_u2u31
  real ac_transformed_pencil_u3u12
  real ac_transformed_pencil_u3u21
  real3 ac_transformed_pencil_udropav
  real ac_transformed_pencil_uglntt
  real ac_transformed_pencil_uglnrho
  real ac_transformed_pencil_ugrho
  real ac_transformed_pencil_ugss
  real3 ac_transformed_pencil_ugu
  real ac_transformed_pencil_ugu2
  Matrix ac_transformed_pencil_uij
  Matrix ac_transformed_pencil_uij5
  real3 ac_transformed_pencil_uij5glnrho
  Matrix ac_transformed_pencil_uij6
  Tensor ac_transformed_pencil_uijk
  Matrix ac_transformed_pencil_utotij
  real3 ac_transformed_pencil_uu
  real3 ac_transformed_pencil_uu0
  real3 ac_transformed_pencil_uu_advec
  real3 ac_transformed_pencil_uu_sph
  real ac_transformed_pencil_uuadvec_glnrho
  real ac_transformed_pencil_uuadvec_grho
  real ac_transformed_pencil_uuadvec_gss
  real3 ac_transformed_pencil_uuadvec_guu
  real3 ac_transformed_pencil_uud[AC_ndustspec__mod__cparam]
  real3 ac_transformed_pencil_uun
  real3 ac_transformed_pencil_uutot
  real3 ac_transformed_pencil_uxb
  real ac_transformed_pencil_va2
  real3 ac_transformed_pencil_velx
  real ac_transformed_pencil_visc_heat
  real ac_transformed_pencil_watermixingratio
  real ac_transformed_pencil_x_mn
  real ac_transformed_pencil_yh
  real ac_transformed_pencil_y_mn
  real ac_transformed_pencil_z_mn
  bool lfirst
  int ilntt_table
  real lntt
  int ierr
  real lnrho__1_3
  real yh__1_3
  real lntt__1_3
  real tt_1_3
  real fractions_1_3
  real exponent_1_3
  real lnrho_0_1_3
  real yh_0_1_3
  real lntt_0_1_3
  real r_0_1_3
  real dlnttdy_0_1_3
  real drdy_0_1_3
  real temp_0_1_3
  real dlnppdlnrho_0_1_3
  real fractions_0_1_3
  real fractions1_0_1_3
  real dlnppdss_0_1_3
  real tt1_0_1_3
  real lnrho__2_3
  real yh__2_3
  real lntt__2_3
  real tt_2_3
  real fractions_2_3
  real exponent_2_3
  real lnrho_0_2_3
  real yh_0_2_3
  real lntt_0_2_3
  real r_0_2_3
  real dlnttdy_0_2_3
  real drdy_0_2_3
  real temp_0_2_3
  real dlnppdlnrho_0_2_3
  real fractions_0_2_3
  real fractions1_0_2_3
  real dlnppdss_0_2_3
  real tt1_0_2_3
  real lnrho__1_4
  real yh__1_4
  real lntt__1_4
  real tt_1_4
  real fractions_1_4
  real exponent_1_4
  real lnrho_0_1_4
  real yh_0_1_4
  real lntt_0_1_4
  real r_0_1_4
  real dlnttdy_0_1_4
  real drdy_0_1_4
  real temp_0_1_4
  real dlnppdlnrho_0_1_4
  real fractions_0_1_4
  real fractions1_0_1_4
  real dlnppdss_0_1_4
  real tt1_0_1_4
  real lnrho__2_4
  real yh__2_4
  real lntt__2_4
  real tt_2_4
  real fractions_2_4
  real exponent_2_4
  real lnrho_0_2_4
  real yh_0_2_4
  real lntt_0_2_4
  real r_0_2_4
  real dlnttdy_0_2_4
  real drdy_0_2_4
  real temp_0_2_4
  real dlnppdlnrho_0_2_4
  real fractions_0_2_4
  real fractions1_0_2_4
  real dlnppdss_0_2_4
  real tt1_0_2_4
  real lnrho__1_5
  real yh__1_5
  real lntt__1_5
  real tt_1_5
  real fractions_1_5
  real exponent_1_5
  real lnrho_0_1_5
  real yh_0_1_5
  real lntt_0_1_5
  real r_0_1_5
  real dlnttdy_0_1_5
  real drdy_0_1_5
  real temp_0_1_5
  real dlnppdlnrho_0_1_5
  real fractions_0_1_5
  real fractions1_0_1_5
  real dlnppdss_0_1_5
  real tt1_0_1_5
  real lnrho__2_5
  real yh__2_5
  real lntt__2_5
  real tt_2_5
  real fractions_2_5
  real exponent_2_5
  real lnrho_0_2_5
  real yh_0_2_5
  real lntt_0_2_5
  real r_0_2_5
  real dlnttdy_0_2_5
  real drdy_0_2_5
  real temp_0_2_5
  real dlnppdlnrho_0_2_5
  real fractions_0_2_5
  real fractions1_0_2_5
  real dlnppdss_0_2_5
  real tt1_0_2_5
  Matrix aij_6
  real3 aa_6
  real3 bb_6
  real b2_6
  Matrix aij_7
  real3 aa_7
  real3 bb_7
  real b2_7
  Matrix uij_8
  real3 uu_8
  real3 oo_8
  real o2_8
  if(AC_enum_source_function_type__mod__radiation == enum_lte_string) {
    if (AC_lcutoff_opticallythin__mod__radiation) {
      if(mx == nx) {
        lnrho__1_3=value(Field(AC_ilnrho__mod__cdata-1))
        yh__1_3   =value(Field(AC_iyh__mod__cdata-1))
        lntt__1_3 =value(Field(AC_ilntt__mod__cdata-1))
        if (present(lntt)) {
          lntt=lntt__1_3
        }
        if (false) {
          lnrho_0_1_3=value(Field(AC_ilnrho__mod__cdata-1))
          yh_0_1_3=value(Field(AC_iyh__mod__cdata-1))
          lntt_0_1_3=value(Field(AC_ilntt__mod__cdata-1))
          tt1_0_1_3=exp(-lntt_0_1_3)
          fractions_0_1_3=1+yh_0_1_3+AC_xhe__mod__equationofstate
          fractions1_0_1_3=1/fractions_0_1_3
          r_0_1_3=AC_lnrho_e__mod__equationofstate-lnrho_0_1_3+1.5*(lntt_0_1_3-AC_lntt_ion__mod__equationofstate)-AC_tt_ion__mod__equationofstate*tt1_0_1_3+log(1-yh_0_1_3+epsi)-2*log(yh_0_1_3)
          dlnttdy_0_1_3=(2*(-r_0_1_3-AC_tt_ion__mod__equationofstate*tt1_0_1_3)-3)/3*fractions1_0_1_3
          drdy_0_1_3=dlnttdy_0_1_3*(1.5+AC_tt_ion__mod__equationofstate*tt1_0_1_3)-1/(1-yh_0_1_3+epsi)-2/yh_0_1_3
          temp_0_1_3=(dlnttdy_0_1_3+fractions1_0_1_3)/drdy_0_1_3
          dlnppdlnrho_0_1_3=(5-2*AC_tt_ion__mod__equationofstate*tt1_0_1_3*temp_0_1_3)/3
          dlnppdss_0_1_3=AC_ss_ion1__mod__equationofstate*fractions1_0_1_3*(dlnppdlnrho_0_1_3-temp_0_1_3-1)
        }
        if (false || false || false) {
          tt_1_3=exp(lntt__1_3)
        }
        if (false || false) {
          fractions_1_3=1+yh__1_3+AC_xhe__mod__equationofstate
        }
        if (false) {
          exponent_1_3 = (2*lnrho__1_3-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__1_3)+AC_tt_ion___mod__equationofstate/tt_1_3)
          exponent_1_3 = min(exponent_1_3,log(huge1)-5.)
        }
      }
      else if(mx == mx)   {
        lnrho__2_3=value(Field(AC_ilnrho__mod__cdata-1))
        yh__2_3   =value(Field(AC_iyh__mod__cdata-1))
        lntt__2_3 =value(Field(AC_ilntt__mod__cdata-1))
        if (present(lntt)) {
          lntt=lntt__2_3
        }
        if (false) {
          lnrho_0_2_3=value(Field(AC_ilnrho__mod__cdata-1))
          yh_0_2_3=value(Field(AC_iyh__mod__cdata-1))
          lntt_0_2_3=value(Field(AC_ilntt__mod__cdata-1))
          tt1_0_2_3=exp(-lntt_0_2_3)
          fractions_0_2_3=1+yh_0_2_3+AC_xhe__mod__equationofstate
          fractions1_0_2_3=1/fractions_0_2_3
          r_0_2_3=AC_lnrho_e__mod__equationofstate-lnrho_0_2_3+1.5*(lntt_0_2_3-AC_lntt_ion__mod__equationofstate)-AC_tt_ion__mod__equationofstate*tt1_0_2_3+log(1-yh_0_2_3+epsi)-2*log(yh_0_2_3)
          dlnttdy_0_2_3=(2*(-r_0_2_3-AC_tt_ion__mod__equationofstate*tt1_0_2_3)-3)/3*fractions1_0_2_3
          drdy_0_2_3=dlnttdy_0_2_3*(1.5+AC_tt_ion__mod__equationofstate*tt1_0_2_3)-1/(1-yh_0_2_3+epsi)-2/yh_0_2_3
          temp_0_2_3=(dlnttdy_0_2_3+fractions1_0_2_3)/drdy_0_2_3
          dlnppdlnrho_0_2_3=(5-2*AC_tt_ion__mod__equationofstate*tt1_0_2_3*temp_0_2_3)/3
          dlnppdss_0_2_3=AC_ss_ion1__mod__equationofstate*fractions1_0_2_3*(dlnppdlnrho_0_2_3-temp_0_2_3-1)
        }
        if (false || false || false) {
          tt_2_3=exp(lntt__2_3)
        }
        if (false || false) {
          fractions_2_3=1+yh__2_3+AC_xhe__mod__equationofstate
        }
        if (false) {
          exponent_2_3 = (2*lnrho__2_3-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__2_3)+AC_tt_ion___mod__equationofstate/tt_2_3)
          exponent_2_3 = min(exponent_2_3,log(huge1)-5.)
        }
      }
      else {
      }
      srad__mod__radiation=AC_arad__mod__radiation*exp(4*lntt)*AC_scalefactor_srad__mod__radiation[inu-1]
      if (AC_z__mod__cdata[AC_n__mod__cdata-1] > AC_z_cutoff__mod__radiation) {
        srad__mod__radiation=srad__mod__radiation*0.5*(1.-tanh((lntt-log(1.0e4))/log(2.0)))
      }
    }
    else {
      if(mx == nx) {
        lnrho__1_4=value(Field(AC_ilnrho__mod__cdata-1))
        yh__1_4   =value(Field(AC_iyh__mod__cdata-1))
        lntt__1_4 =value(Field(AC_ilntt__mod__cdata-1))
        if (present(lntt)) {
          lntt=lntt__1_4
        }
        if (false) {
          lnrho_0_1_4=value(Field(AC_ilnrho__mod__cdata-1))
          yh_0_1_4=value(Field(AC_iyh__mod__cdata-1))
          lntt_0_1_4=value(Field(AC_ilntt__mod__cdata-1))
          tt1_0_1_4=exp(-lntt_0_1_4)
          fractions_0_1_4=1+yh_0_1_4+AC_xhe__mod__equationofstate
          fractions1_0_1_4=1/fractions_0_1_4
          r_0_1_4=AC_lnrho_e__mod__equationofstate-lnrho_0_1_4+1.5*(lntt_0_1_4-AC_lntt_ion__mod__equationofstate)-AC_tt_ion__mod__equationofstate*tt1_0_1_4+log(1-yh_0_1_4+epsi)-2*log(yh_0_1_4)
          dlnttdy_0_1_4=(2*(-r_0_1_4-AC_tt_ion__mod__equationofstate*tt1_0_1_4)-3)/3*fractions1_0_1_4
          drdy_0_1_4=dlnttdy_0_1_4*(1.5+AC_tt_ion__mod__equationofstate*tt1_0_1_4)-1/(1-yh_0_1_4+epsi)-2/yh_0_1_4
          temp_0_1_4=(dlnttdy_0_1_4+fractions1_0_1_4)/drdy_0_1_4
          dlnppdlnrho_0_1_4=(5-2*AC_tt_ion__mod__equationofstate*tt1_0_1_4*temp_0_1_4)/3
          dlnppdss_0_1_4=AC_ss_ion1__mod__equationofstate*fractions1_0_1_4*(dlnppdlnrho_0_1_4-temp_0_1_4-1)
        }
        if (false || false || false) {
          tt_1_4=exp(lntt__1_4)
        }
        if (false || false) {
          fractions_1_4=1+yh__1_4+AC_xhe__mod__equationofstate
        }
        if (false) {
          exponent_1_4 = (2*lnrho__1_4-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__1_4)+AC_tt_ion___mod__equationofstate/tt_1_4)
          exponent_1_4 = min(exponent_1_4,log(huge1)-5.)
        }
      }
      else if(mx == mx)   {
        lnrho__2_4=value(Field(AC_ilnrho__mod__cdata-1))
        yh__2_4   =value(Field(AC_iyh__mod__cdata-1))
        lntt__2_4 =value(Field(AC_ilntt__mod__cdata-1))
        if (present(lntt)) {
          lntt=lntt__2_4
        }
        if (false) {
          lnrho_0_2_4=value(Field(AC_ilnrho__mod__cdata-1))
          yh_0_2_4=value(Field(AC_iyh__mod__cdata-1))
          lntt_0_2_4=value(Field(AC_ilntt__mod__cdata-1))
          tt1_0_2_4=exp(-lntt_0_2_4)
          fractions_0_2_4=1+yh_0_2_4+AC_xhe__mod__equationofstate
          fractions1_0_2_4=1/fractions_0_2_4
          r_0_2_4=AC_lnrho_e__mod__equationofstate-lnrho_0_2_4+1.5*(lntt_0_2_4-AC_lntt_ion__mod__equationofstate)-AC_tt_ion__mod__equationofstate*tt1_0_2_4+log(1-yh_0_2_4+epsi)-2*log(yh_0_2_4)
          dlnttdy_0_2_4=(2*(-r_0_2_4-AC_tt_ion__mod__equationofstate*tt1_0_2_4)-3)/3*fractions1_0_2_4
          drdy_0_2_4=dlnttdy_0_2_4*(1.5+AC_tt_ion__mod__equationofstate*tt1_0_2_4)-1/(1-yh_0_2_4+epsi)-2/yh_0_2_4
          temp_0_2_4=(dlnttdy_0_2_4+fractions1_0_2_4)/drdy_0_2_4
          dlnppdlnrho_0_2_4=(5-2*AC_tt_ion__mod__equationofstate*tt1_0_2_4*temp_0_2_4)/3
          dlnppdss_0_2_4=AC_ss_ion1__mod__equationofstate*fractions1_0_2_4*(dlnppdlnrho_0_2_4-temp_0_2_4-1)
        }
        if (false || false || false) {
          tt_2_4=exp(lntt__2_4)
        }
        if (false || false) {
          fractions_2_4=1+yh__2_4+AC_xhe__mod__equationofstate
        }
        if (false) {
          exponent_2_4 = (2*lnrho__2_4-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__2_4)+AC_tt_ion___mod__equationofstate/tt_2_4)
          exponent_2_4 = min(exponent_2_4,log(huge1)-5.)
        }
      }
      else {
      }
      srad__mod__radiation=AC_arad__mod__radiation*exp(4*lntt)*AC_scalefactor_srad__mod__radiation[inu-1]
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_twozcolored_string)   {
    if(mx == nx) {
      lnrho__1_5=value(Field(AC_ilnrho__mod__cdata-1))
      yh__1_5   =value(Field(AC_iyh__mod__cdata-1))
      lntt__1_5 =value(Field(AC_ilntt__mod__cdata-1))
      if (present(lntt)) {
        lntt=lntt__1_5
      }
      if (false) {
        lnrho_0_1_5=value(Field(AC_ilnrho__mod__cdata-1))
        yh_0_1_5=value(Field(AC_iyh__mod__cdata-1))
        lntt_0_1_5=value(Field(AC_ilntt__mod__cdata-1))
        tt1_0_1_5=exp(-lntt_0_1_5)
        fractions_0_1_5=1+yh_0_1_5+AC_xhe__mod__equationofstate
        fractions1_0_1_5=1/fractions_0_1_5
        r_0_1_5=AC_lnrho_e__mod__equationofstate-lnrho_0_1_5+1.5*(lntt_0_1_5-AC_lntt_ion__mod__equationofstate)-AC_tt_ion__mod__equationofstate*tt1_0_1_5+log(1-yh_0_1_5+epsi)-2*log(yh_0_1_5)
        dlnttdy_0_1_5=(2*(-r_0_1_5-AC_tt_ion__mod__equationofstate*tt1_0_1_5)-3)/3*fractions1_0_1_5
        drdy_0_1_5=dlnttdy_0_1_5*(1.5+AC_tt_ion__mod__equationofstate*tt1_0_1_5)-1/(1-yh_0_1_5+epsi)-2/yh_0_1_5
        temp_0_1_5=(dlnttdy_0_1_5+fractions1_0_1_5)/drdy_0_1_5
        dlnppdlnrho_0_1_5=(5-2*AC_tt_ion__mod__equationofstate*tt1_0_1_5*temp_0_1_5)/3
        dlnppdss_0_1_5=AC_ss_ion1__mod__equationofstate*fractions1_0_1_5*(dlnppdlnrho_0_1_5-temp_0_1_5-1)
      }
      if (false || false || false) {
        tt_1_5=exp(lntt__1_5)
      }
      if (false || false) {
        fractions_1_5=1+yh__1_5+AC_xhe__mod__equationofstate
      }
      if (false) {
        exponent_1_5 = (2*lnrho__1_5-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__1_5)+AC_tt_ion___mod__equationofstate/tt_1_5)
        exponent_1_5 = min(exponent_1_5,log(huge1)-5.)
      }
    }
    else if(mx == mx)   {
      lnrho__2_5=value(Field(AC_ilnrho__mod__cdata-1))
      yh__2_5   =value(Field(AC_iyh__mod__cdata-1))
      lntt__2_5 =value(Field(AC_ilntt__mod__cdata-1))
      if (present(lntt)) {
        lntt=lntt__2_5
      }
      if (false) {
        lnrho_0_2_5=value(Field(AC_ilnrho__mod__cdata-1))
        yh_0_2_5=value(Field(AC_iyh__mod__cdata-1))
        lntt_0_2_5=value(Field(AC_ilntt__mod__cdata-1))
        tt1_0_2_5=exp(-lntt_0_2_5)
        fractions_0_2_5=1+yh_0_2_5+AC_xhe__mod__equationofstate
        fractions1_0_2_5=1/fractions_0_2_5
        r_0_2_5=AC_lnrho_e__mod__equationofstate-lnrho_0_2_5+1.5*(lntt_0_2_5-AC_lntt_ion__mod__equationofstate)-AC_tt_ion__mod__equationofstate*tt1_0_2_5+log(1-yh_0_2_5+epsi)-2*log(yh_0_2_5)
        dlnttdy_0_2_5=(2*(-r_0_2_5-AC_tt_ion__mod__equationofstate*tt1_0_2_5)-3)/3*fractions1_0_2_5
        drdy_0_2_5=dlnttdy_0_2_5*(1.5+AC_tt_ion__mod__equationofstate*tt1_0_2_5)-1/(1-yh_0_2_5+epsi)-2/yh_0_2_5
        temp_0_2_5=(dlnttdy_0_2_5+fractions1_0_2_5)/drdy_0_2_5
        dlnppdlnrho_0_2_5=(5-2*AC_tt_ion__mod__equationofstate*tt1_0_2_5*temp_0_2_5)/3
        dlnppdss_0_2_5=AC_ss_ion1__mod__equationofstate*fractions1_0_2_5*(dlnppdlnrho_0_2_5-temp_0_2_5-1)
      }
      if (false || false || false) {
        tt_2_5=exp(lntt__2_5)
      }
      if (false || false) {
        fractions_2_5=1+yh__2_5+AC_xhe__mod__equationofstate
      }
      if (false) {
        exponent_2_5 = (2*lnrho__2_5-AC_lnrho_e___mod__equationofstate+1.5*(AC_lntt_ion___mod__equationofstate-lntt__2_5)+AC_tt_ion___mod__equationofstate/tt_2_5)
        exponent_2_5 = min(exponent_2_5,log(huge1)-5.)
      }
    }
    else {
    }
    ilntt_table=max(min(int((lntt-AC_lntt_table0__mod__radiation)*AC_dlntt_table__mod__radiation),AC_nlntt_table__mod__radiation-1),1)
    srad__mod__radiation=exp(AC_lnss_table__mod__radiation[ilntt_table-1][inu-1]  +(AC_lnss_table__mod__radiation[1+ilntt_table-1][inu-1]-AC_lnss_table__mod__radiation[ilntt_table-1][inu-1])  *(lntt-AC_lntt_table__mod__radiation[ilntt_table-1])/(AC_lntt_table__mod__radiation[1+ilntt_table-1]-AC_lntt_table__mod__radiation[ilntt_table-1]))/AC_unit_flux__mod__cdata
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_blob_string)   {
    if (lfirst) {
      srad__mod__radiation=AC_srad_const__mod__radiation+AC_amplsrad__mod__radiation*exp(-((AC_x__mod__cdata[vertexIdx.x]/AC_radius_srad__mod__radiation)*(AC_x__mod__cdata[vertexIdx.x]/AC_radius_srad__mod__radiation)))  *exp(-((AC_y__mod__cdata[vertexIdx.y]/AC_radius_srad__mod__radiation)*(AC_y__mod__cdata[vertexIdx.y]/AC_radius_srad__mod__radiation)))  *exp(-((AC_z__mod__cdata[vertexIdx.z]/AC_radius_srad__mod__radiation)*(AC_z__mod__cdata[vertexIdx.z]/AC_radius_srad__mod__radiation)))
      lfirst=false
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_cos_string)   {
    if (lfirst) {
      srad__mod__radiation=AC_srad_const__mod__radiation+AC_amplsrad__mod__radiation*cos(AC_kx_srad__mod__radiation*AC_x__mod__cdata[vertexIdx.x])  *cos(AC_ky_srad__mod__radiation*AC_y__mod__cdata[vertexIdx.y])  *cos(AC_kz_srad__mod__radiation*AC_z__mod__cdata[vertexIdx.z])
      lfirst=false
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_b2_string)   {
    if (AC_iaa__mod__cdata==0) {
    }
    else {
      aa_6=value(F_AVEC)
      aij_6 = gradient_tensor((Field3){Field(AC_iaa__mod__cdata-1), Field(AC_iaa__mod__cdata), Field(AC_iaa__mod__cdata+1)})
      bb_6=curl(aij_6,aa_6)
      b2_6 = dot(bb_6,bb_6)
      srad__mod__radiation=b2_6
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_b2zw2_string)   {
    if (inu==1) {
      if (AC_iaa__mod__cdata==0) {
      }
      else {
        aa_7=value(F_AVEC)
        aij_7 = gradient_tensor((Field3){Field(AC_iaa__mod__cdata-1), Field(AC_iaa__mod__cdata), Field(AC_iaa__mod__cdata+1)})
        bb_7=curl(aij_7,aa_7)
        b2_7 = dot(bb_7,bb_7)
        srad__mod__radiation=b2_7
      }
    }
    else if (inu==2) {
      if (AC_iuu__mod__cdata==0) {
      }
      else {
        uu_8=value(F_UVEC)
        uij_8 = gradient_tensor((Field3){Field(AC_iuu__mod__cdata-1), Field(AC_iuu__mod__cdata), Field(AC_iuu__mod__cdata+1)})
        oo_8=curl(uij_8,uu_8)
        o2_8 = dot(oo_8,oo_8)
        srad__mod__radiation=o2_8
      }
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_read_file_string)   {
  }
  else if(AC_enum_source_function_type__mod__radiation == enum_nothing_string)   {
    srad__mod__radiation=0.0
  }
  else {
  }

  write(SRAD,srad__mod__radiation)
}
#endif
