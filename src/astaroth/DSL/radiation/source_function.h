Kernel source_function(int inu){
  real srad__mod__radiation
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
  if(AC_enum_source_function_type__mod__radiation == AC_enum_lte_string__mod__cparam) {
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
  else if(AC_enum_source_function_type__mod__radiation == AC_enum_twozcolored_string__mod__cparam) {
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
  else if(AC_enum_source_function_type__mod__radiation == AC_enum_blob_string__mod__cparam) {
    if (lfirst) {
      srad__mod__radiation=AC_srad_const__mod__radiation+AC_amplsrad__mod__radiation*exp(-((AC_x__mod__cdata[vertexIdx.x]/AC_radius_srad__mod__radiation)*(AC_x__mod__cdata[vertexIdx.x]/AC_radius_srad__mod__radiation)))  *exp(-((AC_y__mod__cdata[vertexIdx.y]/AC_radius_srad__mod__radiation)*(AC_y__mod__cdata[vertexIdx.y]/AC_radius_srad__mod__radiation)))  *exp(-((AC_z__mod__cdata[vertexIdx.z]/AC_radius_srad__mod__radiation)*(AC_z__mod__cdata[vertexIdx.z]/AC_radius_srad__mod__radiation)))
      lfirst=false
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == AC_enum_cos_string__mod__cparam) {
    if (lfirst) {
      srad__mod__radiation=AC_srad_const__mod__radiation+AC_amplsrad__mod__radiation*cos(AC_kx_srad__mod__radiation*AC_x__mod__cdata[vertexIdx.x])  *cos(AC_ky_srad__mod__radiation*AC_y__mod__cdata[vertexIdx.y])  *cos(AC_kz_srad__mod__radiation*AC_z__mod__cdata[vertexIdx.z])
      lfirst=false
    }
  }
  else if(AC_enum_source_function_type__mod__radiation == AC_enum_b2_string__mod__cparam) {
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
  else if(AC_enum_source_function_type__mod__radiation == AC_enum_b2zw2_string__mod__cparam) {
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
  else if(AC_enum_source_function_type__mod__radiation == AC_enum_read_file_string__mod__cparam) {
  }
  else if(AC_enum_source_function_type__mod__radiation == AC_enum_nothing_string__mod__cparam) {
    srad__mod__radiation=0.0
  }
  else {
  }
  write(SRAD,srad__mod__radiation)
}
