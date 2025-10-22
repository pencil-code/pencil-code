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
  else if(AC_enum_source_function_type__mod__radiation == AC_enum_b2zw2_string) {
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
