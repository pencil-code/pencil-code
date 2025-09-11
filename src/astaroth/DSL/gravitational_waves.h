#if LGRAVITATIONAL_WAVES_HTXK
Kernel gravitational_waves_solve_and_stress(real AC_t__mod__cdata, real AC_dt__mod__cdata){
    real pij[6]
    real kij[6]
    real e_t[6]
    real e_x[6]
    real sij_re[6]
    real sij_im[6]
    real delij[6]
    real3 e1
    real3 e2
    real3 kvec
    int i
    int j
    int p
    int q
    int ik
    int stat
    int ij
    int pq
    int ip
    int jq
    int jstress_ij
    real fact
    real delkt
    real om2_min
    real kmin
    real ksqr
    real one_over_k2
    real k1
    real k2
    real k3
    real k1sqr
    real k2sqr
    real k3sqr
    real ksqrt
    real hhtre
    real hhtim
    real hhxre
    real hhxim
    real coefare
    real coefaim
    real ggtre
    real ggtim
    real ggxre
    real ggxim
    real coefbre
    real coefbim
    real e_ij_t
    real e_ij_x
    real cosot
    real sinot
    real sinot_minus
    real om12
    real om
    real om1
    real om2
    real dt1
    real ett
    real etx
    real ext
    real exx
    real discrim2
    real om_rat_lam
    real om_rat_mat
    real om_rat_matt
    real om_rat_tot1
    real ds_t_re
    real ds_t_im
    real ds_x_re
    real ds_x_im
    complex coefa
    complex coefb
    complex om_cmplx
    complex hcomplex_new
    complex gcomplex_new
    complex discrim
    complex det1
    complex lam1
    complex lam2
    complex explam1t
    complex explam2t
    complex cosoth
    complex cosotg
    complex sinoth
    complex sinotg 
    bool lsign_om2
    real DF_HHT = value(F_HHT)
    real DF_HHX = value(F_HHX)
    real DF_HHTIM = value(F_HHTIM)
    real DF_HHXIM = value(F_HHXIM)
    real DF_GGT = value(F_GGT)
    real DF_GGX = value(F_GGX)
    real DF_GGTIM = value(F_GGTIM)
    real DF_GGXIM = value(F_GGXIM)


    delkt=AC_delk__mod__gravitational_waves_htxk
    if (AC_ldelkt__mod__gravitational_waves_htxk) {
      if(AC_enum_idelkt__mod__gravitational_waves_htxk == AC_enum_jump_string__mod__cparam) {
        if (AC_t__mod__cdata>AC_tdelk__mod__gravitational_waves_htxk) {
          delkt=0.
        }
      }
      else if(AC_enum_idelkt__mod__gravitational_waves_htxk == AC_enum_exponential_string__mod__cparam) {
        if (AC_t__mod__cdata>AC_tdelk__mod__gravitational_waves_htxk) {
          delkt=exp(-(AC_t__mod__cdata-AC_tdelk__mod__gravitational_waves_htxk)/AC_tau_delk__mod__gravitational_waves_htxk)
        }
      }
      else {
      }
    }
    real scale_factor__mod__gravitational_waves_htxk
    if (AC_lgpu__mod__cparam) {
      if (AC_lread_scl_factor_file__mod__cdata) {
        lgt_current_0=alog10(AC_t__mod__cdata)+AC_lgt_ini__mod__gravitational_waves_htxk
        int it_file_0=int((lgt_current_0-AC_lgt0__mod__gravitational_waves_htxk)/AC_dlgt__mod__gravitational_waves_htxk)+1
        lgt1_0=AC_lgt_file__mod__gravitational_waves_htxk[it_file_0-1]
        lgt2_0=AC_lgt_file__mod__gravitational_waves_htxk[1+it_file_0-1]
        lgf1_0=AC_lgff__mod__gravitational_waves_htxk[it_file_0-1]
        lgf2_0=AC_lgff__mod__gravitational_waves_htxk[1+it_file_0-1]
        lgf_0=lgf1_0+(lgt_current_0-lgt1_0)*(lgf2_0-lgf1_0)/(lgt2_0-lgt1_0)
        scl_factor_target__mod__cdata=pow(10,lgf_0)/AC_a_ini__mod__gravitational_waves_htxk
        scale_factor__mod__gravitational_waves_htxk=pow(10,lgf_0)/AC_a_ini__mod__gravitational_waves_htxk
      }
      else {
        if (AC_lreheating_gw__mod__gravitational_waves_htxk) {
          scale_factor__mod__gravitational_waves_htxk=0.25*((AC_t__mod__cdata+1.)*(AC_t__mod__cdata+1.))
        }
        else if (AC_lmatter_gw__mod__gravitational_waves_htxk) {
          scale_factor__mod__gravitational_waves_htxk=(AC_t__mod__cdata*AC_t__mod__cdata)/AC_t_equality__mod__gravitational_waves_htxk
        }
        else if (AC_ldark_energy_gw__mod__gravitational_waves_htxk) {
          scale_factor__mod__gravitational_waves_htxk=(AC_t_acceleration__mod__gravitational_waves_htxk*AC_t_acceleration__mod__gravitational_waves_htxk*AC_t_acceleration__mod__gravitational_waves_htxk)/(AC_t__mod__cdata*AC_t_equality__mod__gravitational_waves_htxk)
        }
        else if (AC_lscalar__mod__gravitational_waves_htxk) {
          scale_factor__mod__gravitational_waves_htxk=exp(AC_f_ode__mod__cdata[AC_ilna__mod__gravitational_waves_htxk-1])
        }
        else {
          if (AC_t__mod__cdata+AC_tshift__mod__gravitational_waves_htxk==0.) {
            scale_factor__mod__gravitational_waves_htxk=1.
          }
          else {
            scale_factor__mod__gravitational_waves_htxk=pow((AC_t__mod__cdata+AC_tshift__mod__gravitational_waves_htxk),AC_nscale_factor_conformal__mod__gravitational_waves_htxk)
          }
        }
      }
    }
    real hp_target__mod__cdata
    real appa_target__mod__cdata
    if (AC_lread_scl_factor_file__mod__cdata) {
      lgt_current_2=alog10(AC_t__mod__cdata)+AC_lgt_ini__mod__gravitational_waves_htxk
      int it_file_2=int((lgt_current_2-AC_lgt0__mod__gravitational_waves_htxk)/AC_dlgt__mod__gravitational_waves_htxk)+1
      lgt1_2=AC_lgt_file__mod__gravitational_waves_htxk[it_file_2-1]
      lgt2_2=AC_lgt_file__mod__gravitational_waves_htxk[1+it_file_2-1]
      lgf1_2=AC_lgff2__mod__gravitational_waves_htxk[it_file_2-1]
      lgf2_2=AC_lgff2__mod__gravitational_waves_htxk[1+it_file_2-1]
      lgf_2=lgf1_2+(lgt_current_2-lgt1_2)*(lgf2_2-lgf1_2)/(lgt2_2-lgt1_2)
      hp_target__mod__cdata=pow(10,lgf_2)/AC_hp_ini__mod__gravitational_waves_htxk
      lgf1_2=AC_lgff3__mod__gravitational_waves_htxk[it_file_2-1]
      lgf2_2=AC_lgff3__mod__gravitational_waves_htxk[1+it_file_2-1]
      lgf_2=lgf1_2+(lgt_current_2-lgt1_2)*(lgf2_2-lgf1_2)/(lgt2_2-lgt1_2)
      appa_target__mod__cdata=pow(10,lgf_2)/(AC_hp_ini__mod__gravitational_waves_htxk*AC_hp_ini__mod__gravitational_waves_htxk)
    }
    real horndeski_alpt_eff__mod__gravitational_waves_htxk
    real horndeski_alpm_eff__mod__gravitational_waves_htxk
    real horndeski_alpm_eff2__mod__gravitational_waves_htxk
    real horndeski_alpm_eff3__mod__gravitational_waves_htxk
    if (AC_lhorndeski__mod__gravitational_waves_htxk || AC_lhorndeski_xi__mod__gravitational_waves_htxk) {
      if(AC_enum_ihorndeski_time__mod__gravitational_waves_htxk == AC_enum_const_string__mod__cparam) {
        horndeski_alpt_eff__mod__gravitational_waves_htxk=AC_horndeski_alpt__mod__gravitational_waves_htxk
        horndeski_alpm_eff__mod__gravitational_waves_htxk=AC_horndeski_alpm__mod__gravitational_waves_htxk
      }
      else if(AC_enum_ihorndeski_time__mod__gravitational_waves_htxk == enum_tanh_string) {
        horndeski_alpt_eff__mod__gravitational_waves_htxk=AC_horndeski_alpt__mod__gravitational_waves_htxk*tanh(1.-pow((scale_factor__mod__gravitational_waves_htxk/AC_scale_factor0__mod__gravitational_waves_htxk),AC_horndeski_alpt_exp__mod__gravitational_waves_htxk))
      }
      else if(AC_enum_ihorndeski_time__mod__gravitational_waves_htxk == enum_exp_string) {
        horndeski_alpt_eff__mod__gravitational_waves_htxk=AC_horndeski_alpt__mod__gravitational_waves_htxk*exp(-pow((scale_factor__mod__gravitational_waves_htxk/AC_scale_factor0__mod__gravitational_waves_htxk),AC_horndeski_alpt_exp__mod__gravitational_waves_htxk))
      }
      else if(AC_enum_ihorndeski_time__mod__gravitational_waves_htxk == enum_scale_factor_power_string) {
        horndeski_alpt_eff__mod__gravitational_waves_htxk=AC_horndeski_alpt__mod__gravitational_waves_htxk
        horndeski_alpm_eff__mod__gravitational_waves_htxk=AC_horndeski_alpm__mod__gravitational_waves_htxk*pow((scale_factor__mod__gravitational_waves_htxk*AC_a_ini__mod__gravitational_waves_htxk/AC_scale_factor0__mod__gravitational_waves_htxk),AC_horndeski_alpm_exp__mod__gravitational_waves_htxk)
      }
      else if(AC_enum_ihorndeski_time__mod__gravitational_waves_htxk == enum_matter_string) {
        horndeski_alpt_eff__mod__gravitational_waves_htxk=AC_horndeski_alpt__mod__gravitational_waves_htxk
        if (AC_lread_scl_factor_file__mod__cdata && AC_lread_scl_factor_file_exists__mod__gravitational_waves_htxk) {
          om_rat_matt=pow((scale_factor__mod__gravitational_waves_htxk*AC_a_ini__mod__gravitational_waves_htxk/AC_scale_factor0__mod__gravitational_waves_htxk),(-3))*AC_omm0__mod__gravitational_waves_htxk
          om_rat_tot1=((AC_a_ini__mod__gravitational_waves_htxk*AC_h0__mod__gravitational_waves_htxk*scale_factor__mod__gravitational_waves_htxk/hp_target__mod__cdata/AC_hp_ini__mod__gravitational_waves_htxk)*(AC_a_ini__mod__gravitational_waves_htxk*AC_h0__mod__gravitational_waves_htxk*scale_factor__mod__gravitational_waves_htxk/hp_target__mod__cdata/AC_hp_ini__mod__gravitational_waves_htxk))
          horndeski_alpm_eff__mod__gravitational_waves_htxk=AC_horndeski_alpm__mod__gravitational_waves_htxk*(1-om_rat_matt*om_rat_tot1)/(1-AC_omm0__mod__gravitational_waves_htxk)
        }
        else {
        }
      }
      else if(AC_enum_ihorndeski_time__mod__gravitational_waves_htxk == enum_dark_energy_string) {
        horndeski_alpt_eff__mod__gravitational_waves_htxk=AC_horndeski_alpt__mod__gravitational_waves_htxk
        if (AC_lread_scl_factor_file__mod__cdata && AC_lread_scl_factor_file_exists__mod__gravitational_waves_htxk) {
          om_rat_tot1=((AC_a_ini__mod__gravitational_waves_htxk*AC_h0__mod__gravitational_waves_htxk*scale_factor__mod__gravitational_waves_htxk/hp_target__mod__cdata/AC_hp_ini__mod__gravitational_waves_htxk)*(AC_a_ini__mod__gravitational_waves_htxk*AC_h0__mod__gravitational_waves_htxk*scale_factor__mod__gravitational_waves_htxk/hp_target__mod__cdata/AC_hp_ini__mod__gravitational_waves_htxk))
          horndeski_alpm_eff__mod__gravitational_waves_htxk=AC_horndeski_alpm__mod__gravitational_waves_htxk*om_rat_tot1
        }
        else {
        }
      }
      else {
      }
      if (AC_lread_scl_factor_file__mod__cdata && AC_lread_scl_factor_file_exists__mod__gravitational_waves_htxk) {
        if (AC_lhorndeski__mod__gravitational_waves_htxk) {
          horndeski_alpm_eff__mod__gravitational_waves_htxk=horndeski_alpm_eff__mod__gravitational_waves_htxk*hp_target__mod__cdata
          horndeski_alpm_eff2__mod__gravitational_waves_htxk=horndeski_alpm_eff__mod__gravitational_waves_htxk*hp_target__mod__cdata
        }
        else {
          horndeski_alpm_eff2__mod__gravitational_waves_htxk=(1+0.5*horndeski_alpm_eff__mod__gravitational_waves_htxk)*(hp_target__mod__cdata*hp_target__mod__cdata)
          horndeski_alpm_eff2__mod__gravitational_waves_htxk=horndeski_alpm_eff2__mod__gravitational_waves_htxk*0.5*horndeski_alpm_eff__mod__gravitational_waves_htxk
          horndeski_alpm_eff3__mod__gravitational_waves_htxk=0.5*AC_horndeski_alpm_prime__mod__gravitational_waves_htxk*hp_target__mod__cdata
          horndeski_alpm_eff__mod__gravitational_waves_htxk=1.+0.5*horndeski_alpm_eff__mod__gravitational_waves_htxk
        }
      }
      else {
        if (AC_lhorndeski__mod__gravitational_waves_htxk) {
          horndeski_alpm_eff__mod__gravitational_waves_htxk=horndeski_alpm_eff__mod__gravitational_waves_htxk/scale_factor__mod__gravitational_waves_htxk
          horndeski_alpm_eff2__mod__gravitational_waves_htxk=horndeski_alpm_eff__mod__gravitational_waves_htxk/scale_factor__mod__gravitational_waves_htxk
        }
        else {
          horndeski_alpm_eff2__mod__gravitational_waves_htxk=(1+0.5*horndeski_alpm_eff__mod__gravitational_waves_htxk)/(scale_factor__mod__gravitational_waves_htxk*scale_factor__mod__gravitational_waves_htxk)
          horndeski_alpm_eff2__mod__gravitational_waves_htxk=horndeski_alpm_eff2__mod__gravitational_waves_htxk*0.5*horndeski_alpm_eff__mod__gravitational_waves_htxk
          horndeski_alpm_eff3__mod__gravitational_waves_htxk=0.5*AC_horndeski_alpm_prime__mod__gravitational_waves_htxk/scale_factor__mod__gravitational_waves_htxk
          horndeski_alpm_eff__mod__gravitational_waves_htxk=1.+0.5*horndeski_alpm_eff__mod__gravitational_waves_htxk
        }
      }
    }
    real appa_om__mod__gravitational_waves_htxk
    if (AC_lread_scl_factor_file__mod__cdata && AC_lread_scl_factor_file_exists__mod__gravitational_waves_htxk) {
      appa_om__mod__gravitational_waves_htxk=appa_target__mod__cdata
    }
    if (AC_lhorndeski_xi__mod__gravitational_waves_htxk) {
      appa_om__mod__gravitational_waves_htxk=appa_om__mod__gravitational_waves_htxk*horndeski_alpm_eff__mod__gravitational_waves_htxk+horndeski_alpm_eff2__mod__gravitational_waves_htxk
      appa_om__mod__gravitational_waves_htxk=appa_om__mod__gravitational_waves_htxk+horndeski_alpm_eff3__mod__gravitational_waves_htxk
    }
    if(
	!(AC_lread_scl_factor_file__mod__cdata && AC_lread_scl_factor_file_exists__mod__gravitational_waves_htxk) && !AC_lhorndeski_xi__mod__gravitational_waves_htxk
      )
    {
	    appa_om__mod__gravitational_waves_htxk = AC_appa_om_init__mod__gravitational_waves_htxk
    }

    kmin=2*pi/sqrt((AC_lx__mod__cdata*AC_lx__mod__cdata)+(AC_ly__mod__cdata*AC_ly__mod__cdata)+(AC_lz__mod__cdata*AC_lz__mod__cdata))
    om2_min=((1e-4*kmin)*(1e-4*kmin))
    s_t_re=0.
    s_t_im=0.
    s_x_re=0.
    s_x_im=0.
    k1=AC_kx_fft__mod__fourier[ikx+(AC_ipx__mod__cdata*nx)-1]
    k2=AC_ky_fft__mod__fourier[iky+(AC_ipy__mod__cdata*ny)-1]
    k3=AC_kz_fft__mod__fourier[ikz+(AC_ipz__mod__cdata*nz)-1]
    k1sqr=(k1*k1)
    k2sqr=(k2*k2)
    k3sqr=(k3*k3)
    ksqr=k1sqr+k2sqr+k3sqr
    ksqrt = sqrt(ksqr)
    if (AC_lroot__mod__cdata && ikx==1 && iky==1 && ikz==1) {
      e1.x = 0.
      e1.y = 0.
      e1.z = 0.
      e2.x = 0.
      e2.y = 0.
      e2.z = 0.
      pij=0.
      kij=0.
      om=0.
      om2=0.
    }
    else {
      one_over_k2=1./ksqr
      if (AC_linflation__mod__gravitational_waves_htxk) {
        om2=4.*ksqr-2./(AC_t__mod__cdata*AC_t__mod__cdata)
        lsign_om2=(om2 >= 0.)
        om=sqrt(abs(om2))
      }
      else if (AC_lreheating_gw__mod__gravitational_waves_htxk) {
        om2=ksqr-2./((AC_t__mod__cdata+1.)*(AC_t__mod__cdata+1.))
        lsign_om2=(om2 >= 0.)
        om=sqrt(abs(om2))
      }
      else if (AC_lscalar__mod__gravitational_waves_htxk) {
        om2=ksqr-0.0
        lsign_om2=(om2 >= 0.)
        om=sqrt(abs(om2))
      }
      else if (AC_lmatter_gw__mod__gravitational_waves_htxk  ||  AC_ldark_energy_gw__mod__gravitational_waves_htxk) {
        om2=ksqr-2./(AC_t__mod__cdata*AC_t__mod__cdata)
        lsign_om2=(om2 >= 0.)
        om=sqrt(abs(om2))
      }
      else {
        if (delkt!=0.  ||  AC_lhorndeski__mod__gravitational_waves_htxk) {
          if (AC_lhorndeski__mod__gravitational_waves_htxk) {
            om2=(1.+horndeski_alpt_eff__mod__gravitational_waves_htxk)*ksqr+(delkt*delkt)-horndeski_alpm_eff2__mod__gravitational_waves_htxk-appa_om__mod__gravitational_waves_htxk
            om_cmplx=sqrt(cmplx(om2,0.))
            om=AC_impossible__mod__cparam
          }
          else if (AC_lhorndeski_xi__mod__gravitational_waves_htxk) {
            om2=(1.+horndeski_alpt_eff__mod__gravitational_waves_htxk)*ksqr+(delkt*delkt)-appa_om__mod__gravitational_waves_htxk
          }
          else {
            om2=ksqr+(delkt*delkt)-appa_om__mod__gravitational_waves_htxk
            om=sqrt(om2)
          }
        }
        else {
          om2=ksqr-appa_om__mod__gravitational_waves_htxk
          om=sqrt(om2)
        }
        lsign_om2=true
      }
      if(abs(k1)<abs(k2)) {
        if(abs(k1)<abs(k3)) {
          e1=real3(0., -k3, +k2)
          e2=real3(k2sqr+k3sqr, -k2*k1, -k3*k1)
        }
        else {
          e1=real3(k2, -k1, 0.)
          e2=real3(k1*k3, k2*k3, -(k1sqr+k2sqr))
        }
      }
      else {
        if(abs(k2)<abs(k3)) {
          e1=real3(-k3, 0., +k1)
          e2=real3(+k1*k2, -(k1sqr+k3sqr), +k3*k2)
        }
        else {
          e1=real3(k2, -k1, 0.)
          e2=real3(k1*k3, k2*k3, -(k1sqr+k2sqr))
        }
      }
      e1=e1/sqrt((e1.x*e1.x)+(e1.y*e1.y)+(e1.z*e1.z))
      e2=e2/sqrt((e2.x*e2.x)+(e2.y*e2.y)+(e2.z*e2.z))
      pij[1-1]=1.-k1sqr*one_over_k2
      pij[2-1]=1.-k2sqr*one_over_k2
      pij[3-1]=1.-k3sqr*one_over_k2
      pij[4-1]=-k1*k2*one_over_k2
      pij[5-1]=-k2*k3*one_over_k2
      pij[6-1]=-k3*k1*one_over_k2
      if (AC_llighthill__mod__gravitational_waves_htxk) {
        kij[1-1]=-k1sqr
        kij[2-1]=-k2sqr
        kij[3-1]=-k3sqr
        kij[4-1]=-k1*k2
        kij[5-1]=-k2*k3
        kij[6-1]=-k3*k1
      }
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    e_t[ij-1]=e1.x*e1.x-e2.x*e2.x
    e_x[ij-1]=e1.x*e2.x+e2.x*e1.x
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    e_t[ij-1]=e1.y*e1.x-e2.y*e2.x
    e_x[ij-1]=e1.y*e2.x+e2.y*e1.x
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    e_t[ij-1]=e1.z*e1.x-e2.z*e2.x
    e_x[ij-1]=e1.z*e2.x+e2.z*e1.x
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    e_t[ij-1]=e1.x*e1.y-e2.x*e2.y
    e_x[ij-1]=e1.x*e2.y+e2.x*e1.y
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    e_t[ij-1]=e1.y*e1.y-e2.y*e2.y
    e_x[ij-1]=e1.y*e2.y+e2.y*e1.y
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    e_t[ij-1]=e1.z*e1.y-e2.z*e2.y
    e_x[ij-1]=e1.z*e2.y+e2.z*e1.y
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    e_t[ij-1]=e1.x*e1.z-e2.x*e2.z
    e_x[ij-1]=e1.x*e2.z+e2.x*e1.z
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    e_t[ij-1]=e1.y*e1.z-e2.y*e2.z
    e_x[ij-1]=e1.y*e2.z+e2.y*e1.z
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    e_t[ij-1]=e1.z*e1.z-e2.z*e2.z
    e_x[ij-1]=e1.z*e2.z+e2.z*e1.z
    if (AC_lswitch_sign_e_x__mod__gravitational_waves_htxk) {
      if (k3<0.) {
        e_x=-e_x
      }
      else if (k3==0.) {
        if (k2<0.) {
          e_x=-e_x
        }
        else if (k2==0.) {
          if (k1<0.) {
            e_x=-e_x
          }
        }
      }
    }
    sij_re=0.
    sij_im=0.
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    pq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    ip=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    jq=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_re__mod__gravitational_waves_htxk[pq-1]
    sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_tpq_im__mod__gravitational_waves_htxk[pq-1]
    if (AC_lnonlinear_source__mod__gravitational_waves_htxk && AC_lnonlinear_tpq_trans__mod__gravitational_waves_htxk) {
      sij_re[ij-1]=sij_re[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_re__mod__gravitational_waves_htxk[pq-1]
      sij_im[ij-1]=sij_im[ij-1]+(pij[ip-1]*pij[jq-1]-0.5*pij[ij-1]*pij[pq-1])*AC_nonlinear_tpq_im__mod__gravitational_waves_htxk[pq-1]
    }
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
    s_t_re=s_t_re+0.5*e_t[ij-1]*sij_re[ij-1]
    s_t_im=s_t_im+0.5*e_t[ij-1]*sij_im[ij-1]
    s_x_re=s_x_re+0.5*e_x[ij-1]*sij_re[ij-1]
    s_x_im=s_x_im+0.5*e_x[ij-1]*sij_im[ij-1]
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
    s_t_re=s_t_re+0.5*e_t[ij-1]*sij_re[ij-1]
    s_t_im=s_t_im+0.5*e_t[ij-1]*sij_im[ij-1]
    s_x_re=s_x_re+0.5*e_x[ij-1]*sij_re[ij-1]
    s_x_im=s_x_im+0.5*e_x[ij-1]*sij_im[ij-1]
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
    s_t_re=s_t_re+0.5*e_t[ij-1]*sij_re[ij-1]
    s_t_im=s_t_im+0.5*e_t[ij-1]*sij_im[ij-1]
    s_x_re=s_x_re+0.5*e_x[ij-1]*sij_re[ij-1]
    s_x_im=s_x_im+0.5*e_x[ij-1]*sij_im[ij-1]
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
    s_t_re=s_t_re+0.5*e_t[ij-1]*sij_re[ij-1]
    s_t_im=s_t_im+0.5*e_t[ij-1]*sij_im[ij-1]
    s_x_re=s_x_re+0.5*e_x[ij-1]*sij_re[ij-1]
    s_x_im=s_x_im+0.5*e_x[ij-1]*sij_im[ij-1]
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
    s_t_re=s_t_re+0.5*e_t[ij-1]*sij_re[ij-1]
    s_t_im=s_t_im+0.5*e_t[ij-1]*sij_im[ij-1]
    s_x_re=s_x_re+0.5*e_x[ij-1]*sij_re[ij-1]
    s_x_im=s_x_im+0.5*e_x[ij-1]*sij_im[ij-1]
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
    s_t_re=s_t_re+0.5*e_t[ij-1]*sij_re[ij-1]
    s_t_im=s_t_im+0.5*e_t[ij-1]*sij_im[ij-1]
    s_x_re=s_x_re+0.5*e_x[ij-1]*sij_re[ij-1]
    s_x_im=s_x_im+0.5*e_x[ij-1]*sij_im[ij-1]
    ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
    s_t_re=s_t_re+0.5*e_t[ij-1]*sij_re[ij-1]
    s_t_im=s_t_im+0.5*e_t[ij-1]*sij_im[ij-1]
    s_x_re=s_x_re+0.5*e_x[ij-1]*sij_re[ij-1]
    s_x_im=s_x_im+0.5*e_x[ij-1]*sij_im[ij-1]
    ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
    s_t_re=s_t_re+0.5*e_t[ij-1]*sij_re[ij-1]
    s_t_im=s_t_im+0.5*e_t[ij-1]*sij_im[ij-1]
    s_x_re=s_x_re+0.5*e_x[ij-1]*sij_re[ij-1]
    s_x_im=s_x_im+0.5*e_x[ij-1]*sij_im[ij-1]
    ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
    s_t_re=s_t_re+0.5*e_t[ij-1]*sij_re[ij-1]
    s_t_im=s_t_im+0.5*e_t[ij-1]*sij_im[ij-1]
    s_x_re=s_x_re+0.5*e_x[ij-1]*sij_re[ij-1]
    s_x_im=s_x_im+0.5*e_x[ij-1]*sij_im[ij-1]
    if (AC_llighthill__mod__gravitational_waves_htxk) {
      ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][1-1]
      s_t_re=s_t_re+kij[ij-1]*AC_tpq_re__mod__gravitational_waves_htxk[ij-1]
      s_t_im=s_t_im+kij[ij-1]*AC_tpq_im__mod__gravitational_waves_htxk[ij-1]
      ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][1-1]
      s_t_re=s_t_re+kij[ij-1]*AC_tpq_re__mod__gravitational_waves_htxk[ij-1]
      s_t_im=s_t_im+kij[ij-1]*AC_tpq_im__mod__gravitational_waves_htxk[ij-1]
      ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][1-1]
      s_t_re=s_t_re+kij[ij-1]*AC_tpq_re__mod__gravitational_waves_htxk[ij-1]
      s_t_im=s_t_im+kij[ij-1]*AC_tpq_im__mod__gravitational_waves_htxk[ij-1]
      ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][2-1]
      s_t_re=s_t_re+kij[ij-1]*AC_tpq_re__mod__gravitational_waves_htxk[ij-1]
      s_t_im=s_t_im+kij[ij-1]*AC_tpq_im__mod__gravitational_waves_htxk[ij-1]
      ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][2-1]
      s_t_re=s_t_re+kij[ij-1]*AC_tpq_re__mod__gravitational_waves_htxk[ij-1]
      s_t_im=s_t_im+kij[ij-1]*AC_tpq_im__mod__gravitational_waves_htxk[ij-1]
      ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][2-1]
      s_t_re=s_t_re+kij[ij-1]*AC_tpq_re__mod__gravitational_waves_htxk[ij-1]
      s_t_im=s_t_im+kij[ij-1]*AC_tpq_im__mod__gravitational_waves_htxk[ij-1]
      ij=AC_ij_table__mod__gravitational_waves_htxk[1-1][3-1]
      s_t_re=s_t_re+kij[ij-1]*AC_tpq_re__mod__gravitational_waves_htxk[ij-1]
      s_t_im=s_t_im+kij[ij-1]*AC_tpq_im__mod__gravitational_waves_htxk[ij-1]
      ij=AC_ij_table__mod__gravitational_waves_htxk[2-1][3-1]
      s_t_re=s_t_re+kij[ij-1]*AC_tpq_re__mod__gravitational_waves_htxk[ij-1]
      s_t_im=s_t_im+kij[ij-1]*AC_tpq_im__mod__gravitational_waves_htxk[ij-1]
      ij=AC_ij_table__mod__gravitational_waves_htxk[3-1][3-1]
      s_t_re=s_t_re+kij[ij-1]*AC_tpq_re__mod__gravitational_waves_htxk[ij-1]
      s_t_im=s_t_im+kij[ij-1]*AC_tpq_im__mod__gravitational_waves_htxk[ij-1]
    }
    if (AC_lnophase_in_stress__mod__gravitational_waves_htxk) {
      if (AC_lconstmod_in_stress__mod__gravitational_waves_htxk) {
        s_t_re=exp(-ksqr/(AC_k_in_stress__mod__gravitational_waves_htxk*AC_k_in_stress__mod__gravitational_waves_htxk))
        s_x_re=exp(-ksqr/(AC_k_in_stress__mod__gravitational_waves_htxk*AC_k_in_stress__mod__gravitational_waves_htxk))
      }
      else {
        if (ksqr==0.) {
          s_t_re=0.
          s_x_re=0.
        }
        else {
          s_t_re=sqrt((s_t_re*s_t_re)+(s_t_im*s_t_im))
          s_x_re=sqrt((s_x_re*s_x_re)+(s_x_im*s_x_im))
        }
      }
      s_t_im=0.
      s_x_im=0.
      if (AC_llinphase_in_stress__mod__gravitational_waves_htxk) {
        s_t_re=s_t_re*cos(AC_slope_linphase_in_stress__mod__gravitational_waves_htxk*AC_t__mod__cdata)
        s_t_im=s_t_re*sin(AC_slope_linphase_in_stress__mod__gravitational_waves_htxk*AC_t__mod__cdata)
        s_x_re=s_x_re*cos(AC_slope_linphase_in_stress__mod__gravitational_waves_htxk*AC_t__mod__cdata)
        s_x_im=s_x_re*sin(AC_slope_linphase_in_stress__mod__gravitational_waves_htxk*AC_t__mod__cdata)
      }
    }

    hhtre=value(F_HHT)
    hhxre=value(F_HHX)
    hhtim=value(F_HHTIM)
    hhxim=value(F_HHXIM)
    ggtre=value(F_GGT)
    ggxre=value(F_GGX)
    ggtim=value(F_GGTIM)
    ggxim=value(F_GGXIM)

    om12=1./om2
    if (AC_lhorndeski__mod__gravitational_waves_htxk) {
      discrim2=(horndeski_alpm_eff__mod__gravitational_waves_htxk*horndeski_alpm_eff__mod__gravitational_waves_htxk)-4.*om2
      if (discrim2==0.) {
        discrim2=AC_tini__mod__cparam
      }
      discrim=sqrt(cmplx(discrim2,0.))
      lam1=0.5*(-horndeski_alpm_eff__mod__gravitational_waves_htxk+discrim)
      lam2=0.5*(-horndeski_alpm_eff__mod__gravitational_waves_htxk-discrim)
      explam1t=exp(lam1*AC_dt__mod__cdata)
      explam2t=exp(lam2*AC_dt__mod__cdata)
      det1=1./discrim
      cosoth=det1*(lam1*explam2t-lam2*explam1t)
      cosotg=det1*(lam1*explam1t-lam2*explam2t)
      sinoth=-det1*(     explam2t-     explam1t)*om_cmplx
      sinotg=+det1*(     explam2t-     explam1t)/om_cmplx*lam1*lam2
    }
    else {
      if (lsign_om2) {
        cosot=cos(om*AC_dt__mod__cdata)
        sinot=sin(om*AC_dt__mod__cdata)
        sinot_minus=-sinot
      }
      else {
        cosot=cosh(om*AC_dt__mod__cdata)
        sinot=sinh(om*AC_dt__mod__cdata)
        sinot_minus=+sinot
      }
    }
    if (AC_lhorndeski__mod__gravitational_waves_htxk) {
      coefa=cmplx(hhtre-om12*s_t_re,hhtim-om12*s_t_im)
      coefb=cmplx(ggtre                         ,ggtim    )/om_cmplx
      hcomplex_new= cosoth*coefa+sinoth*coefb+om12*cmplx(s_t_re,s_t_im)
      gcomplex_new=(sinotg*coefa+cosotg*coefb)*om_cmplx
      DF_HHT= hcomplex_new.x
      DF_HHTIM=aimag(hcomplex_new)
      DF_GGT= gcomplex_new.x
      DF_GGTIM=aimag(gcomplex_new)
    }
    else {
      om1=1./om
      coefare=(hhtre-om12*s_t_re)
      coefaim=(hhtim-om12*s_t_im)
      coefbre=ggtre*om1
      coefbim=ggtim*om1
      DF_HHT=coefare*cosot+coefbre*sinot+om12*s_t_re
      DF_HHTIM=coefaim*cosot+coefbim*sinot+om12*s_t_im
      DF_GGT=coefbre*cosot*om+coefare*om*sinot_minus
      DF_GGTIM=coefbim*cosot*om+coefaim*om*sinot_minus
      if (AC_itorder_gw__mod__gravitational_waves_htxk==2) {
        if (AC_dt__mod__cdata==0.) {
          dt1=0.
        }
        else {
          dt1=1./AC_dt__mod__cdata
        }
        ds_t_re=s_t_re-value(F_STRESST)
        ds_t_im=s_t_im-value(F_STRESSTIM)
        DF_HHT   +=   ds_t_re*om12*(1.-om1*dt1*sinot)
        DF_HHTIM +=   ds_t_im*om12*(1.-om1*dt1*sinot)
        DF_GGT   +=   ds_t_re*om12*dt1*(1.-cosot)
        DF_GGTIM +=   ds_t_im*om12*dt1*(1.-cosot)
      }
    }
    if (AC_lhorndeski__mod__gravitational_waves_htxk) {
      coefa=cmplx(hhxre-om12*s_x_re,hhxim-om12*s_x_im)
      coefb=cmplx(ggxre                         ,ggxim    )/om_cmplx
      hcomplex_new= cosoth*coefa+sinoth*coefb+om12*cmplx(s_x_re,s_x_im)
      gcomplex_new=(sinotg*coefa+cosotg*coefb)*om_cmplx
      DF_HHX= hcomplex_new.x
      DF_HHXIM=aimag(hcomplex_new)
      DF_GGX= gcomplex_new.x
      DF_GGXIM=aimag(gcomplex_new)
    }
    else {
      coefare=(hhxre-om12*s_x_re)
      coefaim=(hhxim-om12*s_x_im)
      coefbre=ggxre*om1
      coefbim=ggxim*om1
      DF_HHX=coefare*cosot+coefbre*sinot+om12*s_x_re
      DF_HHXIM=coefaim*cosot+coefbim*sinot+om12*s_x_im
      DF_GGX=coefbre*cosot*om+coefare*om*sinot_minus
      DF_GGXIM=coefbim*cosot*om+coefaim*om*sinot_minus
      if (AC_itorder_gw__mod__gravitational_waves_htxk==2) {
        ds_x_re=s_x_re-value(F_STRESSX)
        ds_x_im=s_x_im-value(F_STRESSXIM)
        DF_HHX    += ds_x_re*om12*(AC_dt__mod__cdata-om1*sinot)
        DF_HHXIM  += ds_x_im*om12*(AC_dt__mod__cdata-om1*sinot)
        DF_GGX    += ds_x_re*om12*(1.-cosot)
        DF_GGXIM  += ds_x_im*om12*(1.-cosot)
      }
    }
    if(om2 <= om2_min)
    {
      DF_HHT    = 0. 
      DF_HHTIM  = 0. 
      DF_GGT    = 0. 
      DF_GGTIM  = 0. 
      DF_HHX    = 0. 
      DF_HHXIM  = 0. 
      DF_GGX    = 0. 
      DF_GGXIM  = 0. 
    }
    if (AC_itorder_gw__mod__gravitational_waves_htxk==2)
    {
        write(F_STRESST  ,s_t_re)
        write(F_STRESSTIM,s_t_im)
        write(F_STRESSX  ,s_x_re)
        write(F_STRESSXIM,s_x_im)
    }

    write(F_HHT,DF_HHT)
    write(F_HHX,DF_HHX)
    write(F_HHTIM,DF_HHTIM)
    write(F_HHXIM,DF_HHXIM)

    write(F_GGT,DF_GGT)
    write(F_GGX,DF_GGX)
    write(F_GGTIM,DF_GGTIM)
    write(F_GGXIM,DF_GGXIM)
}

#else
Kernel gravitational_waves_solve_and_stress(real AC_t__mod__cdata, real AC_dt__mod__cdata){
	suppress_unused_warning(AC_t__mod__cdata)
	suppress_unused_warning(AC_dt__mod__cdata)
}
#endif

