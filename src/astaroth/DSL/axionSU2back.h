#if LAXIONSU2BACK
global output real AC_grand_sum
global output real AC_dgrant_sum

global output real AC_trdoteff2km_sum
global output real AC_trdoteff2m_sum
global output real AC_treff2km_sum
global output real AC_treff2m_sum
global output real AC_tldoteff2km_sum
global output real AC_tldoteff2m_sum
global output real AC_tleff2km_sum
global output real AC_tleff2m_sum

field_order(AC_iaxi_psi__mod__axionsu2back-1) Field F_AXI_PSI    
field_order(AC_iaxi_psidot__mod__axionsu2back-1) Field F_AXI_PSIDOT 
field_order(AC_iaxi_tr__mod__axionsu2back-1) Field F_AXI_TR     
field_order(AC_iaxi_trdot__mod__axionsu2back-1) Field F_AXI_TRDOT  

field_order(AC_iaxi_impsi__mod__axionsu2back-1) Field F_AXI_IMPSI  
field_order(AC_iaxi_impsidot__mod__axionsu2back-1) Field F_AXI_IMPSIDOT 
field_order(AC_iaxi_imtr__mod__axionsu2back-1) Field F_AXI_IMTR    
field_order(AC_iaxi_imtrdot__mod__axionsu2back-1) Field F_AXI_IMTRDOT 

field_order(AC_iaxi_psil__mod__axionsu2back-1) Field F_AXI_PSIL   
field_order(AC_iaxi_psildot__mod__axionsu2back-1) Field F_AXI_PSILDOT
field_order(AC_iaxi_tl__mod__axionsu2back-1) Field F_AXI_TL    
field_order(AC_iaxi_tldot__mod__axionsu2back-1) Field F_AXI_TLDOT 
field_order(AC_iaxi_impsil__mod__axionsu2back-1) Field F_AXI_IMPSIL   
field_order(AC_iaxi_impsildot__mod__axionsu2back-1) Field F_AXI_IMPSILDOT
field_order(AC_iaxi_imtl__mod__axionsu2back-1) Field F_AXI_IMTL    
field_order(AC_iaxi_imtldot__mod__axionsu2back-1) Field F_AXI_IMTLDOT 

field_order(AC_iaxi_ur__mod__axionsu2back-1)      Field F_AXI_UR
field_order(AC_iaxi_urdot__mod__axionsu2back-1)   Field F_AXI_URDOT
field_order(AC_iaxi_imur__mod__axionsu2back-1)    Field F_AXI_IMUR
field_order(AC_iaxi_imurdot__mod__axionsu2back-1) Field F_AXI_IMURDOT

field_order(AC_iaxi_ul__mod__axionsu2back-1)      Field F_AXI_UL
field_order(AC_iaxi_uldot__mod__axionsu2back-1)   Field F_AXI_ULDOT
field_order(AC_iaxi_imul__mod__axionsu2back-1)    Field F_AXI_IMUL
field_order(AC_iaxi_imuldot__mod__axionsu2back-1) Field F_AXI_IMULDOT

Kernel calc_axion_integral(real AC_t__mod__cdata){
  real trdoteff2m
  real trdoteff2km
  real treff2m
  real treff2km

  real tldoteff2m
  real tldoteff2km
  real tleff2m
  real tleff2km

  real tr
  real trdot
  real imtr
  real imtrdot
  real treff2
  real trdoteff2
  real tl
  real tldot
  real imtl
  real imtldot
  real tleff2
  real tldoteff2
  real psi
  real psidot
  real impsi
  real impsidot
  real trpsi
  real trpsik
  real trpsidot
  real trdotpsi
  real xi
  real chidot
  real mq
  real qdot
  real a
  real inflaton
  real a_return_value_0
  real epsilon_sr_0
  real inflaton_0
  real get_mq_return_value_1

  real h__mod__axionsu2back
  real grand__mod__axionsu2back
  real dgrant__mod__axionsu2back

  if (! AC_lconf_time__mod__axionsu2back  &&  AC_lhubble_var__mod__axionsu2back) {
    inflaton=AC_inflaton_ini__mod__axionsu2back-sqrt(2./3.)*AC_m_inflaton__mod__axionsu2back*AC_t__mod__cdata
    h__mod__axionsu2back=0.41*AC_m_inflaton__mod__axionsu2back*inflaton
  }
  else if (lgpu &&  ! AC_lhubble__mod__axionsu2back) {
    h__mod__axionsu2back=AC_h_init__mod__axionsu2back
  }
  if (AC_lconf_time__mod__axionsu2back) {
    if (AC_lhubble_var__mod__axionsu2back) {
      epsilon_sr_0=0.8*(1+tanh(0.3*(log(-1/(h__mod__axionsu2back*AC_t__mod__cdata))-18)))*0.5
      a_return_value_0=-1./(h__mod__axionsu2back*AC_t__mod__cdata*(1-epsilon_sr_0))
    }
    else {
      a_return_value_0=-1./(h__mod__axionsu2back*AC_t__mod__cdata)
    }
  }
  else if (AC_lhubble_var__mod__axionsu2back) {
    inflaton_0=AC_inflaton_ini__mod__axionsu2back-sqrt(2./3.)*AC_m_inflaton__mod__axionsu2back*AC_t__mod__cdata
    a_return_value_0=exp(((AC_inflaton_ini__mod__axionsu2back*AC_inflaton_ini__mod__axionsu2back)-(inflaton_0*inflaton_0))*0.25)
  }
  else {
    if (AC_lhubble__mod__axionsu2back) {
      a_return_value_0 = exp(AC_f_ode__mod__cdata[AC_iaxi_lna__mod__axionsu2back-1])
    }
    else {
      a_return_value_0 = exp(h__mod__axionsu2back*AC_t__mod__cdata)
    }
  }
  a = a_return_value_0
  chidot=AC_f_ode__mod__cdata[AC_iaxi_chidot__mod__axionsu2back-1]
  qdot  =AC_f_ode__mod__cdata[AC_iaxi_qdot__mod__axionsu2back-1]
  if (AC_lconf_time__mod__axionsu2back) {
    if (AC_lhubble_var__mod__axionsu2back) {
      xi=AC_lamf__mod__axionsu2back*chidot*(0.5/(a*h__mod__axionsu2back))
    }
    else {
      xi=AC_lamf__mod__axionsu2back*chidot*(-0.5*AC_t__mod__cdata)
    }
  }
  else {
    xi=AC_lamf__mod__axionsu2back*chidot/(2.*h__mod__axionsu2back)
  }
  if (AC_lkeep_mq_const__mod__axionsu2back) {
    get_mq_return_value_1=AC_g__mod__axionsu2back*AC_q0__mod__axionsu2back/h__mod__axionsu2back
  }
  else {
    get_mq_return_value_1=AC_g__mod__axionsu2back*AC_f_ode__mod__cdata[AC_iaxi_q__mod__axionsu2back-1]/h__mod__axionsu2back
  }
  mq = get_mq_return_value_1
  psi   =value(Field(AC_iaxi_psi__mod__axionsu2back-1))
  psidot=value(Field(AC_iaxi_psidot__mod__axionsu2back-1))
  tr   =value(Field(AC_iaxi_tr__mod__axionsu2back-1))
  trdot=value(Field(AC_iaxi_trdot__mod__axionsu2back-1))
  if (AC_lim_psi_tr__mod__axionsu2back) {
    impsi   =value(Field(AC_iaxi_impsi__mod__axionsu2back-1))
    impsidot=value(Field(AC_iaxi_impsidot__mod__axionsu2back-1))
    imtr   =value(Field(AC_iaxi_imtr__mod__axionsu2back-1))
    imtrdot=value(Field(AC_iaxi_imtrdot__mod__axionsu2back-1))
  }
  treff2=(tr*tr)
  trdoteff2=tr*trdot
  trpsi=tr*psi
  trpsidot=tr*psidot
  trdotpsi=trdot*psi
  if (AC_lim_psi_tr__mod__axionsu2back) {
    treff2=treff2+(imtr*imtr)
    trdoteff2=trdoteff2+imtr*imtrdot
    trpsi=trpsi+imtr*impsi
    trpsidot=trpsidot+imtr*impsidot
    trdotpsi=trdotpsi+imtrdot*impsi
  }
  if (AC_horizon_factor__mod__axionsu2back==0.) {
    if (treff2<1./(2.*a*h__mod__axionsu2back)) {
      treff2=0.
      trdoteff2=0.
      trpsi=0.
      trpsidot=0.
      trdotpsi=0.
    }
  }
  else if (AC_horizon_factor__mod__axionsu2back>0.) {
    if (AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]>(a*h__mod__axionsu2back*AC_horizon_factor__mod__axionsu2back)) {
      treff2=0.
      trdoteff2=0.
      trpsi=0.
      trpsidot=0.
      trdotpsi=0.
    }
  }
  if (AC_lleft_psil_tl__mod__axionsu2back) {
    tl   =value(Field(AC_iaxi_tl__mod__axionsu2back-1))
    tldot=value(Field(AC_iaxi_tldot__mod__axionsu2back-1))
    imtl   =value(Field(AC_iaxi_imtl__mod__axionsu2back-1))
    imtldot=value(Field(AC_iaxi_imtldot__mod__axionsu2back-1))
    tleff2=(tl*tl)
    tldoteff2=tl*tldot
    tleff2=tleff2+(imtl*imtl)
    tldoteff2=tldoteff2+imtl*imtldot
    if (AC_horizon_factor__mod__axionsu2back==0.) {
      if (tleff2<1./(2.*a*h__mod__axionsu2back)) {
        tleff2=0.
        tldoteff2=0.
      }
    }
    else if (AC_horizon_factor__mod__axionsu2back>0.) {
      if (AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]>(a*h__mod__axionsu2back*AC_horizon_factor__mod__axionsu2back)) {
        tleff2=0.
        tldoteff2=0.
      }
    }
  }
  if (AC_llnk_spacing_adjustable__mod__axionsu2back  ||  AC_llnk_spacing__mod__axionsu2back) {
    trpsim=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*trpsi
    trpsikm=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*trpsi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)
    trpsidotm=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*trpsidot
    trdotpsim=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*trdotpsi
    trdoteff2km=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*trdoteff2*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)
    trdoteff2m=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*trdoteff2
    treff2km=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*treff2*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)
    treff2m=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*treff2
    grand__mod__axionsu2back=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*(xi*h__mod__axionsu2back-AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)*treff2*(+   AC_g__mod__axionsu2back/(3.*(a*a)))/(twopi*twopi*twopi)
    grant__mod__axionsu2back=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*(mq*h__mod__axionsu2back-AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)*treff2*(-AC_lamf__mod__axionsu2back/(2.*(a*a)))/(twopi*twopi*twopi)
    if (AC_lleft_psil_tl__mod__axionsu2back) {
      tldoteff2km=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*tldoteff2*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)
      tldoteff2m=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*tldoteff2
      tleff2km=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*tleff2*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)
      tleff2m=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*tleff2
      grand__mod__axionsu2back=grand__mod__axionsu2back+(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*(xi*h__mod__axionsu2back+AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)*tleff2*(+   AC_g__mod__axionsu2back/(3.*(a*a)))/(twopi*twopi*twopi)
      grant__mod__axionsu2back=grant__mod__axionsu2back+(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*(mq*h__mod__axionsu2back+AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)*tleff2*(-AC_lamf__mod__axionsu2back/(2.*(a*a)))/(twopi*twopi*twopi)
    }
    if (AC_lconf_time__mod__axionsu2back) {
      dgrant__mod__axionsu2back=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*(-AC_lamf__mod__axionsu2back/(2.*(a*a*a)))*(  (a*mq*(h__mod__axionsu2back*h__mod__axionsu2back)+AC_g__mod__axionsu2back*qdot)*treff2+(mq*h__mod__axionsu2back-AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)*2*trdoteff2  )/(twopi*twopi*twopi)
      if (AC_lleft_psil_tl__mod__axionsu2back) {
        dgrant__mod__axionsu2back=dgrant__mod__axionsu2back+(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*(-AC_lamf__mod__axionsu2back/(2.*(a*a*a)))*(  (a*mq*(h__mod__axionsu2back*h__mod__axionsu2back)+AC_g__mod__axionsu2back*qdot)*tleff2+(mq*h__mod__axionsu2back+AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)*2*tldoteff2  )/(twopi*twopi*twopi)
      }
    }
    else {
      dgrant__mod__axionsu2back=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*(-AC_lamf__mod__axionsu2back/(2.*(a*a*a)))*(  (a*mq*(h__mod__axionsu2back*h__mod__axionsu2back)+a*AC_g__mod__axionsu2back*qdot)*treff2+(a*mq*h__mod__axionsu2back-AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*2*trdoteff2  )/(twopi*twopi*twopi)
      if (AC_lleft_psil_tl__mod__axionsu2back) {
        dgrant__mod__axionsu2back=dgrant__mod__axionsu2back+(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionsu2back)*(-AC_lamf__mod__axionsu2back/(2.*(a*a*a)))*(  (a*mq*(h__mod__axionsu2back*h__mod__axionsu2back)+a*AC_g__mod__axionsu2back*qdot)*tleff2+(a*mq*h__mod__axionsu2back-AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*2*tldoteff2  )/(twopi*twopi*twopi)
      }
    }
  }
  else {
    trpsidotm=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionsu2back)*trpsidot
    trdotpsim=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionsu2back)*trdotpsi
    trdoteff2km=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionsu2back)*trdoteff2*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)
    trdoteff2m=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionsu2back)*trdoteff2
    treff2km=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionsu2back)*treff2*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)
    treff2m=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionsu2back)*treff2
    grand__mod__axionsu2back=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionsu2back)*(xi*h__mod__axionsu2back-AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)*treff2*(+   AC_g__mod__axionsu2back/(3.*(a*a)))/(twopi*twopi*twopi)
    grant__mod__axionsu2back=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionsu2back)*(mq*h__mod__axionsu2back-AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)*treff2*(-AC_lamf__mod__axionsu2back/(2.*(a*a)))/(twopi*twopi*twopi)
    if (AC_lconf_time__mod__axionsu2back) {
      dgrant__mod__axionsu2back=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionsu2back)*(-AC_lamf__mod__axionsu2back/(2.*(a*a*a)))*(  (a*mq*(h__mod__axionsu2back*h__mod__axionsu2back)+AC_g__mod__axionsu2back*qdot)*treff2+(mq*h__mod__axionsu2back-AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]/a)*2*trdoteff2  )/(twopi*twopi*twopi)
    }
    else {
      dgrant__mod__axionsu2back=(4.*pi*(AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionsu2back)*(-AC_lamf__mod__axionsu2back/(2.*(a*a*a)))*(  (a*mq*(h__mod__axionsu2back*h__mod__axionsu2back)+a*AC_g__mod__axionsu2back*qdot)*treff2+(a*mq*h__mod__axionsu2back-AC_k__mod__axionsu2back[vertexIdx.x-NGHOST_VAL])*2*trdoteff2  )/(twopi*twopi*twopi)
    }
  }
  reduce_sum(grand__mod__axionsu2back, AC_grand_sum)
  reduce_sum(dgrant__mod__axionsu2back,AC_dgrant_sum)
  reduce_sum(trdoteff2km,AC_trdoteff2km_sum)
  reduce_sum(trdoteff2m,AC_trdoteff2m_sum)
  reduce_sum(treff2km,AC_treff2km_sum)
  reduce_sum(treff2m,AC_treff2m_sum)
  reduce_sum(tldoteff2km,AC_tldoteff2km_sum)
  reduce_sum(tldoteff2m ,AC_tldoteff2m_sum)
  reduce_sum(tleff2km,   AC_tleff2km_sum)
  reduce_sum(tleff2m    ,AC_tleff2m_sum)
}
#else
#if LAXIONU1BACK
global output real AC_edotb_sum__mod__axionu1back
global output real AC_rhoe__mod__axionu1back
global output real AC_rhob__mod__axionu1back


field_order(AC_iaxi_ar__mod__axionu1back-1)    Field F_AXI_AR
field_order(AC_iaxi_ardot__mod__axionu1back-1) Field F_AXI_ARDOT 
field_order(AC_iaxi_al__mod__axionu1back-1)    Field F_AXI_AL
field_order(AC_iaxi_aldot__mod__axionu1back-1) Field F_AXI_ALDOT

field_order(AC_iaxi_imar__mod__axionu1back-1)    Field F_AXI_IMAR
field_order(AC_iaxi_imardot__mod__axionu1back-1) Field F_AXI_IMARDOT 
field_order(AC_iaxi_imal__mod__axionu1back-1)    Field F_AXI_IMAL
field_order(AC_iaxi_imaldot__mod__axionu1back-1) Field F_AXI_IMALDOT

Kernel calc_axion_integral(real t){
  suppress_unused_warning(t)
  real ar
  real ardot
  real imar
  real imardot
  real areff2
  real arardoteff
  real imarardoteff
  real ardoteff2
  real imareff2
  real imardoteff2
  real al
  real aldot
  real imal
  real imaldot
  real alaldoteff
  real imalaldoteff
  real aleff2
  real aldoteff2
  real imaleff2
  real imaldoteff2
  real psi
  real psidot
  real impsi
  real impsidot
  real a
  real phi_0
  real phidot_0
  real v_0
  real beta_0
  real vprime_0
  real phiddot_0
  real tmp_0
  real alpha_gmssm_0
  real h__mod__axionu1back
  if (AC_lhubble__mod__axionu1back) {
    phi_0=AC_f_ode__mod__cdata[AC_iaxi_phi__mod__axionu1back-1]
    phidot_0=AC_f_ode__mod__cdata[AC_iaxi_phidot__mod__axionu1back-1]
    if(AC_enum_v_choice__mod__axionu1back == enum_alpha_attractors_string) {
      beta_0=sqrt(2./(3.*AC_alpha__mod__axionu1back))
      v_0=AC_alpha__mod__axionu1back*AC_m_alpha__mod__axionu1back*pow(((tanh(beta_0*phi_0/2)*tanh(beta_0*phi_0/2))),AC_n_alpha__mod__axionu1back)
      vprime_0=AC_alpha__mod__axionu1back*AC_m_alpha__mod__axionu1back*AC_n_alpha__mod__axionu1back*beta_0*tanh(beta_0*phi_0/2)*(1./(cosh(beta_0*phi_0/2)*cosh(beta_0*phi_0/2)))
    }
    else if(AC_enum_v_choice__mod__axionu1back == enum_quadratic_string)   {
      v_0=0.5*(AC_m_phi__mod__axionu1back*AC_m_phi__mod__axionu1back)*(phi_0*phi_0)
      vprime_0=(AC_m_phi__mod__axionu1back*AC_m_phi__mod__axionu1back)*phi_0
    }
    else if(AC_enum_v_choice__mod__axionu1back == enum_gmssm_string)   {
      tmp_0=phi_0/AC_phi_0__mod__axionu1back
      alpha_gmssm_0=1-AC_alpha1_gmssm__mod__axionu1back
      v_0=(AC_lambda_gmssm__mod__axionu1back*AC_lambda_gmssm__mod__axionu1back*AC_lambda_gmssm__mod__axionu1back*AC_lambda_gmssm__mod__axionu1back)*(0.5*(tmp_0*tmp_0)-(alpha_gmssm_0*onethird)*(tmp_0*tmp_0*tmp_0*tmp_0*tmp_0*tmp_0)+ (alpha_gmssm_0/10.)*(tmp_0*tmp_0*tmp_0*tmp_0*tmp_0*tmp_0*tmp_0*tmp_0*tmp_0*tmp_0))
      vprime_0=(AC_lambda_gmssm__mod__axionu1back*AC_lambda_gmssm__mod__axionu1back*AC_lambda_gmssm__mod__axionu1back*AC_lambda_gmssm__mod__axionu1back)*(tmp_0-2*alpha_gmssm_0*(tmp_0*tmp_0*tmp_0*tmp_0*tmp_0)+ alpha_gmssm_0*(tmp_0*tmp_0*tmp_0*tmp_0*tmp_0*tmp_0*tmp_0*tmp_0*tmp_0))/AC_phi_0__mod__axionu1back
    }
    else {
    }
    if (AC_lbackreact__mod__axionu1back) {
      h__mod__axionu1back=sqrt(8.*pi*onethird*(1/AC_mpl2__mod__axionu1back)*(0.5*(phidot_0*phidot_0)+v_0+AC_rhoe__mod__axionu1back+AC_rhob__mod__axionu1back))
    }
    else {
      h__mod__axionu1back=sqrt(8.*pi*onethird*(1/AC_mpl2__mod__axionu1back)*(0.5*(phidot_0*phidot_0)+v_0))
    }
    phiddot_0=-3.*h__mod__axionu1back*phidot_0-vprime_0
  }
  else if (lgpu)   {
    h__mod__axionu1back=AC_h_init__mod__axionu1back
  }
  if (AC_lhubble__mod__axionu1back) {
    a = exp(AC_f_ode__mod__cdata[AC_iaxi_lna__mod__axionu1back-1])
  }
  else {
    a = 1.
  }
  ar   =value(Field(AC_iaxi_ar__mod__axionu1back-1))
  ardot=value(Field(AC_iaxi_ardot__mod__axionu1back-1))
  al   =value(Field(AC_iaxi_al__mod__axionu1back-1))
  aldot=value(Field(AC_iaxi_aldot__mod__axionu1back-1))
  imar   =value(Field(AC_iaxi_imar__mod__axionu1back-1))
  imardot=value(Field(AC_iaxi_imardot__mod__axionu1back-1))
  imal   =value(Field(AC_iaxi_imal__mod__axionu1back-1))
  imaldot=value(Field(AC_iaxi_imaldot__mod__axionu1back-1))
  areff2=ar*ar
  aleff2=al*al
  ardoteff2=ardot*ardot
  aldoteff2=aldot*aldot
  imareff2=imar*imar
  imaleff2=imal*imal
  imardoteff2=imardot*imardot
  imaldoteff2=imaldot*imaldot
  arardoteff=ar*ardot
  imarardoteff=imar*imardot
  alaldoteff=al*aldot
  imalaldoteff=imal*imaldot
  if (AC_lquant_filter__mod__axionu1back) {
    if (AC_horizon_factor__mod__axionu1back==0.) {
      if (areff2<1./(2.*a*h__mod__axionu1back)) {
        areff2=0.
        aleff2=0.
        imareff2=0.
        imaleff2=0.
        aldoteff2=0.
        ardoteff2=0.
        imardoteff2=0.
        imaldoteff2=0.
        arardoteff=0.
        imarardoteff=0.
        alaldoteff=0.
        imalaldoteff=0.
      }
    }
    else if (AC_horizon_factor__mod__axionu1back>0.) {
      if (AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]>(a*h__mod__axionu1back*AC_horizon_factor__mod__axionu1back)) {
        areff2=0.
        aleff2=0.
        imareff2=0.
        imaleff2=0.
        aldoteff2=0.
        ardoteff2=0.
        imardoteff2=0.
        imaldoteff2=0.
        arardoteff=0.
        imarardoteff=0.
        alaldoteff=0.
        imalaldoteff=0.
      }
    }
  }
  b2__mod__axionu1back    = 0.0
  e2__mod__axionu1back    = 0.0
  edotb__mod__axionu1back = 0.0
  if (AC_llnk_spacing_adjustable__mod__axionu1back  ||  AC_llnk_spacing__mod__axionu1back) {
    b2__mod__axionu1back=0.5*(4.*pi*(AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionu1back)*pow(a,(-4))*(areff2+imareff2+aleff2+imaleff2)/(twopi*twopi*twopi)
    e2__mod__axionu1back=0.5*(4.*pi*(AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionu1back)*pow(a,(-2))*(ardoteff2+imardoteff2+aldoteff2+imaldoteff2)/(twopi*twopi*twopi)
    edotb__mod__axionu1back=(4.*pi*(AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__axionu1back)*a*(arardoteff+imarardoteff-alaldoteff-imalaldoteff)/(twopi*twopi*twopi)
  }
  else {
    b2__mod__axionu1back=0.5*(4.*pi*(AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionu1back)*(areff2+imareff2+aleff2+imaleff2)/(twopi*twopi*twopi)
    e2__mod__axionu1back=0.5*(4.*pi*(AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionu1back)*(a*a)*(ardoteff2+imardoteff2+aldoteff2+imaldoteff2)/(twopi*twopi*twopi)
    edotb__mod__axionu1back=(4.*pi*(AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL]*AC_k__mod__axionu1back[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__axionu1back)*a*(arardoteff+imarardoteff-alaldoteff-imalaldoteff)/(twopi*twopi*twopi)
  }
  reduce_sum(b2__mod__axionu1back,AC_rhob__mod__axionu1back)
  reduce_sum(e2__mod__axionu1back,AC_rhoe__mod__axionu1back)
  reduce_sum(edotb__mod__axionu1back, AC_edotb_sum__mod__axionu1back)
}

#else
Kernel calc_axion_integral(real t){suppress_unused_warning(t)}
#endif
#endif
