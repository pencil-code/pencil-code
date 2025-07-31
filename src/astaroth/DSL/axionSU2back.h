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

  real h__mod__special
  real grand__mod__special
  real dgrant__mod__special

  if (! AC_lconf_time__mod__special  &&  AC_lhubble_var__mod__special) {
    inflaton=AC_inflaton_ini__mod__special-sqrt(2./3.)*AC_m_inflaton__mod__special*AC_t__mod__cdata
    h__mod__special=0.41*AC_m_inflaton__mod__special*inflaton
  }
  else if (AC_lgpu__mod__cparam  &&  ! AC_lhubble__mod__special) {
    h__mod__special=AC_h_init__mod__special
  }
  if (AC_lconf_time__mod__special) {
    if (AC_lhubble_var__mod__special) {
      epsilon_sr_0=0.8*(1+tanh(0.3*(log(-1/(h__mod__special*AC_t__mod__cdata))-18)))*0.5
      a_return_value_0=-1./(h__mod__special*AC_t__mod__cdata*(1-epsilon_sr_0))
    }
    else {
      a_return_value_0=-1./(h__mod__special*AC_t__mod__cdata)
    }
  }
  else if (AC_lhubble_var__mod__special) {
    inflaton_0=AC_inflaton_ini__mod__special-sqrt(2./3.)*AC_m_inflaton__mod__special*AC_t__mod__cdata
    a_return_value_0=exp(((AC_inflaton_ini__mod__special*AC_inflaton_ini__mod__special)-(inflaton_0*inflaton_0))*0.25)
  }
  else {
    if (AC_lhubble__mod__special) {
      a_return_value_0 = exp(AC_f_ode__mod__cdata[AC_iaxi_lna__mod__special-1])
    }
    else {
      a_return_value_0 = exp(h__mod__special*AC_t__mod__cdata)
    }
  }
  a = a_return_value_0
  chidot=AC_f_ode__mod__cdata[AC_iaxi_chidot__mod__special-1]
  qdot  =AC_f_ode__mod__cdata[AC_iaxi_qdot__mod__special-1]
  if (AC_lconf_time__mod__special) {
    if (AC_lhubble_var__mod__special) {
      xi=AC_lamf__mod__special*chidot*(0.5/(a*h__mod__special))
    }
    else {
      xi=AC_lamf__mod__special*chidot*(-0.5*AC_t__mod__cdata)
    }
  }
  else {
    xi=AC_lamf__mod__special*chidot/(2.*h__mod__special)
  }
  if (AC_lkeep_mq_const__mod__special) {
    get_mq_return_value_1=AC_g__mod__special*AC_q0__mod__special/h__mod__special
  }
  else {
    get_mq_return_value_1=AC_g__mod__special*AC_f_ode__mod__cdata[AC_iaxi_q__mod__special-1]/h__mod__special
  }
  mq = get_mq_return_value_1
  psi   =value(Field(AC_iaxi_psi__mod__special-1))
  psidot=value(Field(AC_iaxi_psidot__mod__special-1))
  tr   =value(Field(AC_iaxi_tr__mod__special-1))
  trdot=value(Field(AC_iaxi_trdot__mod__special-1))
  if (AC_lim_psi_tr__mod__special) {
    impsi   =value(Field(AC_iaxi_impsi__mod__special-1))
    impsidot=value(Field(AC_iaxi_impsidot__mod__special-1))
    imtr   =value(Field(AC_iaxi_imtr__mod__special-1))
    imtrdot=value(Field(AC_iaxi_imtrdot__mod__special-1))
  }
  treff2=(tr*tr)
  trdoteff2=tr*trdot
  trpsi=tr*psi
  trpsidot=tr*psidot
  trdotpsi=trdot*psi
  if (AC_lim_psi_tr__mod__special) {
    treff2=treff2+(imtr*imtr)
    trdoteff2=trdoteff2+imtr*imtrdot
    trpsi=trpsi+imtr*impsi
    trpsidot=trpsidot+imtr*impsidot
    trdotpsi=trdotpsi+imtrdot*impsi
  }
  if (AC_horizon_factor__mod__special==0.) {
    if (treff2<1./(2.*a*h__mod__special)) {
      treff2=0.
      trdoteff2=0.
      trpsi=0.
      trpsidot=0.
      trdotpsi=0.
    }
  }
  else if (AC_horizon_factor__mod__special>0.) {
    if (AC_k__mod__special[vertexIdx.x-NGHOST_VAL]>(a*h__mod__special*AC_horizon_factor__mod__special)) {
      treff2=0.
      trdoteff2=0.
      trpsi=0.
      trpsidot=0.
      trdotpsi=0.
    }
  }
  if (AC_lleft_psil_tl__mod__special) {
    tl   =value(Field(AC_iaxi_tl__mod__special-1))
    tldot=value(Field(AC_iaxi_tldot__mod__special-1))
    imtl   =value(Field(AC_iaxi_imtl__mod__special-1))
    imtldot=value(Field(AC_iaxi_imtldot__mod__special-1))
    tleff2=(tl*tl)
    tldoteff2=tl*tldot
    tleff2=tleff2+(imtl*imtl)
    tldoteff2=tldoteff2+imtl*imtldot
    if (AC_horizon_factor__mod__special==0.) {
      if (tleff2<1./(2.*a*h__mod__special)) {
        tleff2=0.
        tldoteff2=0.
      }
    }
    else if (AC_horizon_factor__mod__special>0.) {
      if (AC_k__mod__special[vertexIdx.x-NGHOST_VAL]>(a*h__mod__special*AC_horizon_factor__mod__special)) {
        tleff2=0.
        tldoteff2=0.
      }
    }
  }
  if (AC_llnk_spacing_adjustable__mod__special  ||  AC_llnk_spacing__mod__special) {
    trpsim=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*trpsi
    trpsikm=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*trpsi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)
    trpsidotm=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*trpsidot
    trdotpsim=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*trdotpsi
    trdoteff2km=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*trdoteff2*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)
    trdoteff2m=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*trdoteff2
    treff2km=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*treff2*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)
    treff2m=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*treff2
    grand__mod__special=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*(xi*h__mod__special-AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)*treff2*(+   AC_g__mod__special/(3.*(a*a)))/(twopi*twopi*twopi)
    grant__mod__special=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*(mq*h__mod__special-AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)*treff2*(-AC_lamf__mod__special/(2.*(a*a)))/(twopi*twopi*twopi)
    if (AC_lleft_psil_tl__mod__special) {
      tldoteff2km=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*tldoteff2*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)
      tldoteff2m=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*tldoteff2
      tleff2km=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*tleff2*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)
      tleff2m=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*tleff2
      grand__mod__special=grand__mod__special+(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*(xi*h__mod__special+AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)*tleff2*(+   AC_g__mod__special/(3.*(a*a)))/(twopi*twopi*twopi)
      grant__mod__special=grant__mod__special+(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*(mq*h__mod__special+AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)*tleff2*(-AC_lamf__mod__special/(2.*(a*a)))/(twopi*twopi*twopi)
    }
    if (AC_lconf_time__mod__special) {
      dgrant__mod__special=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*(-AC_lamf__mod__special/(2.*(a*a*a)))*(  (a*mq*(h__mod__special*h__mod__special)+AC_g__mod__special*qdot)*treff2+(mq*h__mod__special-AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)*2*trdoteff2  )/(twopi*twopi*twopi)
      if (AC_lleft_psil_tl__mod__special) {
        dgrant__mod__special=dgrant__mod__special+(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*(-AC_lamf__mod__special/(2.*(a*a*a)))*(  (a*mq*(h__mod__special*h__mod__special)+AC_g__mod__special*qdot)*tleff2+(mq*h__mod__special+AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)*2*tldoteff2  )/(twopi*twopi*twopi)
      }
    }
    else {
      dgrant__mod__special=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*(-AC_lamf__mod__special/(2.*(a*a*a)))*(  (a*mq*(h__mod__special*h__mod__special)+a*AC_g__mod__special*qdot)*treff2+(a*mq*h__mod__special-AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*2*trdoteff2  )/(twopi*twopi*twopi)
      if (AC_lleft_psil_tl__mod__special) {
        dgrant__mod__special=dgrant__mod__special+(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dlnk__mod__special)*(-AC_lamf__mod__special/(2.*(a*a*a)))*(  (a*mq*(h__mod__special*h__mod__special)+a*AC_g__mod__special*qdot)*tleff2+(a*mq*h__mod__special-AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*2*tldoteff2  )/(twopi*twopi*twopi)
      }
    }
  }
  else {
    trpsidotm=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__special)*trpsidot
    trdotpsim=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__special)*trdotpsi
    trdoteff2km=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__special)*trdoteff2*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)
    trdoteff2m=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__special)*trdoteff2
    treff2km=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__special)*treff2*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)
    treff2m=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__special)*treff2
    grand__mod__special=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__special)*(xi*h__mod__special-AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)*treff2*(+   AC_g__mod__special/(3.*(a*a)))/(twopi*twopi*twopi)
    grant__mod__special=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__special)*(mq*h__mod__special-AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)*treff2*(-AC_lamf__mod__special/(2.*(a*a)))/(twopi*twopi*twopi)
    if (AC_lconf_time__mod__special) {
      dgrant__mod__special=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__special)*(-AC_lamf__mod__special/(2.*(a*a*a)))*(  (a*mq*(h__mod__special*h__mod__special)+AC_g__mod__special*qdot)*treff2+(mq*h__mod__special-AC_k__mod__special[vertexIdx.x-NGHOST_VAL]/a)*2*trdoteff2  )/(twopi*twopi*twopi)
    }
    else {
      dgrant__mod__special=(4.*pi*(AC_k__mod__special[vertexIdx.x-NGHOST_VAL]*AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*AC_dk__mod__special)*(-AC_lamf__mod__special/(2.*(a*a*a)))*(  (a*mq*(h__mod__special*h__mod__special)+a*AC_g__mod__special*qdot)*treff2+(a*mq*h__mod__special-AC_k__mod__special[vertexIdx.x-NGHOST_VAL])*2*trdoteff2  )/(twopi*twopi*twopi)
    }
  }
  reduce_sum(grand__mod__special, AC_grand_sum)
  reduce_sum(dgrant__mod__special,AC_dgrant_sum)
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
Kernel calc_axion_integral(real t){suppress_unused_warning(t)}
#endif
