#if LBACKREACT_INFL
output global real AC_e2m_all__mod__backreact_infl
output global real AC_b2m_all__mod__backreact_infl
output global real AC_sige1m_all_nonaver__mod__backreact_infl
output global real AC_sigb1m_all_nonaver__mod__backreact_infl
output global real AC_a2rhom_all__mod__backreact_infl
output global real AC_a2rhopm_all__mod__backreact_infl
output global real AC_a2rhophim_all__mod__backreact_infl
output global real AC_a2rhogphim_all__mod__backreact_infl
output global real AC_edotbm_all__mod__backreact_infl
output global real AC_ddotam_all__mod__backreact_infl

field_order(AC_iinfl_phi__mod__backreact_infl-1)  Field F_INFL_PHI
field_order(AC_iinfl_dphi__mod__backreact_infl-1) Field F_INFL_DPHI

Kernel prep_ode_right(){
  real AC_a21__mod__backreact_infl
  real AC_a2__mod__backreact_infl
  if(AC_lflrw__mod__backreact_infl)
  {
    lnascale__mod__backreact_infl=AC_f_ode__mod__cdata[AC_iinfl_lna__mod__backreact_infl-1]
    ascale__mod__backreact_infl  =exp(lnascale__mod__backreact_infl)
    AC_a2__mod__backreact_infl      =ascale__mod__backreact_infl*ascale__mod__backreact_infl
    AC_a21__mod__backreact_infl     =1.0/AC_a2__mod__backreact_infl
  }
  real3 el
  real3 bb
  real3 gphi
  real e2
  real b2
  real gphi2
  real dphi
  real a2rhop
  real a2rho
  real ddota
  real phi
  real vpotential
  real edotb
  real sige1
  real sigb1
  real boost
  real gam_eb
  real eprime
  real bprime
  real jprime1
  phi=value(Field(AC_iinfl_phi__mod__backreact_infl-1))
  dphi=value(Field(AC_iinfl_dphi__mod__backreact_infl-1))
  real a2rhophim__mod__backreact_infl
  if (AC_lphi_hom__mod__disp_current) {
    a2rhop=(dphi*dphi)
    a2rho=0.5*(dphi*dphi)
    a2rhophim__mod__backreact_infl=a2rho
  }
  else {
    gphi = gradient(Field(AC_iinfl_phi__mod__backreact_infl-1))
    gphi2 = dot(gphi,gphi)
    a2rhogphim__mod__backreact_infl=0.5*gphi2
    reduce_sum(a2rhogphim__mod__backreact_infl/nwgrid,AC_a2rhogphim_all__mod__backreact_infl)
    a2rhop=(dphi*dphi)+onethird*gphi2
    a2rho=0.5*((dphi*dphi)+gphi2)
    a2rhophim__mod__backreact_infl=a2rho
  }
  if (AC_iex__mod__cdata!=0  &&  AC_lem_backreact__mod__backreact_infl) {
    el=value(F_EVEC)
    bb = curl((Field3){Field(AC_iaa__mod__cdata-1), Field(AC_iaa__mod__cdata), Field(AC_iaa__mod__cdata+1)})
    b2 = dot(bb,bb)
    e2 = dot(el,el)
    a2rhop=a2rhop+(0.5*fourthird)*(e2+b2)*AC_a21__mod__backreact_infl
    if (! AC_lphi_linear_regime__mod__disp_current) {
      a2rho=a2rho+0.5*(e2+b2)*AC_a21__mod__backreact_infl
    }
    if (AC_lrho_chi__mod__backreact_infl) {
      a2rho=a2rho+AC_scale_rho_chi_heqn__mod__backreact_infl*AC_a2__mod__backreact_infl*AC_f_ode__mod__cdata[AC_iinfl_rho_chi__mod__backreact_infl-1]
    }
  }
  real a2rhopm__mod__backreact_infl=a2rhop
  if(AC_enum_vprime_choice__mod__backreact_infl == enum_quadratic_string) {
    vpotential=0.5*AC_axionmass2__mod__backreact_infl*(phi*phi)
  }
  else if(AC_enum_vprime_choice__mod__backreact_infl == enum_quartic_string) {
    vpotential=AC_axionmass2__mod__backreact_infl*phi+(AC_lambda_axion__mod__backreact_infl/6.)*(phi*phi*phi)
  }
  else if(AC_enum_vprime_choice__mod__backreact_infl == enum_coszprofile_string) {
    vpotential=AC_axionmass2__mod__backreact_infl*AC_lambda_axion__mod__backreact_infl*sin(AC_lambda_axion__mod__backreact_infl*phi)
  }
  else {
  }
  if (AC_lphi_hom__mod__disp_current) {
    ddota=-(dphi*dphi)+4.*AC_a2__mod__backreact_infl*vpotential
  }
  else {
    ddota=-(dphi*dphi)-gphi2+4.*AC_a2__mod__backreact_infl*vpotential
  }
  reduce_sum(ddota*(four_pi_over_three/nwgrid),AC_ddotam_all__mod__backreact_infl)
  a2rho=a2rho+AC_a2__mod__backreact_infl*vpotential
  a2rhom__mod__backreact_infl=a2rho
  if (lmagnetic &&  AC_lem_backreact__mod__backreact_infl) {
    if (AC_lphi_hom__mod__disp_current  ||  AC_lrho_chi__mod__backreact_infl  ||  AC_lnoncollinear_eb__mod__disp_current  ||  AC_lnoncollinear_eb_aver__mod__disp_current   ||  AC_lcollinear_eb__mod__disp_current  ||  AC_lcollinear_eb_aver__mod__disp_current) {
      edotb = dot(el,bb)
      reduce_sum(edotb/nwgrid,AC_edotbm_all__mod__backreact_infl)
      if (AC_lnoncollinear_eb__mod__disp_current) {
        boost=sqrt(((e2-b2)*(e2-b2))+4.*(edotb*edotb))
        gam_eb=sqrt21*sqrt(1.+(e2+b2)/boost)
        eprime=sqrt21*sqrt(e2-b2+boost)
        bprime=sqrt21*sqrt(b2-e2+boost)*sign(1.,edotb)
        if (AC_lallow_bprime_zero__mod__disp_current) {
          if (eprime!=0.) {
            if (bprime!=0.) {
              jprime1=1./(6.*(pi*pi))*eprime*abs(bprime)/tanh(pi*abs(bprime)/eprime)
            }
            else {
              jprime1=1./(6.*(pi*pi*pi))*(eprime*eprime)
            }
            sige1=abs(jprime1)*eprime/(gam_eb*boost)
            sigb1=abs(jprime1)*edotb/(eprime*gam_eb*boost)
          }
          else {
            sige1=0.
            sigb1=0.
          }
        }
        else {
          if (eprime!=0.  &&  bprime!=0.) {
            jprime1=1./(6.*(pi*pi))*eprime*abs(bprime)/tanh(pi*abs(bprime)/eprime)
            sige1=abs(jprime1)*eprime/(gam_eb*boost)
            sigb1=abs(jprime1)*edotb/(eprime*gam_eb*boost)
          }
          else {
            sige1=0.
            sigb1=0.
          }
        }
      }
      if (AC_lcollinear_eb__mod__disp_current) {
        eprime=sqrt(e2)
        bprime=sqrt(b2)
        if (eprime!=0.  &&  bprime!=0.) {
          sige1=1./(6.*(pi*pi))*bprime/tanh(pi*bprime/eprime)
          sigb1=0.
        }
        else {
          sige1=0.
          sigb1=0.
        }
      }
    }
  }
  if ((lmagnetic && AC_lem_backreact__mod__backreact_infl)  &&  (AC_lrho_chi__mod__backreact_infl)) {
    if (AC_lnoncollinear_eb__mod__disp_current  ||  AC_lnoncollinear_eb_aver__mod__disp_current  ||   AC_lcollinear_eb__mod__disp_current  ||  AC_lcollinear_eb_aver__mod__disp_current) {
      e2m__mod__backreact_infl=e2
      b2m__mod__backreact_infl=b2
      if ((AC_lnoncollinear_eb__mod__disp_current  ||  AC_lcollinear_eb__mod__disp_current)) {
        sige1m__mod__backreact_infl=sige1
        sigb1m__mod__backreact_infl=sigb1
	reduce_sum(sige1/nwgrid,AC_sige1m_all_nonaver__mod__backreact_infl)
	reduce_sum(sigb1/nwgrid,AC_sigb1m_all_nonaver__mod__backreact_infl)
      }
      reduce_sum(e2m__mod__backreact_infl/nwgrid,AC_e2m_all__mod__backreact_infl)
      reduce_sum(b2m__mod__backreact_infl/nwgrid,AC_b2m_all__mod__backreact_infl)
    }
  }
  a2rhom__mod__backreact_infl    /= nwgrid
  a2rhopm__mod__backreact_infl   /= nwgrid
  a2rhophim__mod__backreact_infl /= nwgrid
  reduce_sum(a2rhom__mod__backreact_infl,AC_a2rhom_all__mod__backreact_infl)
  reduce_sum(a2rhopm__mod__backreact_infl,AC_a2rhopm_all__mod__backreact_infl)
  reduce_sum(a2rhophim__mod__backreact_infl,AC_a2rhophim_all__mod__backreact_infl)
}
#else

#if LKLEIN_GORDON
output global real AC_e2m_all__mod__klein_gordon
output global real AC_b2m_all__mod__klein_gordon
output global real AC_sige1m_all_nonaver__mod__klein_gordon
output global real AC_sigb1m_all_nonaver__mod__klein_gordon
output global real AC_a2rhom_all__mod__klein_gordon
output global real AC_a2rhopm_all__mod__klein_gordon
output global real AC_a2rhophim_all__mod__klein_gordon
output global real AC_a2rhogphim_all__mod__klein_gordon
output global real AC_edotbm_all__mod__klein_gordon
output global real AC_ddotam_all__mod__klein_gordon

field_order(AC_iphi__mod__klein_gordon-1)  Field F_PHI
field_order(AC_idphi__mod__klein_gordon-1) Field F_DPHI
#define F_PHI_UP_RE  F_PHI
#define F_DPHI_UP_RE F_DPHI
field_order(AC_iphi_up_im__mod__klein_gordon-1)   Field F_PHI_UP_IM
field_order(AC_iphi_down_re__mod__klein_gordon-1) Field F_PHI_DOWN_RE
field_order(AC_iphi_down_im__mod__klein_gordon-1) Field F_PHI_DOWN_IM

field_order(AC_idphi_up_im__mod__klein_gordon-1)   Field F_DPHI_UP_IM
field_order(AC_idphi_down_re__mod__klein_gordon-1) Field F_DPHI_DOWN_RE
field_order(AC_idphi_down_im__mod__klein_gordon-1) Field F_DPHI_DOWN_IM

Kernel prep_ode_right(){
  real AC_a21__mod__klein_gordon
  real AC_a2__mod__klein_gordon
  if(AC_lflrw__mod__klein_gordon)
  {
    lnascale__mod__klein_gordon=AC_f_ode__mod__cdata[AC_ilna__mod__klein_gordon-1]
    ascale__mod__klein_gordon  =exp(lnascale__mod__klein_gordon)
    AC_a2__mod__klein_gordon      =ascale__mod__klein_gordon*ascale__mod__klein_gordon
    AC_a21__mod__klein_gordon     =1.0/AC_a2__mod__klein_gordon
  }
  real3 el
  real3 bb
  real3 gphi
  real e2
  real b2
  real gphi2
  real dphi
  real a2rhop
  real a2rho
  real ddota
  real phi
  real vpotential
  real edotb
  real sige1
  real sigb1
  real boost
  real gam_eb
  real eprime
  real bprime
  real jprime1
  phi=value(Field(AC_iphi__mod__klein_gordon-1))
  dphi=value(Field(AC_idphi__mod__klein_gordon-1))
  real a2rhophim__mod__klein_gordon
  if (AC_lphi_hom__mod__disp_current) {
    a2rhop=(dphi*dphi)
    a2rho=0.5*(dphi*dphi)
    a2rhophim__mod__klein_gordon=a2rho
  }
  else {
    gphi = gradient(Field(AC_iphi__mod__klein_gordon-1))
    gphi2 = dot(gphi,gphi)
    a2rhogphim__mod__klein_gordon=0.5*gphi2
    reduce_sum(a2rhogphim__mod__klein_gordon/nwgrid,AC_a2rhogphim_all__mod__klein_gordon)
    a2rhop=(dphi*dphi)+onethird*gphi2
    a2rho=0.5*((dphi*dphi)+gphi2)
    a2rhophim__mod__klein_gordon=a2rho
  }
  if (AC_iex__mod__cdata!=0  &&  AC_lem_backreact__mod__klein_gordon) {
    el=value(F_EVEC)
    bb = curl((Field3){Field(AC_iaa__mod__cdata-1), Field(AC_iaa__mod__cdata), Field(AC_iaa__mod__cdata+1)})
    b2 = dot(bb,bb)
    e2 = dot(el,el)
    a2rhop=a2rhop+(0.5*fourthird)*(e2+b2)*AC_a21__mod__klein_gordon
    if (! AC_lphi_linear_regime__mod__disp_current) {
      a2rho=a2rho+0.5*(e2+b2)*AC_a21__mod__klein_gordon
    }
    if (AC_lrho_chi__mod__klein_gordon) {
      a2rho=a2rho+AC_scale_rho_chi_heqn__mod__klein_gordon*AC_a2__mod__klein_gordon*AC_f_ode__mod__cdata[AC_iinfl_rho_chi__mod__klein_gordon-1]
    }
  }
  real a2rhopm__mod__klein_gordon=a2rhop
  if(AC_enum_vprime_choice__mod__klein_gordon == enum_quadratic_string) {
    vpotential=0.5*AC_phimass2__mod__klein_gordon*(phi*phi)
  }
  else if(AC_enum_vprime_choice__mod__klein_gordon == enum_quartic_string) {
    vpotential=AC_phimass2__mod__klein_gordon*phi+(AC_lambda_phi__mod__klein_gordon/6.)*(phi*phi*phi)
  }
  else if(AC_enum_vprime_choice__mod__klein_gordon == enum_coszprofile_string) {
    vpotential=AC_phimass2__mod__klein_gordon*AC_lambda_phi__mod__klein_gordon*sin(AC_lambda_phi__mod__klein_gordon*phi)
  }
  else {
  }
  if (AC_lphi_hom__mod__disp_current) {
    ddota=-(dphi*dphi)+4.*AC_a2__mod__klein_gordon*vpotential
  }
  else {
    ddota=-(dphi*dphi)-gphi2+4.*AC_a2__mod__klein_gordon*vpotential
  }
  reduce_sum(ddota*(four_pi_over_three/nwgrid),AC_ddotam_all__mod__klein_gordon)
  a2rho=a2rho+AC_a2__mod__klein_gordon*vpotential
  a2rhom__mod__klein_gordon=a2rho
  if (lmagnetic &&  AC_lem_backreact__mod__klein_gordon) {
    if (AC_lphi_hom__mod__disp_current  ||  AC_lrho_chi__mod__klein_gordon  ||  AC_lnoncollinear_eb__mod__disp_current  ||  AC_lnoncollinear_eb_aver__mod__disp_current   ||  AC_lcollinear_eb__mod__disp_current  ||  AC_lcollinear_eb_aver__mod__disp_current) {
      edotb = dot(el,bb)
      reduce_sum(edotb/nwgrid,AC_edotbm_all__mod__klein_gordon)
      if (AC_lnoncollinear_eb__mod__disp_current) {
        boost=sqrt(((e2-b2)*(e2-b2))+4.*(edotb*edotb))
        gam_eb=sqrt21*sqrt(1.+(e2+b2)/boost)
        eprime=sqrt21*sqrt(e2-b2+boost)
        bprime=sqrt21*sqrt(b2-e2+boost)*sign(1.,edotb)
        if (AC_lallow_bprime_zero__mod__disp_current) {
          if (eprime!=0.) {
            if (bprime!=0.) {
              jprime1=1./(6.*(pi*pi))*eprime*abs(bprime)/tanh(pi*abs(bprime)/eprime)
            }
            else {
              jprime1=1./(6.*(pi*pi*pi))*(eprime*eprime)
            }
            sige1=abs(jprime1)*eprime/(gam_eb*boost)
            sigb1=abs(jprime1)*edotb/(eprime*gam_eb*boost)
          }
          else {
            sige1=0.
            sigb1=0.
          }
        }
        else {
          if (eprime!=0.  &&  bprime!=0.) {
            jprime1=1./(6.*(pi*pi))*eprime*abs(bprime)/tanh(pi*abs(bprime)/eprime)
            sige1=abs(jprime1)*eprime/(gam_eb*boost)
            sigb1=abs(jprime1)*edotb/(eprime*gam_eb*boost)
          }
          else {
            sige1=0.
            sigb1=0.
          }
        }
      }
      if (AC_lcollinear_eb__mod__disp_current) {
        eprime=sqrt(e2)
        bprime=sqrt(b2)
        if (eprime!=0.  &&  bprime!=0.) {
          sige1=1./(6.*(pi*pi))*bprime/tanh(pi*bprime/eprime)
          sigb1=0.
        }
        else {
          sige1=0.
          sigb1=0.
        }
      }
    }
  }
  if ((lmagnetic && AC_lem_backreact__mod__klein_gordon)  &&  (AC_lrho_chi__mod__klein_gordon)) {
    if (AC_lnoncollinear_eb__mod__disp_current  ||  AC_lnoncollinear_eb_aver__mod__disp_current  ||   AC_lcollinear_eb__mod__disp_current  ||  AC_lcollinear_eb_aver__mod__disp_current) {
      e2m__mod__klein_gordon=e2
      b2m__mod__klein_gordon=b2
      if ((AC_lnoncollinear_eb__mod__disp_current  ||  AC_lcollinear_eb__mod__disp_current)) {
        sige1m__mod__klein_gordon=sige1
        sigb1m__mod__klein_gordon=sigb1
	reduce_sum(sige1/nwgrid,AC_sige1m_all_nonaver__mod__klein_gordon)
	reduce_sum(sigb1/nwgrid,AC_sigb1m_all_nonaver__mod__klein_gordon)
      }
      reduce_sum(e2m__mod__klein_gordon/nwgrid,AC_e2m_all__mod__klein_gordon)
      reduce_sum(b2m__mod__klein_gordon/nwgrid,AC_b2m_all__mod__klein_gordon)
    }
  }
  a2rhom__mod__klein_gordon    /= nwgrid
  a2rhopm__mod__klein_gordon   /= nwgrid
  a2rhophim__mod__klein_gordon /= nwgrid
  reduce_sum(a2rhom__mod__klein_gordon,AC_a2rhom_all__mod__klein_gordon)
  reduce_sum(a2rhopm__mod__klein_gordon,AC_a2rhopm_all__mod__klein_gordon)
  reduce_sum(a2rhophim__mod__klein_gordon,AC_a2rhophim_all__mod__klein_gordon)
}

#else
Kernel prep_ode_right(){}
#endif
#endif
