function rtsafe

end


function pc_eoscalc,var1,var2,pp=pp,ee=ee,tt=tt,lntt=lntt,cs2=cs2, $
            lnrho_lntt=lnrho_lntt, lnrho_ss=lnrho_ss, $
            datadir=datadir,param=param,dim=dim, $
            PRINT=PRINT,QUIET=QUIET,HELP=HELP,ss=ss


  if ( keyword_set(HELP) ) then begin
    print, "Usage: "
    print, ""
    print, "result = pc_eoscalc(var1,var2,                                                "
    print, "             /PRESSURE,/ENERGY,                                               "
    print, "             datadir=datadir, proc=proc, /PRINT, /QUIET, /HELP)                              "
    print, "                                                                                            "
    PRint, "Calculate appropriately (for the current PC eos module in use) a thermodynamic quantity     "
    print, "from two quantities provided                                                                "
    print, "Returns zeros and empty in all variables on failure.                                        "
    print, "                                                                                            "
    print, "           If unspecified data is read for global calculation.                              "
    print, "                                                                                            "
    print, ""
    print, "  datadir: specify the root data directory. Default is './data'                     [string]"
    print, ""
    print, "   /PRINT: instruction to print summary of result to standard output                            "
    print, "   /QUIET: instruction not to print any 'helpful' information                               "
    print, "    /HELP: display this usage information, and exit                                         "
    return,0
  ENDIF

  pc_check_math
; default, datadir, 'data'
  IF (not keyword_set(datadir)) THEN datadir='data'
  ; Allow param to be passed it if already loaded (eg. when called from inside another pc_ routine)
  if n_elements(param) eq 0 then pc_read_param,object=param,datadir=datadir
  ;
  ;  Read the positions of variables in f
  ;  Can't just use `@data/pc_constants', as the data directory may have a different name
  ;
  cmd = 'perl -000 -ne '+"'"+'s/[ \t]+/ /g; print join(" & ",split(/\n/,$_)),"\n"'+"' "+datadir+'/pc_constants.pro'
  spawn, cmd, result
  res = flatten_strings(result)
  if (execute(res) ne 1) then $
    message, 'There was a problem with pc_constants.pro', /INFO
 
  result=0.

  lionization = safe_get_tag(param,'lionization',default=safe_get_tag(param,'leos_ionization',default=0)) 
  lionization_fixed = safe_get_tag(param,'lionization_fixed',default=safe_get_tag(param,'leos_ionization_fixed',default=0)) 

  if (not lionization and not lionization_fixed) then begin

    if (keyword_set(lnrho_lnTT)) then begin

      if keyword_set(pp) then begin
       if (param.ldensity_nolog)     then rho = var1 else rho = exp (var1)
       if (param.ltemperature_nolog) then TT = var2 else TT = exp (var2)
       gamma = param.gamma
       gamma_m1 = gamma - 1.
       cp = param.cp
       result = cp * (gamma - 1.0) / gamma * rho * TT
      endif else if keyword_set(ee) then begin
       if (param.ldensity_nolog)     then rho = var1 else rho = exp (var1)
       if (param.ltemperature_nolog) then TT = var2 else TT = exp (var2)
       gamma = param.gamma
       gamma_m1 = gamma - 1.
       cp = param.cp
       result = cp / gamma * rho * TT
      endif else if keyword_set(ss) then begin
       ;result = fn of  where var1,var2 = lnrho, lnTT
       if (param.ldensity_nolog)     then lnrho=alog(var1) else lnrho=var1
       if (param.ltemperature_nolog) then lnTT =alog(var2) else lnTT=var2 
      ;
       cs20=param.cs0^2
       lnrho0=alog(param.rho0)
       gamma=param.gamma
       gamma_m1=gamma-1.     
       cp=param.cp
       lnTT0=alog(cs20/(cp * gamma_m1))
       result=(lnTT-lnTT0)/gamma-gamma_m1/gamma*(lnrho-lnrho0)
      endif

    endif else if (keyword_set(lnrho_ee)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ee
       message,"Thermodynamic combination not implemented yet: /pp,/lnrho_ee"
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ee
       message,"Thermodynamic combination not implemented yet: /pp,/lnrho_ee"
      endif

    endif else begin ; Standard eos ideal gas with lnrho and ss

      if (param.ldensity_nolog) then lnrho=alog(var1) else lnrho=var1
      ss=var2

;     @data/pc_constants.pro
      cmd = 'perl -000 -ne '+"'"+'s/[ \t]+/ /g; print join(" & ",split(/\n/,$_)),"\n"'+"' "+datadir+'/pc_constants.pro'
      spawn, cmd, result
      res = flatten_strings(result)
      if (execute(res) ne 1) then $
        message, 'There was a problem with pc_constants.pro', /INFO
      impossible=3.9085e37
      pc_check_math,location='pc_eoscalc - ideal gas constants from code'
 
      pc_units,obj=units,param=param,dim=dim
;      if (param.mu eq impossible) then begin
;        mu_tmp=1.+2.97153*xHe
;      endif else begin
;        mu_tmp=param.mu
;      endelse
;      pc_check_math,location='pc_eoscalc - calc mu from composition'

      cs20=param.cs0^2
      lnrho0=alog(param.rho0)
      gamma=param.gamma
      gamma_m1=gamma-1.     
;      R_cgs=8.3144D7

      cp=param.cp
;      if (param.lcalc_cp) then begin
;        ;cp=k_B/(param.mu*m_H)
;        cp=float(R_cgs*gamma/(gamma_m1*mu_tmp)/units.velocity^2)
;        pc_check_math,location='pc_eoscalc - cp given lcalc_cp '
;      endif
      cp1=1./cp
      if (gamma ne 1.) then begin
        lnTT0=alog(cs20/(cp * (gamma-1)))
      endif else begin
        lnTT0=alog(cs20/cp)   
      endelse

      pc_check_math,location='pc_eoscalc - ideal gas coefficients'

      lnTT_=lnTT0+gamma*cp1*ss+gamma_m1*(lnrho-lnrho0)
     ; cs2=cs20*exp((gamma-1.)*(lnrho-lnrho0)+gamma*ss)

      if keyword_set(pp) then begin
       ; result=exp(lnrho)*cs2/gamma
        result=gamma_m1*cp/gamma*exp(lnTT_+lnrho)
        pc_check_math,location='pc_eoscalc - pp from lnrho and ss'
      endif else if keyword_set(ee) then begin
        ;result=cs2/(gamma*(gamma-1.))
        result=cp/gamma*exp(lnTT_)
        pc_check_math,location='pc_eoscalc - ee from lnrho and ss'
      endif else if keyword_set(tt) then begin
        ;result=cs20/(cp*(gamma-1.))*exp(gamma*ss+(gamma-1.)*(lnrho-lnrho0))
        result=exp(lnTT_)
        pc_check_math,location='pc_eoscalc - TT from lnrho and ss'
      endif else if keyword_set(lntt) then begin
        ;result=alog(cs20/(cp*(gamma-1.))*exp(gamma*ss+(gamma-1.)*(lnrho-lnrho0)))
        result=lnTT_
        pc_check_math,location='pc_eoscalc - lnTT from lnrho and ss'
      endif else if keyword_set(cs2) then begin
        result=cs20*exp(gamma*cp1*ss+(gamma-1)*(lnrho-lnrho0))
        pc_check_math,location='pc_eoscalc - cs2 from lnrho and ss'
      endif

    endelse

  endif else if (lionization_fixed) then begin

    yH0=param.yH0
    xHe=param.xHe
    xH2=param.xH2
    gamma=param.gamma
    nabla_ad=(gamma-1.)/gamma

    if (keyword_set(lnrho_lnTT)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
       message,"Thermodynamic combination not implemented yet!"
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
       message,"Thermodynamic combination not implemented yet!"
      endif

    endif else if (keyword_set(lnrho_ee)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
       message,"Thermodynamic combination not implemented yet!"
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
       message,"Thermodynamic combination not implemented yet!"
      endif

    endif else begin

      ;Assume lnrho_ss unless specified

      if (param.ldensity_nolog) then lnrho=alog(var1) else lnrho=var1
      ss=var2

;     @data/pc_constants.pro
      cmd = 'perl -000 -ne '+"'"+'s/[ \t]+/ /g; print join(" & ",split(/\n/,$_)),"\n"'+"' "+datadir+'/pc_constants.pro'
      spawn, cmd, result
      res = flatten_strings(result)
      if (execute(res) ne 1) then $
        message, 'There was a problem with pc_constants.pro', /INFO
      lnTT=lnTTss*ss+lnTTlnrho*lnrho+lnTT0
      TT=exp(lnTT)
    
      if keyword_set(pp) then begin
        result=(1.+yH0+xHe-xH2)*exp(lnrho)*TT*ss_ion
        pc_check_math,location='pc_eoscalc - pp from lnrho and ss'
      endif else if keyword_set(ee) then begin
        result=1.5*(1.+yH0+xHe-xH2)*ss_ion*TT+yH0*ss_ion*TT_ion
        pc_check_math,location='pc_eoscalc - ee from lnrho and ss'
      endif else if keyword_set(tt) then begin
        result=tt
        pc_check_math,location='pc_eoscalc - TT from lnrho and ss'
      endif else if keyword_set(lntt) then begin
        result=lnTT
        pc_check_math,location='pc_eoscalc - lnTT from lnrho and ss'
      endif else if keyword_set(cs2) then begin
        result=gamma*(1.+yH0+xHe-xH2)*ss_ion*TT
        pc_check_math,location='pc_eoscalc - cs2 from lnrho and ss'
      endif

    endelse
  endif else if (lionization) then begin

    if (keyword_set(lnrho_lnTT)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
       message,"Thermodynamic combination not implemented yet!"
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
       message,"Thermodynamic combination not implemented yet!"
      endif

    endif else if (keyword_set(lnrho_ee)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
       message,"Thermodynamic combination not implemented yet!"
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
       message,"Thermodynamic combination not implemented yet!"
      endif

    endif else begin

      ;Assume lnrho_ss unless specified

      if (param.ldensity_nolog) then lnrho=alog(var1) else lnrho=var1
      ss=var2

;     @data/pc_constants.pro
      cmd = 'perl -000 -ne '+"'"+'s/[ \t]+/ /g; print join(" & ",split(/\n/,$_)),"\n"'+"' "+datadir+'/pc_constants.pro'
      spawn, cmd, result
      res = flatten_strings(result)
      if (execute(res) ne 1) then $
        message, 'There was a problem with pc_constants.pro', /INFO
    
      if keyword_set(pp) then begin

      endif else if keyword_set(ee) then begin

      endif else if keyword_set(tt) then begin

      endif else if keyword_set(lntt) then begin

      endif else if keyword_set(cs2) then begin

      endif

    endelse
  endif else begin
    ; Make sure and missing cases say so...
    message,"pc_eoscalc: CASE NOT IMPLEMENTED"
  endelse

  return,result

end

