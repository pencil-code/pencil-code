function rtsafe

end


function pc_eoscalc,var1,var2,pp=pp,ee=ee,tt=tt,cs2=cs2, $
            lnrho_lntt=lnrho_lntt, lnrho_ss=lnrho_ss, $
            datadir=datadir,param=param, $
            PRINT=PRINT,QUIET=QUIET,HELP=HELP

  if ( keyword_set(HELP) ) then begin
    print, "Usage: "
    print, ""
    print, "result = pc_eoscalc(var1,var2,                                                "
    print, "             /PRESSURE,/ENERGY,                                               "
    print, "             datadir=datadir, proc=proc, /PRINT, /QUIET, /HELP)                              "
    print, "                                                                                            "
    print, "Calculate appropriately (for the current PC eos module in use) a thermodynamic quantity     "
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

  default, datadir, 'data'
  ; Allow param to be passed it if already loaded (eg. when called from inside another pc_ routine)
  if n_elements(param) eq 0 then pc_read_param,object=param,datadir=datadir
  
  result=0.

  if (not param.lionization and not param.lionization_fixed) then begin

    if (keyword_set(lnrho_lnTT)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, lnTT
       message,"Thermodynamic combination not implemented yet: /pp,/lnrho_lnTT"
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, lnTT
       message,"Thermodynamic combination not implemented yet: /pp,/lnrho_lnTT"
      endif

    endif else if (keyword_set(lnrho_ee)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ee
       message,"Thermodynamic combination not implemented yet: /pp,/lnrho_ee"
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ee
       message,"Thermodynamic combination not implemented yet: /pp,/lnrho_ee"
      endif

    endif else begin

      lnrho=var1
      ss=var2

      if (param.lcalc_cp) then cp=k_B/(mu*m_H) else cp=1.
      cs20=param.cs20
      lnrho0=param.lnrho0
      
      cs2=cs20*exp((gamma-1.)*(lnrho-lnrho0)+gamma*ss)

      if keyword_set(pp) then begin
        result=exp(lnrho)*cs2/gamma
      endif else if keyword_set(ee) then begin
        result=cs2/(gamma-1.)
      endif else if keyword_set(tt) then begin
        result=cs20/(cp*(gamma-1.))*exp(gamma*ss+(gamma-1.)*(lnrho-lnrho0))
      endif else if keyword_set(cs2) then begin
        result=cs2
      endif

    endelse

  endif else if (param.lionization_fixed) then begin

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

      lnrho=var1
      ss=var2

      @data/pc_constants.pro
      lnTT=lnTTss*ss+lnTTlnrho*lnrho+lnTT0
      TT=exp(lnTT)
    
      if keyword_set(pp) then begin
        result=(1.+yH0+xHe-xH2)*exp(lnrho)*TT*ss_ion
      endif else if keyword_set(ee) then begin
        result=1.5*(1.+yH0+xHe-xH2)*ss_ion*TT+yH0*ss_ion*TT_ion
      endif else if keyword_set(tt) then begin
        result=tt
      endif else if keyword_set(cs2) then begin
        result=gamma*(1.+yH0+xHe-xH2)*ss_ion*TT
      endif

    endelse
  endif else if (param.lionization) then begin

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

      lnrho=var1
      ss=var2

      @data/pc_constants.pro
    
      if keyword_set(pp) then begin

      endif else if keyword_set(ee) then begin

      endif else if keyword_set(tt) then begin

      endif else if keyword_set(cs2) then begin

      endif

    endelse
  endif else begin
    ; Make sure and missing cases say so...
    message,"pc_eoscalc: CASE NOT IMPLEMENTED"
  endelse

  return,result

end

