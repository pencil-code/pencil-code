function rtsafe

end


function pc_eoscalc,var1,var2,pp=pp,ee=ee, $
            lnrho_lntt=lnrho_lntt, lnrho_ss=lnrho_ss, $
            datadir=datadir,param=param, $
            PRINT=PRINT,QUIET=QUIET,HELP=HELP

  IF ( keyword_set(HELP) ) THEN BEGIN
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

  ; Allow param to be passed it if already loaded (eg. when called from inside another pc_ routine)
  if n_elements(param) eq 0 then pc_read_param,object=param,datadir=datadir
  
  

  result=0.

if (param.lionization) then begin

    if (keyword_set(lnrho_lnTT)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, lnTT
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, lnTT
      endelse

    endif else if (keyword_set(lnrho_ee)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ee
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ee
      endelse

    endif else then begin

      ;Assume lnrho_ss unless specified? or should we report an error?
      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
      endelse

    endelse

endif else if (param.lionization_fixed) then begin
    if (keyword_set(lnrho_lnTT)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
      endelse

    endif else if (keyword_set(lnrho_ee)) then begin

      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
      endelse

    endif else then begin

      ;Assume lnrho_ss unless specified?
      if keyword_set(pp) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
      endif else if keyword_set(ee) then begin
       ;result = fn of  where var1,var2 = lnrho, ss
      endelse

    endelse
endif else begin
   ; Make sure and missing cases say so...
   message,"pc_eoscalc: CASE NOT IMPLEMENTED"
endelse



  return,result
end

