;  $Id$
;
;  Calculate interstellar cooling term
;
;  Author: Antony Mee
;
function pc_interstellar_cool,lnrho=lnrho,ss=ss,lnTT=lnTT, $
    param=param,datadir=datadir

  default, datadir, 'data'
  if n_elements(param) eq 0 then pc_read_param,object=param,datadir=datadir
  pc_units,obj=units

  m_p_cgs=1.67262158D-24

  unit_Lambda = units.velocity^2 / units.density / units.time

  cooling_select=safe_get_tag(param,'cooling_select')
  coolingfunction_scalefactor=safe_get_tag(param,'coolingfunction_scalefactor')

  ;default,cooling_function_scalefactor,1.

  if (cooling_select eq 'RB') then begin
     print,'pc_interstellar_cool: RB cooling function'
     coolT_cgs=[ 100.D0,     2000.D0,    8000.D0,    1.D5,    4.D7,     1.D9 ]
     coolH_cgs=[ 2.2380D-32, 1.0012D-30, 4.6240D-36, 1.7800D-18, 3.2217D-27, 0.D0 ] / ( m_p_cgs )^2
     coolB=[ 2.,       1.5,      2.867,    -.65,    0.5,      0. ]
  endif else if (cooling_select eq 'SS') then begin
     ; These are the SS et al (2002) coefficients multiplied by m_proton**2 to obtain 
     ; same units as RB above
     print,'pc_interstellar_cool: SS cooling function'
     coolT_cgs=[ 10.D0,     141.D0,    313.D0,    6102.D0,    1.D5,     1.D9 ]
     coolH_cgs=[ 3.42D16, 9.10D18, 1.11D20, 2.00D8, 0.D0, 0.D0 ]
     coolB=[ 2.12,     1.0,      0.56,     3.67,    -.65 ,      0. ]
  endif else if (cooling_select eq 'off') then begin
     print,'pc_interstellar_cool: no cooling applied'
     coolT_cgs=dblarr(6)
     coolH_cgs=dblarr(6)
     coolB=fltarr(6)
  endif

  lncoolH = float(alog(coolH_cgs / unit_Lambda * (units.temperature^coolB) * coolingfunction_scalefactor)) 
  lncoolT = float(alog(coolT_cgs / units.temperature))

  default,lnTT,pc_eoscalc(lnrho,ss,/lnTT,/lnrho_ss)

help,unit_Lambda
help,coolH_cgs
help,units.temperature
print,lncoolT
  result=lnrho*0.
  for i=0,4 do begin
    coolwhere = where ((lnTT le lncoolT[i]) and (lnTT lt lncoolT[i+1])) 
    if max(coolwhere) ge 0 then $
      result[coolwhere]=result[coolwhere]+ $
        exp(lncoolH(i)+lnrho[coolwhere]+lnTT[coolwhere]*coolB(i))
  endfor

  return,result
end

