;  $Id: pc_dust_aux.pro,v 1.6 2004-06-07 09:04:45 ajohan Exp $
;
;  Calculate auxiliary dust variables such as distribution function f
;
;  Author: Anders Johansen
;
function pc_dust_aux,nd=nd,md=md,cmd=cmd,mi=mi,ad=ad,fd=fd,cfd=cfd,rhod=rhod, $
    ppmon=ppmon,ppsat=ppsat,smon=smon,epsd=epsd,lnrho=lnrho,ss=ss,lncc=lncc,$
    unit_md=unit_md,param=param,datadir=datadir

  default, datadir, 'data'
  if n_elements(param) eq 0 then pc_read_param,object=param,datadir=datadir

  result=0.

  deltamd = param.deltamd
  rhods = param.rhods

  unit_md = 1.
  if (param.dust_chemistry eq 'ice') then begin
    mumon = 18.
    mmon  = mumon*1.6733e-24
    unit_md = mmon
  endif

  if (keyword_set(ad)) then begin

    if (param.dust_geometry eq 'sphere') then begin
      result = (3*md*unit_md/(4.*!pi*param.rhods))^(1/3.)
    endif else begin
      print, 'pc_dust_aux: No valid dust geometry found in param'
    endelse

  endif else if (keyword_set(cmd)) then begin
 
    sized=size(nd)
    ndustspec=sized(sized[0])
    md00 = param.md0
    if (md00 eq 0.) then md00 = 4/3.*!pi*(param.ad0)^3*rhods/unit_md
    md=fltarr(ndustspec)
    mdplus=fltarr(ndustspec)
    mdminus=fltarr(ndustspec)
    for i=0,ndustspec-1 do begin
      mdminus(i) = md00*deltamd^i
      mdplus(i)  = md00*deltamd^(i+1)
      md(i) = 0.5*(mdplus(i)+mdminus(i))
    endfor
    result=md

  endif else if (keyword_set(cfd)) then begin

    sized=size(nd)
    ndustspec=sized(sized[0])
    md00 = param.md0
    if (md00 eq 0.) then md00 = 4/3.*!pi*(param.ad0)^3*rhods/unit_md
    mdminus=fltarr(ndustspec)
    mdplus=fltarr(ndustspec)
    for i=0,ndustspec-1 do begin
      mdminus(i) = md00*deltamd^i
      mdplus(i)  = md00*deltamd^(i+1)
    endfor

    if (sized[0] ne 1) then begin
      delta=spread(mdplus-mdminus,indgen(sized[0]-1),sized[1:sized[0]])
    endif else begin
      delta = mdplus-mdminus
    endelse
    
    result=nd/delta

  endif else if (keyword_set(rhod)) then begin

    sized=size(nd)
    result=total(nd*md,sized[0])

  endif else if (keyword_set(ppmon)) then begin

    pp = pc_eoscalc(lnrho,ss,datadir=datadir,/pp)
    mu = (1.+3.97153*param.xHe)/(1-param.xH2+param.xHe)
    result = pp*exp(lncc)*mu/mumon

  endif else if (keyword_set(ppsat)) then begin

    TT = pc_eoscalc(lnrho,ss,datadir=datadir,/tt)
    result = 6.035e12*exp(-5938./TT)

  endif else if (keyword_set(smon)) then begin

    pp = pc_eoscalc(lnrho,ss,datadir=datadir,/pp)
    TT = pc_eoscalc(lnrho,ss,datadir=datadir,/tt)
    mu = (1.+3.97153*param.xHe)/(1-param.xH2+param.xHe)
    result = pp*lncc*mu/mumon/(6.035e12*exp(-5938./TT))

  endif else if (keyword_set(epsd)) then begin

    sized=size(nd)
    result = total(md*nd,sized[0])*unit_md/exp(lnrho)

  endif else if (keyword_set(unit_md)) then begin

    result=unit_md

  endif else begin
    ; Make sure and missing cases say so...
    message,"pc_dust_aux: CASE NOT IMPLEMENTED"
  endelse

  return,result

end

