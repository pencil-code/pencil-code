;  $Id$
;
;  Calculate auxiliary dust variables such as distribution function f
;
;  Author: Anders Johansen
;
function pc_dust_aux,nd=nd,md=md,mi=mi,lnrho=lnrho,ss=ss,lncc=lncc,var=var, $
    param=param,datadir=datadir

  default, datadir, 'data'
  if n_elements(param) eq 0 then pc_read_param,object=param,datadir=datadir

  cmd='grep ndustspec '+datadir+'/index.pro'
  spawn, cmd, result
  res=execute(result[0])

  ldustdensity_log=safe_get_tag(param,'ldustdensity_log')
  if (ldustdensity_log) then begin
    nd_loc=exp(nd)
  endif else begin
    nd_loc=nd
  endelse

  result=0.

  deltamd = param.deltamd
  rhods = param.rhods

  unit_md = 1.
  if (param.dust_chemistry eq 'ice') then begin
    mumon = 18.
    mmon  = mumon*1.6733e-24
    unit_md = mmon
  endif

  if (var eq 'ad') then begin

    if (param.dust_geometry eq 'sphere') then begin
      result = (3*md*unit_md/(4.*!pi*param.rhods))^(1/3.)
    endif else begin
      print, 'pc_dust_aux: No valid dust geometry found in param'
    endelse

  endif else if (var eq 'md') then begin
 
    md00 = param.md0
    if (md00 eq 0.) then md00 = 4/3.*!pi*(param.ad0)^3*rhods/unit_md
    md=fltarr(ndustspec)
    mdplus=fltarr(ndustspec)
    mdminus=fltarr(ndustspec)
    for i=0,ndustspec-1 do begin
      mdminus[i] = md00*deltamd^i
      mdplus[i]  = md00*deltamd^(i+1)
      md[i] = 0.5*(mdplus[i]+mdminus[i])
    endfor
    result=md

  endif else if (var eq 'fd') then begin

    sized=size(nd_loc)
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
    
    result=nd_loc/delta

  endif else if (var eq 'rhod') then begin

    sized=size(nd_loc)
    sizemd=size(md)
    if (sized[0] ne sizemd[0]) then begin
      if (ndustspec eq 1) then begin
        result=nd_loc*md[0]
      endif else begin
        result=total(nd_loc*spread(md,indgen(sized[0]-1),sized[1:sized[0]]),sized[0])
      endelse
    endif else begin
      if (ndustspec eq 1) then begin
        result=nd_loc*md[0]
      endif else begin
        result=total(nd_loc*md,sized[0])
      endelse
    endelse

  endif else if (var eq 'ppmon') then begin

    pp = pc_eoscalc(lnrho,ss,datadir=datadir,/pp)
    mu = (1.+3.97153*param.xHe)/(1-param.xH2+param.xHe)
    result = pp*exp(lncc)*mu/mumon

  endif else if (var eq 'ppsat') then begin

    TT = pc_eoscalc(lnrho,ss,datadir=datadir,/tt)
    result = 6.035e12*exp(-5938./TT)

  endif else if (var eq 'smon') then begin

    pp = pc_eoscalc(lnrho,ss,datadir=datadir,/pp)
    TT = pc_eoscalc(lnrho,ss,datadir=datadir,/tt)
    mu = (1.+3.97153*param.xHe)/(1-param.xH2+param.xHe)
    result = pp*lncc*mu/mumon/(6.035e12*exp(-5938./TT))

  endif else if (var eq 'epsd') then begin

    sized=size(nd_loc)
    if (not defined(md)) then md=1.0
    sizemd=size(md)

    if (sized[0] ne sizemd[0]) then begin
      if (ndustspec eq 1) then begin
        result=md[0]*nd_loc*unit_md/exp(lnrho)
      endif else begin
        result=total(nd_loc*spread(md,indgen(sized[0]-1),sized[1:sized[0]]),sized[0])*unit_md/exp(lnrho)
      endelse
    endif else begin
      if (ndustspec eq 1) then begin
        result=nd_loc*md[0]*unit_md/exp(lnrho)
      endif else begin
        result=total(nd_loc*md,sized[0])*unit_md/exp(lnrho)
      endelse
    endelse

  endif else if (var eq 'unit_md') then begin

    result=unit_md

  endif else if (var eq 'mdave') then begin

    sized=size(nd_loc)
    sizemd=size(md)
    if (sized[0] ne sizemd[0]) then begin
      result=total(nd_loc*spread(md,indgen(sized[0]-1),sized[1:sized[0]]),sized[0])*unit_md/total(nd_loc,sized[0])
    endif else begin
      result=total(nd_loc*md,sized[0])*unit_md/total(nd_loc,sized[0])
    endelse

  endif else begin
    ; Make sure and missing cases say so...
    message,"pc_dust_aux: CASE NOT IMPLEMENTED"
  endelse

  return,result

end

