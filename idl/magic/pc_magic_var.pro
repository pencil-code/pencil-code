; $Id: pc_magic_var.pro,v 1.6 2004-07-10 19:02:41 brandenb Exp $
;
;  Purpose:
;
;  Example:
;
pro pc_magic_var,variables,tags,param=param,datadir=datadir

  default, datadir, 'data'
; Allow param to be passed it if already loaded (eg. when called from inside another pc_ routine)
  if n_elements(param) eq 0 then pc_read_param,object=param,datadir=datadir

  for iv=0,n_elements(variables)-1 do begin
    if variables[iv] eq 'bb' then begin
      variables[iv]='curl(aa)'
      tags[iv]='bb'
    endif else if variables[iv] eq 'divu' then begin
      variables[iv]='div(uu)'
      tags[iv]='divu'
    endif else if variables[iv] eq 'u2' then begin
      variables[iv]='dot2(uu)'
      tags[iv]='u2'
    endif else if variables[iv] eq 'cs2' then begin
      if (param.lionization and not param.lionization_fixed) then begin
        variables[iv]='pc_eoscalc(lnrho,lnTT,/cs2,/lnrho_lnTT)'
      endif else begin
        variables[iv]='pc_eoscalc(lnrho,ss,/cs2,/lnrho_ss)'
      endelse
      tags[iv]='cs2'
    endif else if variables[iv] eq 'ee' then begin
      if (param.lionization and not param.lionization_fixed) then begin
        variables[iv]='pc_eoscalc(lnrho,lnTT,/ee,/lnrho_lnTT)'
      endif else begin
        variables[iv]='pc_eoscalc(lnrho,ss,/ee,/lnrho_ss)'
      endelse
      tags[iv]='ee'
    endif else if variables[iv] eq 'tt' then begin
      if (param.lionization and not param.lionization_fixed) then begin
        variables[iv]='exp(lnTT)'
      endif else begin
        variables[iv]='pc_eoscalc(lnrho,ss,/tt,/lnrho_ss)'
      endelse
      tags[iv]='tt'
    endif else if variables[iv] eq 'pp' then begin
      if (param.lionization and not param.lionization_fixed) then begin
        variables[iv]='pc_eoscalc(lnrho,lnTT,/pp,/lnrho_lnTT)'
      endif else begin
        variables[iv]='pc_eoscalc(lnrho,ss,/pp,/lnrho_ss)'
      endelse
      tags[iv]='pp'
    endif else if variables[iv] eq 'rhod' then begin
      variables[iv]="pc_dust_aux(nd=nd,md=md,param=param,var='rhod')"
      tags[iv]='rhod'
    endif else if variables[iv] eq 'fd' then begin
      variables[iv]="pc_dust_aux(nd=nd,param=param,var='fd')"
      tags[iv]='fd'
    endif else if variables[iv] eq 'ad' then begin
      variables[iv]="pc_dust_aux(md=md,param=param,var='ad')"
      tags[iv]='ad'
    endif else if variables[iv] eq 'epsd' then begin
      variables[iv]="pc_dust_aux(lnrho=lnrho,nd=nd,md=md,par=param,var='epsd')"
      tags[iv]='epsd'
    endif else if variables[iv] eq 'smon' then begin
      variables[iv]="pc_dust_aux(lnrho=lnrho,ss=ss,nd=nd,md=md," + $
          "param=param,datadir=datadir,var='smon')"
      tags[iv]='smon'
    endif else if variables[iv] eq 'unit_md' then begin
      variables[iv]="pc_dust_aux(param=param,var='unit_md')"
      tags[iv]='unit_md'
    endif else if variables[iv] eq 'mdave' then begin
      variables[iv]="pc_dust_aux(nd=nd,md=md,param=param,var='mdave')"
      tags[iv]='mdave'
    endif
  endfor

end
