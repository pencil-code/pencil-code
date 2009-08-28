;;;;;;;;;;;;;;;;;;;;;
;;;   pvert.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   11-Nov-2001
;;;  $Id$
;;;
;;;  Description:
;;;    Plot vertical profiles of uz, lnrho and entropy.
;;;  PostScript output for the manual:
;;;    psa,file='pvert0.eps' & device,XSIZE=14,YSIZE=7
;;;    pvert_layout=[0,4,1]
;;;    .r thermo
;;;    .r pvert
;;;    pse

default, pvert_layout, [0,2,2]
default, nprofs, 10              ; set to N for only N profiles, 0 for all

sym = texsyms()

nign = 3                        ; Number of close-to-bdry points to ignore

;
; construct vector of vertical pencils to plot
;
Nxmax = mx-2*nign > 1
Nymax = my-2*nign > 1
Nmax = Nxmax*Nymax
if ((nprofs le 0) or (nprofs gt Nmax)) then nprofs = Nmax
ixp=[[mx/2]] & iyp=[[my/2]]     ; case nprofs=1
if (nprofs gt 1) then begin
  Nxp = sqrt(nprofs+1e-5)*Nxmax/Nymax > 1
  Nyp = sqrt(nprofs+1e-5)*Nymax/Nxmax > 1
  Nxyp = Nxp*Nyp
  ;; Correct for case where Nxp or Nyp was <<1 and thus Nxp*Nyp>>nprofs
  if (Nxp eq 1) then begin
    Nyp = Nyp*nprofs/Nxyp > 1
  endif else if (Ny eq 1) then begin
    Nxp = Nxp*nprofs/Nxyp > 1
  endif
  Nxyp = Nxp*Nyp
  ixp = nign + spread( indgen(Nxp)*Nxmax/Nxp, 1, Nyp )
  iyp = nign + spread( indgen(Nyp)*Nymax/Nyp, 0, Nxp )
  ixp = floor(ixp)
  iyp = floor(iyp)
; floor() will have dropped a degenerate trailing dimension (any IDL
; function will, this is crazy..), so let us reform again to be sure:
  ixp = reform(ixp,Nxp,Nyp,/OVERWRITE)
  iyp = reform(iyp,Nxp,Nyp,/OVERWRITE)
endif

save_state

!p.multi = pvert_layout

if (!d.name eq 'X') then begin
  red = 130/256.*!d.table_size   ; brick red for color table 5
endif else begin
  red = !p.color   ; black for PostScript  
endelse

if (!d.name eq 'X') then begin
  !p.charsize = 1. + (max(!p.multi)-1)*0.3
endif

!y.title = '!8z!X'

for ivar = 0,3 do begin

  case ivar of
    0: begin
      var = lnrho
      title = '!6ln '+sym.varrho
      if (ny eq 1) then xr = minmax(var[*,3,*]) else xr = minmax(var)
      if (n_elements(lnrhoinit) gt 0) then xr = minmax([xr,lnrhoinit])
    end
    1: begin
      var = uu[*,*,*,2]
      title = '!8u!Dz!N!X'
      if (ny eq 1) then xr = minmax(var[*,3,*]) else xr = minmax(var)
    end
    2: begin
      var = ss
      title = '!6Entropy !8s!X'
      if (ny eq 1) then xr = minmax(var[*,3,*]) else xr = minmax(var)
      if (n_elements(ssinit) gt 0) then xr = minmax([xr,ssinit])
    end
    3: begin
      var = cs0^2/gamma_m1*exp(gamma*ss+gamma_m1*(lnrho-lnrho0))
      title = '!6Temperature !8T!X'
      if (ny eq 1) then xr = minmax(var[*,3,*]) else xr = minmax(var)
      if (n_elements(Tinit) gt 0) then xr = minmax([xr,Tinit])
    end
  endcase

  plot, z, z, /NODATA, $
      XRANGE=xr, XSTYLE=3, $
      YRANGE=minmax(z), YSTYLE=3,  $
      TITLE=title
  for ix=0,(size(ixp))[1]-1 do begin
    for iy=0,(size(ixp))[2]-1 do begin
      oplot, var[ixp[ix,iy],iyp[ix,iy],*], z
    endfor
  endfor
  ophline, [z0, par.z1, par.z2, z0+Lz]
  if (ivar eq 1) then opvline

;; overplot initial profiles
  if (n_elements(Tinit) le 0) then begin
    if (ivar eq 0) then $
        message, 'No Tinit -- you should run thermo.pro', /INFO
  endif else begin
    case ivar of
      0: oplot, lnrhoinit, z, LINE=2, COLOR=red, THICK=2
      1: ;nothing to overplot
      2: oplot, ssinit, z, LINE=2, COLOR=red, THICK=2
      3: oplot, Tinit, z, LINE=2, COLOR=red, THICK=2
    endcase
  endelse

endfor

if (all(pvert_layout eq [0,2,2])) then begin
  ; don't know where to place otherwise
  xyouts, 0.45,0.5, '!8t!6=' + strtrim(t,2), /NORMAL
endif

restore_state

end
; End of file pvert.pro


