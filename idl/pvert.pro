;;;;;;;;;;;;;;;;;;;;;
;;;   pvert.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   11-Nov-2001
;;;
;;;  Description:
;;;   Plot vertical profiles of uz, lnrho and entropy.

default, pvert_layout, [0,2,2]
default, single, 0              ; set to one for plotting one single profile

nign = 3                        ; Number of close-to-bdry points to ignore

save_state

!p.multi = pvert_layout

if (!d.name eq 'X') then begin
  !p.charsize = 1. + (max(!p.multi)-1)*0.3
endif

!y.title = '!8z!X'

for ivar = 0,3 do begin

  case ivar of
    0: begin
      var = lam
      title = 'log density'
      xr = minmax(var)
      if (n_elements(laminit) gt 0) then xr = minmax([xr,laminit])
    end
    1: begin
      var = uu[*,*,*,2]
      title = '!8u!Dz!N!X'
      xr = minmax(var)
    end
    2: begin
      var = ent
      title = 'Entropy'
      xr = minmax(var)
      if (n_elements(entinit) gt 0) then xr = minmax([xr,entinit])
    end
    3: begin
      var = gamma/gamma1*exp(gamma*ent+gamma1*lam)
      title = 'Temperature'
      xr = minmax(var)
      if (n_elements(Tinit) gt 0) then xr = minmax([xr,Tinit])
    end
  endcase

  plot, z, z, /NODATA, $
      XRANGE=xr, XSTYLE=3, $
      YRANGE=minmax(z), YSTYLE=3,  $
      TITLE=title
  if (not single) then begin
    for ix=nign,nx-nign-1 do begin
      for iy=nign,ny-nign-1 do begin
        oplot, var[ix,iy,*], z
      endfor
    endfor
  endif else begin              ; short version for finit-sized PostScript
    oplot, var[nx/2,ny/2,*], z
  endelse
  ophline, [z0,z1,z2,ztop]
  if (ivar eq 1) then opvline

;; overplot initial profiles
  if (n_elements(Tinit) le 0) then begin
    if (ivar eq 0) then $
        message, 'No Tinit -- you should run thermo.pro', /INFO
  endif else begin
    case ivar of
      0: oplot, laminit, z, LINE=2, COLOR=130, THICK=2
      1: ;nothing to overplot
      2: oplot, entinit, z, LINE=2, COLOR=130, THICK=2
      3: oplot, Tinit, z, LINE=2, COLOR=130, THICK=2
    endcase
  endelse

endfor

restore_state

end
; End of file pvert.pro


