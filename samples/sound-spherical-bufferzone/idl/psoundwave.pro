;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   psoundwave.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   29-Apr-2002
;;;  $Id$
;;;
;;;  Description:
;;;   Plot diagnostics for the oblique sound wave problem.

ix = nx/2-3
iy = ny*0.7
iz = nz/2+5
;; Warn if we are in ghost zones
if ((ix lt nghx) or (ix gt nx-nghx-1)) then print, 'Warning: dangerous ix'
if ((iy lt nghy) or (iy gt ny-nghy-1)) then print, 'Warning: dangerous iy'
if ((iz lt nghz) or (iz gt nz-nghz-1)) then print, 'Warning: dangerous iz'

if (n_elements(lam0) le 0) then begin
  print, 'No variable lam0'
  if (t eq 0) then begin
    print, 'Setting lam0=lam, uu0=uu'
    lam0 = lam
    uu0 = uu
  endif else begin
    message, 'aborting'
  endelse
endif

s = texsyms()

save_state

!p.multi=[0,3,2]
!p.charsize=2
col = 120
!y.range=[-1,1]*0.003

plot, x, lam0[*,iy,iz], XTITLE='!8x!X', YTITLE='!6ln '+s.varrho
oplot,x, lam[*,iy,iz], COLOR=col

plot, y, lam0[ix,*,iz], XTITLE='!8y!X', YTITLE='!6ln '+s.varrho, $
    TITLE='!8t!6='+strtrim(t,2)
oplot, y, lam[ix,*,iz], COLOR=col

plot, z, lam0[ix,iy,*], XTITLE='!8z!X', YTITLE='!6ln '+s.varrho
oplot, z, lam[ix,iy,*], COLOR=col



plot, x, uu0[*,iy,iz,0], XTITLE='!8x!X', YTITLE='!8u!Dx!N!X', TITLE='!8u!Dx!N!6(!8x!6)!X'
oplot, x, uu[*,iy,iz,0], COLOR=col

plot, y, uu0[ix,*,iz,1], XTITLE='!8y!X', YTITLE='!8u!Dy!N!X', TITLE='!8u!Dy!N!6(!8y!6)!X'
oplot, y, uu[ix,*,iz,1], COLOR=col

plot, z, uu0[ix,iy,*,2], XTITLE='!8z!X', YTITLE='!8u!Dz!N!X', TITLE='!8u!Dz!N!6(!8z!6)!X'
oplot, z, uu[ix,iy,*,2], COLOR=col

restore_state

end
; End of file psoundwave.pro
