;;;;;;;;;;;;;;;;;;;;;;;;
;;;   opcircle.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   30-Jun-2000
;;;  Version: 0.13
;;;  CVS $Revision: 1.2 $
;;;  Description:
;;;   Overplot circle(s) of given radius/radii
;;;  Usage:
;;;   opcircle, [r1,r2,r3]
;;;   opcircle, 0.5, NVERT=20     ; specify number of vertices
;;;   opcircle, [r1,r2,r3], XCENT=5., YCENT=-3 ; shift centre
;;;   opcircle, 0.5, FILLCOLOR=!p.color  ; fill circle with foreground colour
;;;   opcircle, 0.5, FILLCOLOR=!p.background  ; fill with background colour
;;;   opcircle, [0.5,1], FILLCOLOR=[30,60]  ; fill with different colours

pro opcircle, rv, $
              XCENT=xcent, YCENT=ycent, $
              NVERT=nvert, $
              FILLCOLOR=fillcolor, $
              NOCLIP=noclip, $
              _EXTRA=_extra

  default, nvert, 200
  default, noclip, 0
  default, xcent, 0.
  default, ycent, 0.
  ncirc = n_elements(rv)
  if (n_elements(xcent) eq 1) then xcent = replicate(xcent,ncirc)
  if (n_elements(ycent) eq 1) then ycent = replicate(ycent,ncirc)

  _phi = findgen(nvert+1)/nvert*2*!pi

  if (n_elements(fillcolor) le 0) then begin
    for i=0,n_elements(rv)-1 do begin
      oplot, xcent[i]+rv[i]*cos(_phi), ycent[i]+rv[i]*sin(_phi), $
          NOCLIP=noclip, _EXTRA=_extra
    endfor
  endif else begin
    ;; Fill circle; (different background colours not yet implemented)
    _phi = _phi[1:*]            ; closes automatically
    polyfill, rv[0]*cos(_phi), rv[0]*sin(_phi), COLOR=fillcolor, $
        NOCLIP=noclip
  endelse

end
; End of file opcircle.pro
