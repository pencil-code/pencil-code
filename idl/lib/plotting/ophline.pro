;;;;;;;;;;;;;;;;;;;;;;;
;;;   ophline.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   20-Mar-2001
;;;  Version: 0.22
;;;  Description:
;;;   Overplot one ore more horizontal line(s).
;;;   Without argument, plot line y=0.

pro ophline, y, $
             LINESTYLE=linest , $
             _EXTRA=extra

  default, y, 0.
  default, linest, 1

  xrange = !x.crange
  if (!x.type eq 1) then xrange = 10^xrange ; (plotting with /xlog)

  for i=0,n_elements(y)-1 do begin
    oplot, xrange, [1,1]*y[i], LINESTYLE=linest, _EXTRA=extra
  endfor

end
; End of file ophline.pro
