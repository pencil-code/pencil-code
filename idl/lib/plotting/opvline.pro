;;;;;;;;;;;;;;;;;;;;;;;
;;;   opvline.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   20-Mar-2001
;;;  Version: 0.22
;;;  Description:
;;;   Overplot one or more vertical line(s).
;;;   Without argument, plot line y=0.

pro opvline, x, $
             LINESTYLE=linest , $
             _EXTRA=extra

  default, x, 0.
  default, linest, 1

  yrange = !y.crange
  if (!y.type eq 1) then yrange = 10^yrange ; (plotting with /ylog)

  for i=0,n_elements(x)-1 do begin
    oplot, [1,1]*x[i], yrange, LINESTYLE=linest, _EXTRA=extra
  endfor

end
; End of file opvline.pro
