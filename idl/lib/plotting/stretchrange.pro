;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   stretchrange.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   25-Jun-2002
;;;
;;;  Description:
;;;   Stretch a given range by some amount
;;;  Intended use:
;;;     plot, x, y, YRANGE=stretchrange(minmax(y),0.1)

function stretchrange, range, fact, LOG=logar

  default, fact, 0.05
  default, logar, 0

  if (logar) then begin
    sgn = sign(range[0])
    r0 = alog(abs(range[0]))
    r1 = alog(abs(range[1]))
  endif else begin
    r0 = range[0]
    r1 = range[1]
  endelse

  rm = 0.5*(r0+r1)

  if (logar) then begin
    return, sgn * exp(rm + (1.+fact)*([r0,r1]-rm))
  endif else begin
    return, rm + (1+fact)*(range-rm)
  endelse

end
; End of file stretchrange.pro
