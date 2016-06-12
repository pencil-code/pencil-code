;;;;;;;;;;;;;;;;;;;;
;;;   rms.pro   ;;;
;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   15-Dec-2000
;;;  Version: 0.22
;;;  Description:
;;;   Calculate rms value of a field, allowing for the same arguments
;;;   as TOTAL().
;;;   Set keyword VECTOR to indicate that F is a vector field and you
;;;   want sqrt(mean(dot2(f))).

function rms, f, i, VECTOR=vect, DOUBLE=double

COMPILE_OPT IDL2,HIDDEN

  tmp = f
  if (keyword_set(vect)) then tmp = sqrt(dot2(abs(tmp)))
  if (n_elements(i) eq 0) then $
      return, sqrt(total(abs(tmp)^2,DOUBLE=double)/n_elements(tmp))
  ;
  s=size(tmp)
  return, sqrt(total(abs(tmp)^2,i,DOUBLE=double)/s(i))

end
; End of file rms.pro
