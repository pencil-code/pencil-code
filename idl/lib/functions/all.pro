;;;;;;;;;;;;;;;;;;;
;;;   all.pro   ;;;
;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   07-Feb-2001
;;;
;;;  Description:
;;;   Mimic the F90 ALL() function, i.e. return true if all elements
;;;   of an array are true, false otherwise.

function all, f
  zeros = where(f eq 0)
  if (zeros[0] lt 0) then res=1 else res=0

  return, res

end
; End of file all.pro
