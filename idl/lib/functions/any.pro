;;;;;;;;;;;;;;;;;;;
;;;   any.pro   ;;;
;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   07-Feb-2001
;;;
;;;  Description:
;;;   Mimic the F90 ANY() function, i.e. return true if any of the elements
;;;   of an array is true, false otherwise.

function any, f
  ones = where(f ne 0)
  if (ones[0] ge 0) then res=1 else res=0

  return, res

end
; End of file any.pro
