;;;;;;;;;;;;;;;;;;;;
;;;   sign.pro   ;;;
;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   21-Jun-2001
;;;
;;;  Description:
;;;   The sign of a number -- either in mathematical (one argument) or
;;;   Fortran syntax (2 arguments).
;;;   In the first case, just an alias to my SGN function.

function sign, x, y

  if (n_elements(y) gt 0) then begin
    a = x
    x = y
  endif

  if ((size(x))[0] eq 0) then begin
    if (x eq 0) then begin
      res = 0*x
    endif else begin
      res = x/abs(x)
    endelse
  endif else begin
    good = where(x ne 0)
    res = make_array(SIZE=size(x),VALUE=0)
    if (good[0] ne -1) then begin
      res[good] = x[good]/abs(x[good])
    endif
  endelse

  if (n_elements(a) gt 0) then begin
    return, a*res
  endif else begin
    return, res
  endelse

end
; End of file sign.pro
