;;;;;;;;;;;;;;;;;;;;;;
;;;   minmax.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   08-Mar-2005
;;;
;;;  Description:
;;;   Minimum and maximum together
;;;  Usage:
;;;   mm = minmax(data)         ; only finite data
;;;   mm = minmax(data,/NAN)    ; if data can contain Infs (NaNs are
;;;                               ignored anyway)
;;;   Lx = minmax(x,/RANGE)     ; return max(x)-min(x)

function minmax, f, RANGE=range, _EXTRA=extra

  compile_opt idl2,hidden       ; suppress compilation message
  on_error,2                    ; return to caller if an error occurs

  if (n_elements(range) le 0) then range=0

  fmax = max(f, MIN=fmin, _EXTRA=extra)
  if (range) then begin
    return, fmax-fmin
  endif else begin
    return, [fmin, fmax]
  endelse
end
