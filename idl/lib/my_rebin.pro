;;;;;;;;;;;;;;;;;;;;;;;;
;;;   my_rebin.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   23-Aug-2000
;;;
;;;  Description:
;;;   Like `rebin', but obtains new dimensions in one array and
;;;   accepts complex numbers.

function my_rebin, arr, dimv, _EXTRA=_extra
COMPILE_OPT IDL2,HIDDEN

  ;; Incredibly stupid: rebin does not work on complex data...
  s = size(arr)
  type = s[s[0]+1]
  if (type eq 6) then begin     ; complex
    return, complex(my_rebin(float    (arr),dimv,_EXTRA=_extra), $
                    my_rebin(imaginary(arr),dimv,_EXTRA=_extra))
  endif else if (type eq 9) then begin ; dcomplex
    return, complex(my_rebin(double   (arr),dimv,_EXTRA=_extra), $
                    my_rebin(imaginary(arr),dimv,_EXTRA=_extra))
  endif

  n = n_elements(dimv)
  if (n eq 0) then begin
    message, 'There are no empty arrays in IDL'
  endif else if (n eq 1) then begin
    return, rebin(arr, dimv[0], _EXTRA=_extra)
  endif else if (n eq 2) then begin
    return, rebin(arr, dimv[0], dimv[1], _EXTRA=_extra)
  endif else if (n eq 3) then begin
    return, rebin(arr, dimv[0], dimv[1], dimv[2], _EXTRA=_extra)
  endif else if (n eq 4) then begin
    return, rebin(arr, dimv[0], dimv[1], dimv[2], dimv[3], _EXTRA=_extra)
  endif else if (n eq 5) then begin
    return, rebin(arr, dimv[0], dimv[1], dimv[2], dimv[3], dimv[4], _EXTRA=_extra)
  endif else if (n eq 6) then begin
    return, rebin(arr, dimv[0], dimv[1], dimv[2], dimv[3], dimv[4], dimv[5], _EXTRA=_extra)
  endif else if (n eq 7) then begin
    return, rebin(arr, dimv[0], dimv[1], dimv[2], dimv[3], dimv[4], dimv[5], dimv[6], _EXTRA=_extra)
  endif else begin
    message,"You ask too much of me"
  endelse


end
; End of file my_rebin.pro
