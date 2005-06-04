;;;;;;;;;;;;;;;;;;;;;;
;;;   spread.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   23-Aug-2000
;;;  Version: 0.3
;;;
;;;  Description:
;;;    Mimic the Fortran function `spread'.
;;;  Usage:
;;;      res = spread([arr], dim, size)
;;;      res = spread([arr], [dims],[sizes])
;;;    Inserts a new DIMth dimension of size SIZE or does so
;;;    consecutively with DIMS and SIZES.
;;;  Examples:
;;;    xx = spread(x, [1,2], [ny,nz])
;;;    yy = spread(y, [0,2], [nx,nz])
;;;    zz = spread(z, [0,1], [nx,ny])

function spread, arr, dim, size 
COMPILE_OPT IDL2,HIDDEN

;  if (n_elements(arr) eq 1) then arr = [arr] ; Treat scalar as 1-element array

  s = size(arr)
  N = n_elements(dim)

  if (N ge 2) then begin
    ;; Recursively process the arrays DIM and SIZE
;    res = spread(arr,dim[0:N-2],size[0:N-2])
;    res = spread(res,dim[N-1],size[N-1])
    res = spread(spread(arr,dim[0:N-2],size[0:N-2]),dim[N-1],size[N-1])
    return, res
  endif
  ;; Scalar arguments
  dim = dim[0]                  ; We want scalars
  size = size[0]
  if (dim gt s[0]) then message,"Array has too few dimensions"
  if (s[0] eq 0) then begin
    oldimv = 1
  endif else begin
    oldimv = s[1:s[0]]
  endelse
  if (dim eq 0) then begin
    newdimf = [1, oldimv]
    newdimb = [size, oldimv]
  endif else if (dim eq s[0]) then begin
    newdimf = [oldimv, 1]
    newdimb = [oldimv, size]
  endif else begin
    newdimf = [oldimv[0:dim-1], 1, oldimv[dim:*]]    
    newdimb = [oldimv[0:dim-1], size, oldimv[dim:*]]    
  endelse
  if (s[0] eq 0) then begin
    res = my_rebin([arr],size)
  endif else begin
    res = my_rebin(reform(arr, newdimf), newdimb)
    ;; Any trailing component of newdimb equal to 1 will be skipped by
    ;; rebin, so just to be safe, we call reform:
    res = reform(res, newdimb)
  endelse
  return, res
end
; End of file spread.pro
