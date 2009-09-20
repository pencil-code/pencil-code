;
; $Id$
;
;  Read ensure 'zero' and 'one' are set in the pc_precision common block.
;
pro pc_set_precision, precision=precision, dim=dim, datadir=datadir, QUIET=QUIET
COMPILE_OPT IDL2,HIDDEN
COMMON pc_precision, zero, one
;
if (n_elements(precision) eq 1) then begin
  if ((precision eq 'S') or (precision eq 's')) then begin
      one = 1.e0
  endif else if ((precision eq 'D') or (precision eq 'd')) then begin
      one = 1.D0
  endif else begin
      message, "precision = `"+ precision+"' makes no sense to me"
  endelse
  zero = 0*one
endif else begin
  if (n_elements(dim) ne 1) then $
      pc_read_dim, object=dim, datadir=datadir, /quiet
  precision=dim.precision
  if ((precision eq 'S') or (precision eq 's')) then begin
      one = 1.e0
  endif else if ((precision eq 'D') or (precision eq 'd')) then begin
      one = 1.D0
  endif else begin
      message, "precision = `"+precision+"' makes no sense to me"
  endelse
  zero = 0*one
endelse
;
end
