; $Id: pc_set_precision.pro,v 1.2 2002-11-27 17:55:59 mee Exp $
;
;
;  Read ensure 'zero' and 'one' are set in the pc_precision common block.
;
pro pc_set_precision, precision=precision
  COMMON pc_precision, zero, one

if keyword_set(precision) then begin
    if ((precision eq 'S') or (precision eq 's')) then begin
        one = 1.e0
    endif else if ((precision eq 'D') or (precision eq 'd')) then begin
        one = 1.D0
    endif else begin
        message, "precision = `", precision, "' makes no sense to me"
    endelse
    zero = 0*one
endif else begin
    if not keyword_set(one) then begin
        pc_read_dim,precision=precision,/QUIET
        if ((precision eq 'S') or (precision eq 's')) then begin
            one = 1.e0
        endif else if ((precision eq 'D') or (precision eq 'd')) then begin
            one = 1.D0
        endif else begin
            message, "precision = `", precision, "' makes no sense to me"
        endelse
        zero = 0*one
    endif
endelse


end


