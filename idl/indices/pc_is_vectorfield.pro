;;
;;  $Id$
;;
function pc_is_vectorfield, varsize, $
    subscripts=subscripts, dim=dim, $
    TRIM=TRIM, NOVECTORS=NOVECTORS, YINYANG=yinyang, _EXTRA=e
COMPILE_OPT IDL2,HIDDEN
;
; Returns 1 for vector fields (or sets of vector fields as test methods variables)
;         2 for number of variables not a multiple of 3
;         0 otherwise
;
  if not keyword_set(yinyang) then yinyang=0 & yinyang=logical_true(yinyang)
;
  if (n_elements(dim) ne 1) then pc_read_dim, obj=dim, _EXTRA=e
;
;  The real number of dimensions for trimmed arrays is read from dim.dat.
;
  if (trim) then begin
    ndim = (dim.nx ne 1) + (dim.ny ne 1) + (dim.nz ne 1)
  endif else begin
    ndim=3
  endelse
;
;  Find out if array is a vector field.
;
  result=0
  if ( (varsize[0] eq 1) and (varsize[1] eq 1) ) then begin
;  Special: 0-D runs are represented in a 1-D array with one element.
    result=0
  endif else begin
    if ( varsize[0] gt ndim+yinyang ) then $
      if varsize[ndim+1] mod 3 eq 0 then result=1 else result=2
  endelse
;
;  Find proper subscripts in case of trimmed/non-trimmed variables.
;
  ;if ( (result eq 1) and (arg_present(subscripts)) ) then begin
  if ( arg_present(subscripts) ) then begin
    subscripts=make_array(varsize[0],/STRING,value='*')
    if (not trim) then begin
      subscripts[0]=string(dim.l1)+':'+string(dim.l2) 
      subscripts[1]=string(dim.m1)+':'+string(dim.m2) 
      subscripts[2]=string(dim.n1)+':'+string(dim.n2) 
    endif
  endif

  ;if ( (result eq 0) and (n_elements(subscripts) ne 0) ) then $
  ;    undefine, subscripts
  
  return, result

end
