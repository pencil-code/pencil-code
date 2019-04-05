; 
;  $Id$
;
;  Calculate dot product of a matrix and a vector
;
function multm2,f
  COMPILE_OPT IDL2,HIDDEN
;
  sf=size(f)
;
  if (sf[0] eq 5) then begin
;
    w=0
    for j=0,2 do begin
    for i=0,2 do begin
      w = w + f[*,*,*,i,j]^2
    endfor
    endfor
;
  endif else begin
    print, 'error: multsv only implemented for f=5-D array'
  endelse
;  
  return, w
;
end
