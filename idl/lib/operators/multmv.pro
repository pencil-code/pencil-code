; 
;  $Id$
;
;  Calculate dot product of a matrix and a vector
;
function multmv,f,g
  COMPILE_OPT IDL2,HIDDEN
;
  sf=size(f)
  sg=size(g)
;
  if (sf[0] eq 5 and sg[0] eq 4) then begin
;
    w=make_array(size=sg,/nozero)
    for j=0,2 do begin
      w[*,*,*,0] = w[*,*,*,0] + f[*,*,*,0,j]*g[*,*,*,j]
      w[*,*,*,1] = w[*,*,*,1] + f[*,*,*,1,j]*g[*,*,*,j]
      w[*,*,*,2] = w[*,*,*,2] + f[*,*,*,2,j]*g[*,*,*,j]
    endfor
;
  endif else begin
    print, 'error: multsv only implemented for f=5-D array and g=4-D array'
  endelse
;  
  return, w
;
end
