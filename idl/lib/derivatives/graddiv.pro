;
;  $Id: graddiv.pro,v 1.1 2007-12-09 17:08:53 ajohan Exp $
;
;  Calculate gradient of the divergence of a vector.
;
function graddiv,f
  COMPILE_OPT IDL2,HIDDEN
;
  s=size(f)
;
  if (s[0] eq 4) then begin
;
    w=make_array(n_elements(f[*,0,0,0]),n_elements(f[0,*,0,0]),n_elements(f[0,0,*,0]),3,/nozero)
    w[*,*,*,0]=   xder2(f[*,*,*,0])+xderyder(f[*,*,*,1])+xderzder(f[*,*,*,2])
    w[*,*,*,1]=xderyder(f[*,*,*,0])+   yder2(f[*,*,*,1])+yderzder(f[*,*,*,2])
    w[*,*,*,2]=xderzder(f[*,*,*,0])+yderzder(f[*,*,*,1])+   zder2(f[*,*,*,2])
;
  endif else begin
    print, 'error: graddiv only implemented for 4-D arrays'
  endelse
;
  return, w
;
end
