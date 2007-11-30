;
;  $Id: derij.pro,v 1.1 2007-11-30 12:24:49 ajohan Exp $
;
;  Calculate second derivative matrix f_l,ij.
;
function derij,f
  COMPILE_OPT IDL2,HIDDEN
;
  s=size(f)
;
  if (s[0] eq 4) then begin
;
    w=make_array(n_elements(f[*,0,0,0]),n_elements(f[0,*,0,0]),n_elements(f[0,0,*,0]),3,3,3,/nozero)
    w[*,*,*,*,0,0]=xder2(f[*,*,*,*])
    w[*,*,*,*,1,1]=yder2(f[*,*,*,*])
    w[*,*,*,*,2,2]=zder2(f[*,*,*,*])
    w[*,*,*,*,0,1]=xderyder(f[*,*,*,*]) & w[*,*,*,*,1,0]=w[*,*,*,*,0,1]
    w[*,*,*,*,0,2]=xderzder(f[*,*,*,*]) & w[*,*,*,*,2,0]=w[*,*,*,*,0,2]
    w[*,*,*,*,1,2]=yderzder(f[*,*,*,*]) & w[*,*,*,*,2,1]=w[*,*,*,*,1,2]
;
  endif else begin
    print, 'error: derij only implemented for 4-D arrays'
  endelse
;
  return, w
;
end
