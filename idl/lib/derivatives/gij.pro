;
;  $Id: gij.pro,v 1.1 2004-12-25 14:30:31 ajohan Exp $
;
;  Calculate derivative matrix
;
function gij,f
  COMPILE_OPT IDL2,HIDDEN
;
  s=size(f)
;
  if (sf[0] eq 4) then begin
;
    w=make_array(n_elements(f[*,0,0,0]),n_elements(f[0,*,0,0]),n_elements(f[0,0,*,0]),3,3,/nozero)
    w[*,*,*,0,0]=xder(f[*,*,*,0])
    w[*,*,*,0,1]=yder(f[*,*,*,0])
    w[*,*,*,0,2]=zder(f[*,*,*,0])
    w[*,*,*,1,0]=xder(f[*,*,*,1])
    w[*,*,*,1,1]=yder(f[*,*,*,1])
    w[*,*,*,1,2]=zder(f[*,*,*,1])
    w[*,*,*,2,0]=xder(f[*,*,*,2])
    w[*,*,*,2,1]=yder(f[*,*,*,2])
    w[*,*,*,2,2]=zder(f[*,*,*,2])
;
  endif else begin
    print, 'error: gij only implemented for 4-D arrays'
  endelse
;
  return, w
;
end
