;;
;;  $Id$
;;
;;  Calculate second derivative matrix f_l,ij.
;;
function derij,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
;  Default values.
;
  default, ghost, 0
;
  s=size(f)
;
  if (s[0] eq 4) then begin
;
    w=make_array(n_elements(f[*,0,0,0]),n_elements(f[0,*,0,0]),n_elements(f[0,0,*,0]),3,3,3)
    w[*,*,*,*,0,0]=xder2(f[*,*,*,*])
    w[*,*,*,*,1,1]=yder2(f[*,*,*,*])
    w[*,*,*,*,2,2]=zder2(f[*,*,*,*])
    w[*,*,*,*,0,1]=xderyder(f[*,*,*,*]) & w[*,*,*,*,1,0]=w[*,*,*,*,0,1]
    w[*,*,*,*,0,2]=xderzder(f[*,*,*,*]) & w[*,*,*,*,2,0]=w[*,*,*,*,0,2]
    w[*,*,*,*,1,2]=yderzder(f[*,*,*,*]) & w[*,*,*,*,2,1]=w[*,*,*,*,1,2]
;
  endif else if (s[0] eq 3) then begin
;
    w=make_array(n_elements(f[*,0,0]),n_elements(f[0,*,0]),n_elements(f[0,0,*]),3,3)
    w[*,*,*,0,0]=xder2(f[*,*,*])
    w[*,*,*,1,1]=yder2(f[*,*,*])
    w[*,*,*,2,2]=zder2(f[*,*,*])
    w[*,*,*,0,1]=xderyder(f[*,*,*]) & w[*,*,*,1,0]=w[*,*,*,0,1]
    w[*,*,*,0,2]=xderzder(f[*,*,*]) & w[*,*,*,2,0]=w[*,*,*,0,2]
    w[*,*,*,1,2]=yderzder(f[*,*,*]) & w[*,*,*,2,1]=w[*,*,*,1,2]
;
  endif else begin
    print, 'error: derij only implemented for 3-D and 4-D arrays'
  endelse
;
;  Set ghost zones.
;
  if (ghost) then w=pc_setghost(w,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, w
;
end
