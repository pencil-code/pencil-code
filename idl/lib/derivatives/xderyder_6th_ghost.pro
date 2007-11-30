;
;  $Id: xderyder_6th_ghost.pro,v 1.1 2007-11-30 12:24:49 ajohan Exp $
;
;  Second derivative d2f/dxdy
;  - 6th-order
;  - with ghost cells
;
;***********************************************************************
function xderyder,f
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
;
;  calculate mx, my, and mz, based on the input array size
;
  s=size(f) & d=make_array(size=s)
  mx=s[1] & my=s[2] & mz=s[3]
;
;  Not implemented for non-equidistant grid.
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (not all(lequidist)) then begin
    print, 'xderyder is not implemented for non equidistant grid'
    stop
  endif
;
;  Determine location of ghost zones - assume nghost=3 for now.
;
  l1=3 & l2=mx-4
  m1=3 & m2=my-4
  n1=3 & n2=mz-4
;
  dx_1=1/(x[4]-x[3])
  dy_1=1/(y[4]-y[3])
;
;  Calculate d2f/dxdy.
;
  if (s[0] eq 4) then begin
    d=fltarr(mx,my,mz,3)
    for n=n1,n2 do begin & for m=m1,m2 do begin
      fac=(1/60.0^2)*dx_1*dy_1
      d[l1:l2,m,n,*]=fac*( $
          45.*( (45.*(f[l1+1:l2+1,m+1,n,*]-f[l1-1:l2-1,m+1,n,*])  $
                 -9.*(f[l1+2:l2+2,m+1,n,*]-f[l1-2:l2-2,m+1,n,*])  $
                    +(f[l1+3:l2+3,m+1,n,*]-f[l1-3:l2-3,m+1,n,*])) $
               -(45.*(f[l1+1:l2+1,m-1,n,*]-f[l1-1:l2-1,m-1,n,*])  $
                 -9.*(f[l1+2:l2+2,m-1,n,*]-f[l1-2:l2-2,m-1,n,*])  $
                    +(f[l1+3:l2+3,m-1,n,*]-f[l1-3:l2-3,m-1,n,*])))$
          -9.*( (45.*(f[l1+1:l2+1,m+2,n,*]-f[l1-1:l2-1,m+2,n,*])  $
                 -9.*(f[l1+2:l2+2,m+2,n,*]-f[l1-2:l2-2,m+2,n,*])  $
                    +(f[l1+3:l2+3,m+2,n,*]-f[l1-3:l2-3,m+2,n,*])) $
               -(45.*(f[l1+1:l2+1,m-2,n,*]-f[l1-1:l2-1,m-2,n,*])  $
                 -9.*(f[l1+2:l2+2,m-2,n,*]-f[l1-2:l2-2,m-2,n,*])  $
                    +(f[l1+3:l2+3,m-2,n,*]-f[l1-3:l2-3,m-2,n,*])))$
             +( (45.*(f[l1+1:l2+1,m+3,n,*]-f[l1-1:l2-1,m+3,n,*])  $
                 -9.*(f[l1+2:l2+2,m+3,n,*]-f[l1-2:l2-2,m+3,n,*])  $
                    +(f[l1+3:l2+3,m+3,n,*]-f[l1-3:l2-3,m+3,n,*])) $
               -(45.*(f[l1+1:l2+1,m-3,n,*]-f[l1-1:l2-1,m-3,n,*])  $
                 -9.*(f[l1+2:l2+2,m-3,n,*]-f[l1-2:l2-2,m-3,n,*])  $
                    +(f[l1+3:l2+3,m-3,n,*]-f[l1-3:l2-3,m-3,n,*]))) )
    endfor & endfor
  endif else begin
    print, 'error: xderyder_6th_ghost not implemented for ', $
        strtrim(s[0],2), '-D arrays'
  endelse
;
  return, d
;
end
