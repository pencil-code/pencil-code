;;
;;  $Id$
;;
;;  Second derivative d2f/dxdy
;;  - 6th-order
;;  - with ghost cells
;;
function xderyder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z
  common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
;
;  Default values.
;
  default, ghost, 0
;
;  Assume nghost=3 for now.
;
  default, nghost, 3
;
;  Calculate mx, my, and mz, based on the input array size
;
  s=size(f) & d=make_array(size=s)
  mx=s[1] & my=s[2] & mz=s[3]
;
;  Determine location of ghost zones
;
  l1=nghost & l2=mx-nghost-1
  m1=nghost & m2=my-nghost-1
  n1=nghost & n2=mz-nghost-1
;
;  Calculate d2f/dxdy.
;
  if (s[0] eq 3) then begin
    if (not ldegenerated[0] and not ldegenerated[1]) then begin
      for n=n1,n2 do begin & for m=m1,m2 do begin
        ;
        ;  take care of nonuniform mesh
        ;
        if (lequidist[0]) then begin
          fac=1/(60.^2*(x[4]-x[3])*(y[4]-y[3]))
        endif else begin
          fac=(1/60.^2)*dx_1[l1:l2]*dy_1[m]
        endelse
        ;
        d[l1:l2,m,n]=fac*( $
            45.*( (45.*(f[l1+1:l2+1,m+1,n]-f[l1-1:l2-1,m+1,n])  $
                   -9.*(f[l1+2:l2+2,m+1,n]-f[l1-2:l2-2,m+1,n])  $
                      +(f[l1+3:l2+3,m+1,n]-f[l1-3:l2-3,m+1,n])) $
                 -(45.*(f[l1+1:l2+1,m-1,n]-f[l1-1:l2-1,m-1,n])  $
                   -9.*(f[l1+2:l2+2,m-1,n]-f[l1-2:l2-2,m-1,n])  $
                      +(f[l1+3:l2+3,m-1,n]-f[l1-3:l2-3,m-1,n])))$
            -9.*( (45.*(f[l1+1:l2+1,m+2,n]-f[l1-1:l2-1,m+2,n])  $
                   -9.*(f[l1+2:l2+2,m+2,n]-f[l1-2:l2-2,m+2,n])  $
                      +(f[l1+3:l2+3,m+2,n]-f[l1-3:l2-3,m+2,n])) $
                 -(45.*(f[l1+1:l2+1,m-2,n]-f[l1-1:l2-1,m-2,n])  $
                   -9.*(f[l1+2:l2+2,m-2,n]-f[l1-2:l2-2,m-2,n])  $
                      +(f[l1+3:l2+3,m-2,n]-f[l1-3:l2-3,m-2,n])))$
               +( (45.*(f[l1+1:l2+1,m+3,n]-f[l1-1:l2-1,m+3,n])  $
                   -9.*(f[l1+2:l2+2,m+3,n]-f[l1-2:l2-2,m+3,n])  $
                      +(f[l1+3:l2+3,m+3,n]-f[l1-3:l2-3,m+3,n])) $
                 -(45.*(f[l1+1:l2+1,m-3,n]-f[l1-1:l2-1,m-3,n])  $
                   -9.*(f[l1+2:l2+2,m-3,n]-f[l1-2:l2-2,m-3,n])  $
                      +(f[l1+3:l2+3,m-3,n]-f[l1-3:l2-3,m-3,n]))) )
      endfor & endfor
    endif
;
  endif else if (s[0] eq 4) then begin
;
    if (not ldegenerated[0] and not ldegenerated[1]) then begin
      for n=n1,n2 do begin & for m=m1,m2 do begin
        ;
        ;  take care of nonuniform mesh
        ;
        if (lequidist[0]) then begin
          fac=1/(60.^2*(x[4]-x[3])*(y[4]-y[3]))
        endif else begin
          fac=(1/60.^2)*dx_1[l1:l2]*dy_1[m]
        endelse
        ;
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
    endif
  endif else begin
    print, 'error: xderyder_6th_ghost not implemented for ', $
        strtrim(s[0],2), '-D arrays'
  endelse
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
