;;
;;  $Id$
;;
;;  Second derivative d2/dxdz
;;  - 6th-order
;;  - with ghost cells
;;
function xderzder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
;  Calculate mx, my, and mz, based on the input array size.
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
;  Calculate d2f/dxdz.
;
  if (s[0] eq 3) then begin
    if (not ldegenerated[1] and not ldegenerated[2]) then begin
      for n=n1,n2 do begin & for m=m1,m2 do begin
        ;
        ;  take care of nonuniform mesh in the x-direction
        ;
        if (lequidist[0]) then begin
          fac=1/60.^2*dx_1[l1]*dz_1[n]
        endif else begin
          fac=1/60.^2*dx_1[l1:l2]*dz_1[n]
        endelse
        ;
        d[l1:l2,m,n]=fac*( $
            45.*( (45.*(f[l1+1:l2+1,m,n+1]-f[l1-1:l2-1,m,n+1])  $
                   -9.*(f[l1+2:l2+2,m,n+1]-f[l1-2:l2-2,m,n+1])  $
                      +(f[l1+3:l2+3,m,n+1]-f[l1-3:l2-3,m,n+1])) $
                 -(45.*(f[l1+1:l2+1,m,n-1]-f[l1-1:l2-1,m,n-1])  $
                   -9.*(f[l1+2:l2+2,m,n-1]-f[l1-2:l2-2,m,n-1])  $
                      +(f[l1+3:l2+3,m,n-1]-f[l1-3:l2-3,m,n-1])))$
            -9.*( (45.*(f[l1+1:l2+1,m,n+2]-f[l1-1:l2-1,m,n+2])  $
                   -9.*(f[l1+2:l2+2,m,n+2]-f[l1-2:l2-2,m,n+2])  $
                      +(f[l1+3:l2+3,m,n+2]-f[l1-3:l2-3,m,n+2])) $
                 -(45.*(f[l1+1:l2+1,m,n-2]-f[l1-1:l2-1,m,n-2])  $
                   -9.*(f[l1+2:l2+2,m,n-2]-f[l1-2:l2-2,m,n-2])  $
                      +(f[l1+3:l2+3,m,n-2]-f[l1-3:l2-3,m,n-2])))$
               +( (45.*(f[l1+1:l2+1,m,n+3]-f[l1-1:l2-1,m,n+3])  $
                   -9.*(f[l1+2:l2+2,m,n+3]-f[l1-2:l2-2,m,n+3])  $
                      +(f[l1+3:l2+3,m,n+3]-f[l1-3:l2-3,m,n+3])) $
                 -(45.*(f[l1+1:l2+1,m,n-3]-f[l1-1:l2-1,m,n-3])  $
                   -9.*(f[l1+2:l2+2,m,n-3]-f[l1-2:l2-2,m,n-3])  $
                      +(f[l1+3:l2+3,m,n-3]-f[l1-3:l2-3,m,n-3]))) )
      endfor & endfor
    endif
;
  endif else if (s[0] eq 4) then begin
;
    if (not ldegenerated[0] and not ldegenerated[2]) then begin
      for n=n1,n2 do begin & for m=m1,m2 do begin
        ;
        ;  take care of nonuniform mesh
        ;
        if (lequidist[0]) then begin
          fac=1/60.^2*dx_1[l1]*dz_1[n]
        endif else begin
          fac=1/60.^2*dx_1[l1:l2]*dz_1[n]
        endelse
        ;
        d[l1:l2,m,n,*]=fac*( $
            45.*( (45.*(f[l1+1:l2+1,m,n+1,*]-f[l1-1:l2-1,m,n+1,*])  $
                   -9.*(f[l1+2:l2+2,m,n+1,*]-f[l1-2:l2-2,m,n+1,*])  $
                      +(f[l1+3:l2+3,m,n+1,*]-f[l1-3:l2-3,m,n+1,*])) $
                 -(45.*(f[l1+1:l2+1,m,n-1,*]-f[l1-1:l2-1,m,n-1,*])  $
                   -9.*(f[l1+2:l2+2,m,n-1,*]-f[l1-2:l2-2,m,n-1,*])  $
                      +(f[l1+3:l2+3,m,n-1,*]-f[l1-3:l2-3,m,n-1,*])))$
            -9.*( (45.*(f[l1+1:l2+1,m,n+2,*]-f[l1-1:l2-1,m,n+2,*])  $
                   -9.*(f[l1+2:l2+2,m,n+2,*]-f[l1-2:l2-2,m,n+2,*])  $
                      +(f[l1+3:l2+3,m,n+2,*]-f[l1-3:l2-3,m,n+2,*])) $
                 -(45.*(f[l1+1:l2+1,m,n-2,*]-f[l1-1:l2-1,m,n-2,*])  $
                   -9.*(f[l1+2:l2+2,m,n-2,*]-f[l1-2:l2-2,m,n-2,*])  $
                      +(f[l1+3:l2+3,m,n-2,*]-f[l1-3:l2-3,m,n-2,*])))$
               +( (45.*(f[l1+1:l2+1,m,n+3,*]-f[l1-1:l2-1,m,n+3,*])  $
                   -9.*(f[l1+2:l2+2,m,n+3,*]-f[l1-2:l2-2,m,n+3,*])  $
                      +(f[l1+3:l2+3,m,n+3,*]-f[l1-3:l2-3,m,n+3,*])) $
                 -(45.*(f[l1+1:l2+1,m,n-3,*]-f[l1-1:l2-1,m,n-3,*])  $
                   -9.*(f[l1+2:l2+2,m,n-3,*]-f[l1-2:l2-2,m,n-3,*])  $
                      +(f[l1+3:l2+3,m,n-3,*]-f[l1-3:l2-3,m,n-3,*]))) )
      endfor & endfor
    endif
  endif else begin
    print, 'error: xderzder_6th_ghost not implemented for ', $
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
