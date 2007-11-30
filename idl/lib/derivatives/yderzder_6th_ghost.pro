;
;  $Id: yderzder_6th_ghost.pro,v 1.2 2007-11-30 14:29:46 ajohan Exp $
;
;  Second derivative d2/dydz
;  - 6th-order
;  - with ghost cells
;
;***********************************************************************
function yderzder,f
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
    print, 'yderzder is not implemented for non equidistant grid'
    stop
  endif
;
;  Determine location of ghost zones - assume nghost=3 for now.
;
  l1=3 & l2=mx-4
  m1=3 & m2=my-4
  n1=3 & n2=mz-4
;
  dy_1=1/(y[4]-y[3])
  dz_1=1/(z[4]-z[3])
;
;  Calculate d2/dydz.
;
  if (s[0] eq 3) then begin
    d=fltarr(mx,my,mz)
    for n=n1,n2 do begin & for m=m1,m2 do begin
      fac=(1/60.0^2)*dy_1*dz_1
      d[l1:l2,m,n,*]=fac*( $
          45.*( (45.*(f[l1:l2,m+1,n+1]-f[l1:l2,m-1,n+1])  $
                 -9.*(f[l1:l2,m+2,n+1]-f[l1:l2,m-2,n+1])  $
                    +(f[l1:l2,m+3,n+1]-f[l1:l2,m-3,n+1])) $
               -(45.*(f[l1:l2,m+1,n-1]-f[l1:l2,m-1,n-1])  $
                 -9.*(f[l1:l2,m+2,n-1]-f[l1:l2,m-2,n-1])  $
                    +(f[l1:l2,m+3,n-1]-f[l1:l2,m-3,n-1])))$
          -9.*( (45.*(f[l1:l2,m+1,n+2]-f[l1:l2,m-1,n+2])  $
                 -9.*(f[l1:l2,m+2,n+2]-f[l1:l2,m-2,n+2])  $
                    +(f[l1:l2,m+3,n+2]-f[l1:l2,m-3,n+2])) $
               -(45.*(f[l1:l2,m+1,n-2]-f[l1:l2,m-1,n-2])  $
                 -9.*(f[l1:l2,m+2,n-2]-f[l1:l2,m-2,n-2])  $
                    +(f[l1:l2,m+3,n-2]-f[l1:l2,m-3,n-2])))$
             +( (45.*(f[l1:l2,m+1,n+3]-f[l1:l2,m-1,n+3])  $
                 -9.*(f[l1:l2,m+2,n+3]-f[l1:l2,m-2,n+3])  $
                    +(f[l1:l2,m+3,n+3]-f[l1:l2,m-3,n+3])) $
               -(45.*(f[l1:l2,m+1,n-3]-f[l1:l2,m-1,n-3])  $
                 -9.*(f[l1:l2,m+2,n-3]-f[l1:l2,m-2,n-3])  $
                    +(f[l1:l2,m+3,n-3]-f[l1:l2,m-3,n-3]))) )
    endfor & endfor
  endif else if (s[0] eq 4) then begin
    d=fltarr(mx,my,mz,3)
    for n=n1,n2 do begin & for m=m1,m2 do begin
      fac=(1/60.0^2)*dy_1*dz_1
      d[l1:l2,m,n,*]=fac*( $
          45.*( (45.*(f[l1:l2,m+1,n+1,*]-f[l1:l2,m-1,n+1,*])  $
                 -9.*(f[l1:l2,m+2,n+1,*]-f[l1:l2,m-2,n+1,*])  $
                    +(f[l1:l2,m+3,n+1,*]-f[l1:l2,m-3,n+1,*])) $
               -(45.*(f[l1:l2,m+1,n-1,*]-f[l1:l2,m-1,n-1,*])  $
                 -9.*(f[l1:l2,m+2,n-1,*]-f[l1:l2,m-2,n-1,*])  $
                    +(f[l1:l2,m+3,n-1,*]-f[l1:l2,m-3,n-1,*])))$
          -9.*( (45.*(f[l1:l2,m+1,n+2,*]-f[l1:l2,m-1,n+2,*])  $
                 -9.*(f[l1:l2,m+2,n+2,*]-f[l1:l2,m-2,n+2,*])  $
                    +(f[l1:l2,m+3,n+2,*]-f[l1:l2,m-3,n+2,*])) $
               -(45.*(f[l1:l2,m+1,n-2,*]-f[l1:l2,m-1,n-2,*])  $
                 -9.*(f[l1:l2,m+2,n-2,*]-f[l1:l2,m-2,n-2,*])  $
                    +(f[l1:l2,m+3,n-2,*]-f[l1:l2,m-3,n-2,*])))$
             +( (45.*(f[l1:l2,m+1,n+3,*]-f[l1:l2,m-1,n+3,*])  $
                 -9.*(f[l1:l2,m+2,n+3,*]-f[l1:l2,m-2,n+3,*])  $
                    +(f[l1:l2,m+3,n+3,*]-f[l1:l2,m-3,n+3,*])) $
               -(45.*(f[l1:l2,m+1,n-3,*]-f[l1:l2,m-1,n-3,*])  $
                 -9.*(f[l1:l2,m+2,n-3,*]-f[l1:l2,m-2,n-3,*])  $
                    +(f[l1:l2,m+3,n-3,*]-f[l1:l2,m-3,n-3,*]))) )
    endfor & endfor
  endif else begin
    print, 'error: yderzder_6th_ghost not implemented for ', $
        strtrim(s[0],2), '-D arrays'
  endelse
;
  return, d
;
end
