;;
;;  $Id$
;;
;;  Sixth derivative d^6/dy^6
;;  - 6th-order (7-point stencil)
;;  - with ghost cells
;;
function yder6,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z
  common cdat_grid,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist,lperi,ldegenerated
;
;  Default values.
;
  default, ghost, 0
;
;  Calculate nx, ny, and nz, based on the input array size.
;
  s=size(f) & d=make_array(size=s)
  nx=s[1] & ny=s[2] & nz=s[3]
;
;  Check for degenerate case (no y-derivative)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (ny eq 1) then return,fltarr(nx,ny,nz)
;
;  Determine location of ghost zones, assume nghost=3 for now.
;
   l1=3 & l2=nx-4 & m1=3 & m2=ny-4 & n1=3 & n2=nz-4
;
  if (lequidist[1]) then begin
    dy6=1./(y[4]-y[3])^6
  endif else begin
;
;  Nonuniform mesh not implemented.
;
    print, 'Nonuniform mesh not implemented for yder6_6th_ghost.pro'
    stop
  endelse
;
  if (s[0] eq 3) then begin
    if (not ldegenerated[1]) then begin
      d[l1:l2,m1:m2,n1:n2]=dy6* $
          ( -20.*f[l1:l2,m1:m2,n1:n2] $
            +15.*(f[l1:l2,m1-1:m2-1,n1:n2]+f[l1:l2,m1+1:m2+1,n1:n2]) $
             -6.*(f[l1:l2,m1-2:m2-2,n1:n2]+f[l1:l2,m1+2:m2+2,n1:n2]) $
             +1.*(f[l1:l2,m1-3:m2-3,n1:n2]+f[l1:l2,m1+3:m2+3,n1:n2]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (not ldegenerated[1]) then begin
      d[l1:l2,m1:m2,n1:n2,*]=dy6* $
          ( -20.*f[l1:l2,m1:m2,n1:n2,*] $
            +15.*(f[l1:l2,m1-1:m2-1,n1:n2,*]+f[l1:l2,m1+1:m2+1,n1:n2,*]) $
             -6.*(f[l1:l2,m1-2:m2-2,n1:n2,*]+f[l1:l2,m1+2:m2+2,n1:n2,*]) $
             +1.*(f[l1:l2,m1-3:m2-3,n1:n2,*]+f[l1:l2,m1+3:m2+3,n1:n2,*]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: yder6_6th_ghost not implemented for ', $
           strtrim(s[0],2), '-D arrays'
    stop
  endelse
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
