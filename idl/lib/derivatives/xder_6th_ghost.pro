;;
;;  $Id: xder_6th_ghost.pro,v 1.13 2008-06-10 13:07:41 ajohan Exp $
;;
;;  First derivative d/dx
;;  - 6th-order
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function xder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
;
;  Default values.
;
  default, ghost, 0
;
;  Calculate nx, ny, and nz, based on the input array size.
;
  s=size(f) & d=make_array(size=s)
  nx=s[1] & ny=s[2] & nz=s[3]
; 26-jun-2007/dintrans: 2-D case only means (x,z) for the moment
  if (s[0] eq 2) then nz=s[2]
;
;  Check for degenerate case (no x-extension).
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (nx eq 1) then return,fltarr(nx,ny,nz)
;
;  Determine location of ghost zones, assume nghost=3 for now.
;
  l1=3 & l2=nx-4
;
  if (lequidist[0]) then begin
    dx2=1./(60.*(x[4]-x[3]))
  endif else begin
    dx2=dx_1[l1:l2]/60.
  endelse
;
  if (s[0] eq 3) then begin
    if (l2 gt l1) then begin
      if (lequidist[0] eq 0) then dx2=spread(dx2,[1,2],[s[2],s[3]])
      d[l1:l2,*,*]=dx2*( +45.*(f[l1+1:l2+1,*,*]-f[l1-1:l2-1,*,*]) $
                          -9.*(f[l1+2:l2+2,*,*]-f[l1-2:l2-2,*,*]) $
                             +(f[l1+3:l2+3,*,*]-f[l1-3:l2-3,*,*]) )
    endif else begin
      d[l1:l2,*,*]=0.
    endelse
;
  endif else if (s[0] eq 2) then begin
;
    if (l2 gt l1) then begin
      if (lequidist[0] eq 0) then dx2=spread(dx2,1,s[2])
      d[l1:l2,*]=dx2*( +45.*(f[l1+1:l2+1,*]-f[l1-1:l2-1,*]) $
                        -9.*(f[l1+2:l2+2,*]-f[l1-2:l2-2,*]) $
                           +(f[l1+3:l2+3,*]-f[l1-3:l2-3,*]) )
    endif else d[l1:l2,*]=0.
;
  endif else if (s[0] eq 4) then begin
;
    if (l2 gt l1) then begin
      if (lequidist[0] eq 0) then dx2=spread(dx2,[1,2,3],[s[2],s[3],s[4]])
      d[l1:l2,*,*,*]=dx2*( +45.*(f[l1+1:l2+1,*,*,*]-f[l1-1:l2-1,*,*,*]) $
                            -9.*(f[l1+2:l2+2,*,*,*]-f[l1-2:l2-2,*,*,*]) $
                               +(f[l1+3:l2+3,*,*,*]-f[l1-3:l2-3,*,*,*]) )
    endif else begin
      d[l1:l2,*,*,*]=0.
    endelse
;
  endif else begin
    print, 'error: xder_6th_ghost not implemented for ', $
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
