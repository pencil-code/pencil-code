;;
;;  $Id$
;;
;;  Second derivative d^2/dy^2
;;  - 6th-order (7-point stencil)
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function yder2,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
;
;  Check for degenerate case (no x-extension).
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (ny eq 1) then return, fltarr(nx,ny,nz)
;
;  Determine location of ghost zones, assume nghost=3 for now.
;
  l1=3 & l2=nx-4 & m1=3 & m2=ny-4 & n1=3 & n2=nz-4
;
  if (lequidist[1]) then begin
    dy2=1./(180.*(y[4]-y[3])^2)
  endif else begin
    dy2=dy_1[m1:m2]^2/180.
;
;  Nonuniform mesh correction.
;  d2f/dy2  = f"*psi'^2 + psi"f', see also the manual.
;
    d1=yder(f)
  endelse
;
  if (s[0] eq 2) then begin
    if (m2 gt m1) then begin
      if (lequidist[1] eq 0) then begin
        dy2 =    spread(dy2,     [0,2],[s[1],s[3]])
        dd  = d1*spread(dy_tilde,[0,2],[s[1],s[3]])
        ; will also work on slices like yder2(ss[10,*,*])
      endif
      d[l1:l2,m1:m2]=dy2* $
          (-490.*f[l1:l2,m1:m2] $
           +270.*(f[l1:l2,m1-1:m2-1]+f[l1:l2,m1+1:m2+1]) $
            -27.*(f[l1:l2,m1-2:m2-2]+f[l1:l2,m1+2:m2+2]) $
             +2.*(f[l1:l2,m1-3:m2-3]+f[l1:l2,m1+3:m2+3]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 3) then begin
    if (m2 gt m1) then begin
      if (lequidist[1] eq 0) then begin
        dy2 =    spread(dy2,     [0,2],[s[1],s[3]])
        dd  = d1*spread(dy_tilde,[0,2],[s[1],s[3]])
        ; will also work on slices like yder2(ss[10,*,*])
      endif
      d[l1:l2,m1:m2,n1:n2]=dy2* $
          (-490.*f[l1:l2,m1:m2,n1:n2] $
           +270.*(f[l1:l2,m1-1:m2-1,n1:n2]+f[l1:l2,m1+1:m2+1,n1:n2]) $
            -27.*(f[l1:l2,m1-2:m2-2,n1:n2]+f[l1:l2,m1+2:m2+2,n1:n2]) $
             +2.*(f[l1:l2,m1-3:m2-3,n1:n2]+f[l1:l2,m1+3:m2+3,n1:n2]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (m2 gt m1) then begin
      if (lequidist[1] eq 0) then begin
        dy2 =    spread(dy2,     [0,2,3],[s[1],s[3],s[4]])
        dd  = d1*spread(dy_tilde,[0,2,3],[s[1],s[3],s[4]])
        ; will also work on slices like yder2(uu[10,*,*,*])
      endif
      d[l1:l2,m1:m2,n1:n2,*]=dy2* $
          (-490.*f[l1:l2,m1:m2,n1:n2,*] $
           +270.*(f[l1:l2,m1-1:m2-1,n1:n2,*]+f[l1:l2,m1+1:m2+1,n1:n2,*]) $
            -27.*(f[l1:l2,m1-2:m2-2,n1:n2,*]+f[l1:l2,m1+2:m2+2,n1:n2,*]) $
             +2.*(f[l1:l2,m1-3:m2-3,n1:n2,*]+f[l1:l2,m1+3:m2+3,n1:n2,*]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: yder2_6th_ghost not implemented for ', $
           strtrim(s[0],2), '-D arrays'
  endelse
;
;  Apply correction only for nonuniform mesh.
;
  if (not lequidist[1]) then d=d+dd
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
