;;
;;  $Id$
;;
;;  Second derivative d^2/dz^2
;;  - 6th-order (7-point stencil)
;;  - with ghost cells
;;  - on potentially non-equidistant grid
;;
function zder2,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
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
  if (n_elements(lequidist) ne 3) then lequidist=[-1,-1,-1]
  if (nz eq 1) then return, fltarr(nx,ny,nz)
;
;  Determine location of ghost zones, assume nghost=3 for now.
;
  l1=3 & l2=nx-4 & m1=3 & m2=ny-4 & n1=3 & n2=nz-4
;
  if (lequidist[2]) then begin
    dz2=1./(180.*(z[4]-z[3])^2)
  endif else begin
    dz2=dz_1[n1:n2]^2/180.
;
;  Nonuniform mesh correction.
;  d2f/dz2 = zeta'^2*f" + zeta"*f', see also the manual.
;
    d1=zder(f)
  endelse
;
  if (s[0] eq 2) then begin
    d[l1:l2,m1:m2]=0.
;
  endif else if (s[0] eq 3) then begin
    if (n2 gt n1) then begin
      if (lequidist[2] eq 0) then begin
        dz2 =    spread(dz2,     [0,0],[s[2],s[1]])
        dd  = d1*spread(dz_tilde,[0,0],[s[2],s[1]])
        ; will also work on slices like zder2(ss[10,20,*])
      endif
      d[l1:l2,m1:m2,n1:n2]=dz2* $
          (-490.*f[l1:l2,m1:m2,n1:n2] $
           +270.*(f[l1:l2,m1:m2,n1-1:n2-1]+f[l1:l2,m1:m2,n1+1:n2+1]) $
            -27.*(f[l1:l2,m1:m2,n1-2:n2-2]+f[l1:l2,m1:m2,n1+2:n2+2]) $
             +2.*(f[l1:l2,m1:m2,n1-3:n2-3]+f[l1:l2,m1:m2,n1+3:n2+3]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (n2 gt n1) then begin
      if (lequidist[2] eq 0) then begin
        dz2 =    spread(dz2,     [0,0,3],[s[2],s[1],s[4]])
        dd  = d1*spread(dz_tilde,[0,0,3],[s[2],s[1],s[4]])
        ; will also work on slices like zder2(uu[10,20,*,*])
      endif
      d[l1:l2,m1:m2,n1:n2,*]=dz2* $
          (-490.*f[l1:l2,m1:m2,n1:n2,*] $
           +270.*(f[l1:l2,m1:m2,n1-1:n2-1,*]+f[l1:l2,m1:m2,n1+1:n2+1,*]) $
            -27.*(f[l1:l2,m1:m2,n1-2:n2-2,*]+f[l1:l2,m1:m2,n1+2:n2+2,*]) $
             +2.*(f[l1:l2,m1:m2,n1-3:n2-3,*]+f[l1:l2,m1:m2,n1+3:n2+3,*]) )
    endif else begin
      d[l1:l2,m1:m2,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: zder2_6th_ghost not implemented for ', $
           strtrim(s[0],2), '-D arrays'
  endelse
;
;  Apply correction only for nonuniform mesh.
;
  if (not lequidist[2]) then d=d+dd
;
;  Set ghost zones.
;
  if (ghost) then d=pc_setghost(d,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
  return, d
;
end
