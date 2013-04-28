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
;  Check for degenerate case (no z-derivative)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (mz eq 1) then return, fltarr(mx,my,mz)
;
  l1=nghost & l2=mx-nghost-1
  m1=nghost & m2=my-nghost-1
  n1=nghost & n2=mz-nghost-1
;
  nx = mx - 2*nghost
  ny = my - 2*nghost
  nz = mz - 2*nghost
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
  if (s[0] eq 3) then begin
    if (not ldegenerated[2]) then begin
      if (not lequidist[2]) then begin
        dz2 =    spread(dz2,     [0,1],[nx,ny])
        dd  = d1*spread(dz_tilde,[0,1],[mx,my])
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
    if (not ldegenerated[2]) then begin
      if (not lequidist[2]) then begin
        dz2 =    spread(dz2,     [0,1,3],[nx,ny,s[4]])
        dd  = d1*spread(dz_tilde,[0,1,3],[mx,my,s[4]])
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
