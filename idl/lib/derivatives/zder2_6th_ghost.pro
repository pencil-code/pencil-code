;
;  $Id: zder2_6th_ghost.pro,v 1.7 2005-10-24 08:19:12 dobler Exp $
;
;  Second derivative d^2/dz^2
;  - 6th-order (7-point stencil)
;  - with ghost cells
;  - on potentially non-equidistant grid
;
;***********************************************************************
function zder2,f
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0 
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
;
;  Check for degenerate case (no x-extension)
;
  if (n_elements(lequidist) ne 3) then lequidist=[-1,-1,-1]
  if (nz eq 1) then return,fltarr(nx,ny,nz)
  s=size(f) & d=make_array(size=s)
;
  n1=3 & n2=nz-4
;
  if (lequidist[2]) then begin
    dz2=1./(180.*(z[4]-z[3])^2)
  endif else begin
    dz2=dz_1[n1:n2]^2/180.
;
;  nonuniform mesh correction
;  d2f/dz2 = zeta'^2*f" + zeta"*f'   see also the manual
;
    d1=zder(f)
  endelse
;
  if (s[0] eq 3) then begin
    if (n2 gt n1) then begin
      if (lequidist[2] eq 0) then begin
        dz2 =    spread(dz2,     [0,0],[s[2],s[1]])
        dd  = d1*spread(dz_tilde,[0,0],[s[2],s[1]])
        ; will also work on slices like zder2(ss[10,20,*])
      endif
      d[*,*,n1:n2]=dz2*(-490.*f[*,*,n1:n2]$
                       +270.*(f[*,*,n1-1:n2-1]+f[*,*,n1+1:n2+1])$
                        -27.*(f[*,*,n1-2:n2-2]+f[*,*,n1+2:n2+2])$
                         +2.*(f[*,*,n1-3:n2-3]+f[*,*,n1+3:n2+3])$
                       )
    endif else begin
      d[*,*,n1:n2]=0.
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
      d[*,*,n1:n2,*]=dz2*(-490.*f[*,*,n1:n2,*]$
                         +270.*(f[*,*,n1-1:n2-1,*]+f[*,*,n1+1:n2+1,*])$
                          -27.*(f[*,*,n1-2:n2-2,*]+f[*,*,n1+2:n2+2,*])$
                           +2.*(f[*,*,n1-3:n2-3,*]+f[*,*,n1+3:n2+3,*])$
                       )
    endif else begin
      d[*,*,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: zder2_6th_ghost not implemented for ', $
           strtrim(s[0],2), '-D arrays'
  endelse
;
; apply correction only for nonuniform mesh
;
  if not lequidist[2] then d=d+dd
;
  return,d
;
end
