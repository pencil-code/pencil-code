;
;  $Id: zder_6th_ghost.pro,v 1.9 2005-10-24 02:12:21 dobler Exp $
;
;  6th-order first derivative in x direction for date with ghost cells and on
;  potentially non-equidistant grid.
;
;***********************************************************************
function zder,f
  COMPILE_OPT IDL2,HIDDEN
;
  common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
;
;  Check if we have a degenerate case (no x-extension)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (nz eq 1) then return,fltarr(nx,ny,nz)
  s=size(f) & d=make_array(size=s)
;
  n1=3 & n2=nz-4
;
  if (lequidist[2]) then begin
    dz2=1./(60.*(z[4]-z[3])) 
  endif else begin
    dz2=spread(dz_1[n1:n2]/60.,[0,1],[nx,ny])
  endelse
;
  if (s[0] eq 3) then begin
    if (n2 gt n1) then begin
      d[*,*,n1:n2]=dz2*( +45.*(f[*,*,n1+1:n2+1]-f[*,*,n1-1:n2-1]) $
                          -9.*(f[*,*,n1+2:n2+2]-f[*,*,n1-2:n2-2]) $
                             +(f[*,*,n1+3:n2+3]-f[*,*,n1-3:n2-3]) )
    endif else begin
      d[*,*,n1:n2]=0.
    endelse
;
  endif else if (s[0] eq 4) then begin
;
    if (n2 gt n1) then begin
      d[*,*,n1:n2,*]=dz2*( +45.*(f[*,*,n1+1:n2+1,*]-f[*,*,n1-1:n2-1,*]) $
                            -9.*(f[*,*,n1+2:n2+2,*]-f[*,*,n1-2:n2-2,*]) $
                               +(f[*,*,n1+3:n2+3,*]-f[*,*,n1-3:n2-3,*]) )
    endif else begin
      d[*,*,n1:n2,*]=0.
    endelse
;
  endif else begin
    print, 'error: zder_6th_ghost not implemented for ', $
           strtrim(s[0],2), '-D arrays'
  endelse
;
  return, d
;
end
