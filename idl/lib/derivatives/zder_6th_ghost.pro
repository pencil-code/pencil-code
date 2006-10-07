;
;  $Id: zder_6th_ghost.pro,v 1.11 2006-10-07 09:57:56 brandenb Exp $
;
;  First derivative d/dz
;  - 6th-order (7-point stencil)
;  - with ghost cells
;  - on potentially non-equidistant grid
;
;***********************************************************************
function zder,f
  COMPILE_OPT IDL2,HIDDEN
;
  ;common cdat,x,y,z,nx,ny,nz,nw,ntmax,date0,time0
  ;AB: chose to read in only x, y, and z, not nx, ny, and nz.
  ;AB: Thus, we can redefine them freely.
  ;AB: For non-uniform meshes dx_1, dy_1, and dz_1 would not be ok.
  common cdat,x,y,z
  common cdat_nonequidist,dx_1,dy_1,dz_1,dx_tilde,dy_tilde,dz_tilde,lequidist
;
;  calculate nx, ny, and nz, based on the input array size
;
  s=size(f) & d=make_array(size=s)
  nx=s[1] & ny=s[2] & nz=s[3]
;
;  Check for degenerate case (no z-extension)
;
  if (n_elements(lequidist) ne 3) then lequidist=[1,1,1]
  if (nz eq 1) then return,fltarr(nx,ny,nz)
;
;  determine location of ghost zones, assume nghost=3 for now.
;
  n1=3 & n2=nz-4
;
  if (lequidist[2]) then begin
    dz2=1./(60.*(z[4]-z[3])) 
  endif else begin
    dz2=dz_1[n1:n2]/60.
  endelse
;
  if (s[0] eq 3) then begin
    if (n2 gt n1) then begin
      if (lequidist[2] eq 0) then dz2=spread(dz2,[0,0],[s[2],s[1]])
      ; will also work on slices like zder(ss[10,20,*])
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
      if (lequidist[2] eq 0) then dz2=spread(dz2,[0,0,3],[s[2],s[1],s[4]])
      ; will also work on slices like zder(uu[10,20,*,*])
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
