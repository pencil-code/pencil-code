; Plot iso-surfaces of the square of a given vector field
; Useful to plot isosurfaces of B^2, and \Omega^2 from
; three dimensional simulations
PRO isoplot, omega, thold
  s = SIZE(omega)
  nx=s[1]
  ny=s[2]
  nz=s[3]
  PRINT, 'the size of array:',nx,ny,nz
  psi=findgen(nx,ny,nz)
  msval = 0.
  for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
      for iz=0,nz-1 do begin
        psi[ix,iy,iz] = omega[ix,iy,iz,0]*omega[ix,iy,iz,0]+omega[ix,iy,iz,1]*omega[ix,iy,iz,1]+ $ 
                                   omega[ix,iy,iz,2]*omega[ix,iy,iz,2]
        msval = msval + psi[ix,iy,iz]
      endfor
    endfor
  endfor
  msval = msval/(nx*ny*nz)
  PRINT, 'mean-square value=',msval
  SHADE_VOLUME, psi, thold*msval, v, p
  SCALE3, XRANGE=[0,S[1]], YRANGE=[0,S[2]], ZRANGE=[0,S[3]], AX=0, AZ=45
  TV, POLYSHADE(v, p, /T3D)
; Define number of views:  
nframes = 20   
FOR i = 0, nframes - 1 DO BEGIN & $  
   ; Translate the center of the (0, 1) unit cube  
   ; to (0,0) and rotate about the x-axis:  
   T3D, TR=[-.5, -.5, -.5], ROT=[0, 360./NFRAMES, 0] & $  
   ; Translate the center back to (0.5, 0.5, 0.5):  
   T3D, TR = [.5, .5, .5] & $  
   ; Show the surface:  
   TV, POLYSHADE(v, p, /T3D) & $  
ENDFOR  
END
