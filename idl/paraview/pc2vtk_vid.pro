pro pc2vtk_vid,variables=vars,file=file,tinit=tinit,tend=tend

default,vars,['rho','uu','bb']
default,file,'work'
default,varfile,'var.dat'

pc_read_grid,o=g,/trimxyz,/quiet
for it = tinit, tend do begin
; reading pc variables and setting dimensions
  var = vars
  print, 'reading ' + 'VAR'+strtrim(it,2)
  pc_read_var,varfile='VAR'+strtrim(it,2),obj=o,variables=var,/magic,/trimall
  print, 'varfile read'

  dim = size(o.rho)
  dimx = dim[1]
  dimy = dim[2]
  dimz = dim[3]
  dx = (max(g.x) - min(g.x))/(dimx)
  dy = (max(g.y) - min(g.y))/(dimy)
  dz = (max(g.z) - min(g.z))/(dimz)
  
  fl = file+str(it)+'.vtk'
  openw, iu, fl, /get_lun
  printf, iu, '# vtk DataFile Version 2.0'
  printf, iu, 'density + magnetic field'
  printf, iu, 'BINARY'
  printf, iu, 'DATASET STRUCTURED_POINTS'
  printf, iu, format="(a,3i9)", 'DIMENSIONS ', dimx,dimy,dimz
  printf, iu, format="(a,3f12.8)", 'ORIGIN ', g.x[0], g.y[0], g.z[0]
  printf, iu, format="(a,3f12.8)", 'SPACING ', dx, dy, dz
  printf, iu, format="(a,i9)", 'POINT_DATA ', dim[5]
  printf, iu, 'SCALARS rho float'
  printf, iu, 'LOOKUP_TABLE default'
  for k = 0, dimz-1 do $
     writeu, iu, swap_endian(o.rho[*,*,k])
  printf, iu, 'VECTORS vfield float'
  for k = 0, dimz-1 do $
    for j = 0, dimy-1 do $
     for i = 0, dimx-1 do $
        writeu, iu, swap_endian(o.uu[i,j,k,0]), swap_endian(o.uu[i,j,k,1]), swap_endian(o.uu[i,j,k,2])
  printf, iu, 'VECTORS bfield float'
  for k = 0, dimz-1 do $
    for j = 0, dimy-1 do $
      for i = 0, dimx-1 do $
        writeu, iu, swap_endian(o.bb[i,j,k,0]), swap_endian(o.bb[i,j,k,1]), swap_endian(o.bb[i,j,k,2])
  close, iu

end
end
