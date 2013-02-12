; $Id$
;+
; NAME:
;	WRITE_VTK
;
; PURPOSE:
;	This procedure converts the Pencil Code data into legacy VTK
;	format.
;
; CATEGORY:
;	The Pencil Code - I/O
;
; CALLING SEQUENCE:
;	WRITE_VTK, Obj [, File] [, X=array] [, Y=array] [, Z=array]
;		[, SPACING=array] [, /GRID]
;
; ARGUMENTS:
;	Obj:	Data structure produced by any of the reading procedure
;		from the Pencil Code.
;	File:	File name including the full path (if any) of the output
;		VTK file.  If this argument is omitted, the file name
;		'work.vtk' is used.
;
; KEYWORDS:
;	X:	A 1D array containing the x-coordinates of the grid.
;	Y:	A 1D array containing the y-coordinates of the grid.
;	Z:	A 1D array containing the z-coordinates of the grid.
;	SPACING:	Three-element array containing the spacing in
;		each dimension.
;	GRID:	Set this keyword for writing non-uniform grid.
;
; MODIFICATION HISTORY:
;       Written by:    	Chao-Chin Yang, February 12, 2013.
;-

pro write_vtk, obj, file, x=x, y=y, z=z, spacing=spacing, grid=grid
  compile_opt idl2

; Determine the coordinates and the precision.
  if n_elements(file) eq 0 then file = 'work.vtk'
  if n_elements(x) eq 0 then x = obj.x
  if n_elements(y) eq 0 then y = obj.y
  if n_elements(z) eq 0 then z = obj.z
  if n_elements(spacing) eq 0 then spacing = [obj.dx, obj.dy, obj.dz]
  if keyword_set(grid) then rectilinear_grid = 1 else rectilinear_grid = 0
  data_type = strlowcase(size(x, /tname))

; Find the dimensions of the data.
  nx = n_elements(x)
  ny = n_elements(y)
  nz = n_elements(z)
  ntot = nx * ny * nz
  dimensions = [nx, ny, nz]
  origin = [x[0], y[0], z[0]]
  ndim = n_elements(where(dimensions gt 1))
  
; Open the VTK file for write.
  print, 'Writing ', strtrim(file), '...'
  openw, lun, file, /get_lun, /swap_if_little_endian

; Write the header information.
  printf, lun, '# vtk DataFile Version 2.0'
  printf, lun, 'Pencil Code Data'
  printf, lun, 'BINARY'
  if rectilinear_grid then begin
    printf, lun, 'DATASET RECTILINEAR_GRID'
    printf, lun, 'DIMENSIONS ', dimensions
    printf, lun, 'X_COORDINATES ', nx, ' ', data_type
    writeu, lun, x
    printf, lun, 'Y_COORDINATES ', ny, ' ', data_type
    writeu, lun, y
    printf, lun, 'Z_COORDINATES ', nz, ' ', data_type
    writeu, lun, z
  endif else begin
    printf, lun, 'DATASET STRUCTURED_POINTS'
    printf, lun, 'DIMENSIONS ', dimensions
    printf, lun, 'ORIGIN ', origin
    printf, lun, 'SPACING ', spacing
  endelse
  printf, lun, 'POINT_DATA ', ntot

; Write out each data field.
  tags = tag_names(obj)
  for i = 0, n_tags(obj) - 1 do begin
    ; Skip non-data tags.
    if tags[i] eq 'X' or tags[i] eq 'Y' or tags[i] eq 'Z' or $
       tags[i] eq 'DX' or tags[i] eq 'DY' or tags[i] eq 'DZ' or $
       tags[i] eq 'T' or tags[i] eq 'DELTAY' then continue
    if (size(obj.(i)))[0] eq ndim then begin
      ; Write scalar field.
      printf, lun, 'SCALARS ', strlowcase(tags[i]), ' ', data_type
      printf, lun, 'LOOKUP_TABLE default'
      for j = 0, ntot - 1 do writeu, lun, (obj.(i))[j]
    endif else begin
      ; Write vector field.
      printf, lun, 'VECTORS ', strlowcase(tags[i]), ' ', data_type
      data = reform(obj.(i), [ntot,3])
      for j = 0, ntot - 1 do writeu, lun, data[j,0], data[j,1], data[j,2]
    endelse
  endfor
  close, lun
  free_lun, lun

end

