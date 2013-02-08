; $Id$
;+
; NAME:
;	PC2VTK
;
; PURPOSE:
;	This procedure converts the Pencil Code data into legacy VTK
;	format.
;
; CATEGORY:
;	The Pencil Code - I/O
;
; CALLING SEQUENCE:
;	PC2VTK [, DATADIR=string] [, VARFILE=string] [, VARIABLES=array]
;		[/BBTOO] [/OOTOO] [/TRIMALL]
;
; KEYWORDS:
;	DATADIR:	Same keyword as in PC_READ_VAR
;	VARFILE:	Same keyword as in PC_READ_VAR
;	VARIABLES:	Same keyword as in PC_READ_VAR
;	BBTOO:		Same keyword as in PC_READ_VAR
;	OOTOO:		Same keyword as in PC_READ_VAR
;	TRIMALL:	Same keyword as in PC_READ_VAR
;
; OUTPUTS:
;	This procedure saves the data in legacy VTK format to file
;	VARFILE.vtk.
;
; MODIFICATION HISTORY:
;       Written by:     tavo.buk, April 15, 2011.
;			Chao-Chin Yang, February 8, 2013.
;-

pro pc2vtk, datadir=datadir, varfile=varfile, variables=variables, bbtoo=bbtoo, ootoo=ootoo, trimall=trimall
  compile_opt idl2

; Read the data.
  if n_elements(datadir) eq 0 then datadir = pc_get_datadir()
  if n_elements(varfile) eq 0 then varfile = 'var.dat'
  print, 'Reading ', datadir, '/', varfile, '...'
  pc_read_dim, obj=dim, /quiet
  pc_read_var, obj=f, datadir=datadir, dim=dim, varfile=varfile, trimall=trimall, /quiet
  if dim.precision eq 'D' then data_type = 'double' else data_type = 'float'

; Set the name of the output VTK file.
  if varfile eq 'var.dat' then vtkfile = 'var.vtk' else vtkfile = strtrim(varfile) + '.vtk'
  vtkfile = strtrim(datadir) + '/' + vtkfile

; Determine the dimensions of the data from the scalar field.
  if trimall then begin
    nx = dim.nxgrid
    ny = dim.nygrid
    nz = dim.nzgrid
  endif else begin
    nx = dim.mxgrid
    ny = dim.mygrid
    nz = dim.mzgrid
  endelse
  dimensions = [nx, ny, nz]
  origin = [f.x[0], f.y[0], f.z[0]]
  spacing = [f.dx[0], f.dy[0], f.dz[0]]
  dim_ok = where(dimensions gt 1)
  ntot = product(dimensions, /integer)

; Open the output file and write the header information.
  print, 'Writing ', strtrim(vtkfile), '...'
  openw, lun, vtkfile, /get_lun, /swap_if_little_endian
  printf, lun, '# vtk DataFile Version 2.0'
  printf, lun, 'Pencil Code Data'
  printf, lun, 'BINARY'
  printf, lun, 'DATASET STRUCTURED_POINTS'
  printf, lun, 'DIMENSIONS ', dimensions[dim_ok]
  printf, lun, 'ORIGIN ', origin[dim_ok]
  printf, lun, 'SPACING ', spacing[dim_ok]
  printf, lun, 'POINT_DATA ', ntot

; Write out each data field.
  ndim = n_elements(dim_ok)
  tags = tag_names(f)
  for i = 0, n_tags(f) - 1 do begin
    if tags[i] eq 'X' or tags[i] eq 'Y' or tags[i] eq 'Z' or $
       tags[i] eq 'DX' or tags[i] eq 'DY' or tags[i] eq 'DZ' or $
       tags[i] eq 'T' or tags[i] eq 'DELTAY' then continue
    if (size(f.(i)))[0] eq ndim then begin
      printf, lun, 'SCALARS ', strlowcase(tags[i]), ' ', data_type
      printf, lun, 'LOOKUP_TABLE default'
      for j = 0, ntot - 1 do writeu, lun, (f.(i))[j]
    endif else begin
      printf, lun, 'VECTORS ', strlowcase(tags[i]), ' ', data_type
      data = reform(f.(i), [ntot,3])
      for j = 0, ntot - 1 do writeu, lun, data[j,0], data[j,1], data[j,2]
    endelse
  endfor
  close, lun
  free_lun, lun

end

