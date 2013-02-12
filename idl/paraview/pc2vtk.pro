; $Id$
;+
; NAME:
;	PC2VTK
;
; PURPOSE:
;	This procedure converts one snapshot into legacy VTK format.
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
; MODIFICATION HISTORY:
;       Written by:    	Chao-Chin Yang, February 12, 2013.
;-

pro pc2vtk, datadir=datadir, varfile=varfile, variables=variables, bbtoo=bbtoo, ootoo=ootoo, trimall=trimall
  compile_opt idl2

; Read the data.
  if n_elements(datadir) eq 0 then datadir = pc_get_datadir()
  if n_elements(varfile) eq 0 then varfile = 'var.dat'
  print, 'Reading ', datadir, '/', varfile, '...'
  pc_read_param, obj=par, datadir=datadir, /quiet
  pc_read_var, obj=f, datadir=datadir, varfile=varfile, variables=variables, bbtoo=bbtoo, ootoo=ootoo, trimall=trimall, /quiet

; Convert and write the VTK file.
  if varfile eq 'var.dat' then vtkfile = 'var.vtk' else vtkfile = strtrim(varfile) + '.vtk'
  vtkfile = strtrim(datadir) + '/' + vtkfile
  grid = ~(par.lequidist[0] && par.lequidist[1] && par.lequidist[2])
  write_vtk, f, vtkfile, grid=grid

end

