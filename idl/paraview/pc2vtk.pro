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
;	PC2VTK [, DATADIR=string] [, VARFILE=string] [, <extra parameters>]
;
; KEYWORDS:
;	DATADIR:	Same keyword as in PC_READ_VAR
;	VARFILE:	Same keyword as in PC_READ_VAR
;       All other (keyword) parameters are assumed to be dedicated to pc_read_var, e.g.
;	VARIABLES, BBTOO, OOTOO, TRIMALL, etc.
;       Parameters not known to pc_read_var are silently ignored.
;
; MODIFICATION HISTORY:
;       Written by:    	Chao-Chin Yang, February 12, 2013.
;       Modified by:    M. Rheinhardt, October 3rd, 2026
;-
pro pc2vtk, datadir=datadir, varfile=varfile, _extra=extra
  compile_opt idl2

  ; Read the data.
  datadir = pc_get_datadir(datadir)

  if is_defined(extra) then $
    tag_ivar=tag_exists(extra,'ivar') $
  else $
    tag_ivar=0

  if is_defined(varfile) then begin
    if tag_ivar then extra=remove_tag(extra,'ivar')
  endif else $
    if tag_ivar then $
      print, 'Reading ', datadir, '/VAR'+strtrim(string(ivar),2) $
    else begin
      varfile = 'var.dat'
      print, 'Reading ', datadir, '/', varfile, '...'
    endelse

  pc_read_param, obj=par, datadir=datadir, /quiet
  pc_read_var, obj=f, datadir=datadir, varfile=varfile, /quiet, _extra=extra
  if varfile eq '' then return

; Convert and write the VTK file.
  if varfile eq 'var.dat' then vtkfile = 'var.vtk' else vtkfile = strtrim(varfile) + '.vtk'
  vtkfile = strtrim(datadir) + '/' + vtkfile
  grid = ~(par.lequidist[0] && par.lequidist[1] && par.lequidist[2])
  write_vtk, f, vtkfile, grid=grid

end

