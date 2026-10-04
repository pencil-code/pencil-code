; $Id$
;+
; NAME:
;	PC2VTK_ALLVAR
;
; PURPOSE:
;	This procedure converts the Pencil Code data of each snapshot
;	into legacy VTK format.
;
; CATEGORY:
;	The Pencil Code - I/O
;
; CALLING SEQUENCE:
;	PC2VTK_ALLVAR [, DATADIR=string] [, VARIABLES=array] [, /DOWN] [, <extra parameters>]
;
; KEYWORDS:
;	DATADIR:	Same keyword as in PC_READ_VAR
;	VARIABLES:	Same keyword as in PC_READ_VAR
;       DOWN:           Read downsampled snapshots.
;       All other (keyword) parameters are assumed to be dedicated to pc_read_var, e.g.
;       BBTOO, OOTOO, TRIMALL, etc.
;       Parameters not known to pc_read_var are silently ignored.
;
; OUTPUTS:
;	This procedure saves each snapshot in legacy VTK format to file VAR*.vtk.
;
; MODIFICATION HISTORY:
;       Written by:     Chao-Chin Yang, February 18, 2013.
;       Modified by:    M. Rheinhardt, October 3rd, 2026
;-
pro pc2vtk_allvar, datadir=datadir, variables=variables, down=down, _extra=extra
  compile_opt idl2

  if is_defined(extra) then begin
    if tag_exists(extra,'varfile') then extra=remove_tag(extra,'varfile')
    if tag_exists(extra,'ivar') then extra=remove_tag(extra,'ivar')
  endif

; Read the list of snapshots.
  varlist = pc_read_varlist(datadir=datadir, nvar=nvar, down=down)
; Process each snapshot.
  for i = 0, nvar - 1 do begin
    if n_elements(variables) eq 0 then undefine, vars else vars = variables
    varfile=varlist[i]
    pc2vtk, datadir=datadir, varfile=varfile, variables=vars, _extra=extra
    if varfile eq '' then return
  endfor

end
