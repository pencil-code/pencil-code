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
;	PC2VTK_ALLVAR [, DATADIR=string] [, VARIABLES=array] [/BBTOO]
;		[/OOTOO] [/TRIMALL]
;
; KEYWORDS:
;	DATADIR:	Same keyword as in PC_READ_VAR
;	VARIABLES:	Same keyword as in PC_READ_VAR
;	BBTOO:		Same keyword as in PC_READ_VAR
;	OOTOO:		Same keyword as in PC_READ_VAR
;	TRIMALL:	Same keyword as in PC_READ_VAR
;
; OUTPUTS:
;	This procedure saves each snapshot in legacy VTK format to file
;	VAR*.vtk.
;
; MODIFICATION HISTORY:
;       Written by:     Chao-Chin Yang, February 18, 2013.
;-

pro pc2vtk_allvar, datadir=datadir, variables=variables, bbtoo=bbtoo, ootoo=ootoo, trimall=trimall
  compile_opt idl2

; Read the list of snapshots.
  varlist = pc_read_varlist(datadir=datadir, nvar=nvar)

; Process each snapshot.
  for i = 0, nvar - 1 do begin
    if n_elements(variables) eq 0 then undefine, vars else vars = variables
    pc2vtk, datadir=datadir, varfile=varlist[i], variables=vars, bbtoo=bbtoo, ootoo=ootoo, trimall=trimall
  endfor

end
