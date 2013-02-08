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
;       Written by:     Chao-Chin Yang, February 8, 2013.
;-

pro pc2vtk_allvar, datadir=datadir, variables=variables, bbtoo=bbtoo, ootoo=ootoo, trimall=trimall
  compile_opt idl2

; Find the number snapshots.
  if n_elements(datadir) eq 0 then datadir = pc_get_datadir()
  list = datadir + '/proc0/varN.list'
  nvarN = file_lines(list)

; Process each snapshot.
  openr, lun, list, /get_lun
  varfile = ''
  for i = 1, nvarN do begin
    readf, lun, varfile
    pc2vtk, datadir=datadir, varfile=varfile, variables=variables, bbtoo=bbtoo, ootoo=ootoo, trimall=trimall
  endfor
  close, lun
  free_lun, lun

end
