; $Id$
;+
; NAME:
;	PC_READ_VARLIST
;
; PURPOSE:
;	This function returns a list of the file names of the snapshots.
;
; CATEGORY:
;	The Pencil Code - I/O
;
; CALLING SEQUENCE:
;	list = PC_READ_VARLIST([, DATADIR=string] [, NVAR=variable]
;		[, /PARTICLES])
;
; KEYWORDS:
;	DATADIR:	Set this keyword to a string containing the full
;		path of the data directory.  If omitted, './data' is
;		assumed.
;	NVAR:	Set this keyword to a variable that will contain the
;		total number of snapshots in the list.
;	PARTICLES:	Set this keyword to read the list of particle
;		data instead of fluid.
;
; MODIFICATION HISTORY:
;       Written by:     Chao-Chin Yang, February 18, 2013.
;-

function pc_read_varlist, datadir=datadir, nvar=nvar, particles=particles, allprocs=allprocs
  compile_opt idl2

  default, procdir, '/proc0/'
  default, list_file, 'varN.list'

; Find the list of snapshots.
  datadir = pc_get_datadir(datadir)
  if (size (allprocs, /type) ne 0) then begin
    if (keyword_set (allprocs)) then procdir = '/allprocs/'
  end else if (file_test (datadir+'/allprocs/'+list_file)) then begin
    procdir = '/allprocs/'
  end
  if (keyword_set (particles)) then list_file = 'pvarN.list'
  list_file = datadir + procdir + 'varN.list'

; Read the file name of each snapshot.
  nvar = file_lines (list_file)
  var_list = strarr (nvar)
  openr, lun, list_file, /get_lun
  for i = 0, nvar-1 do begin
    entry = ''
    readf, lun, entry
    var_list[i] = (stregex(entry,'([A-Z]+[0-9]+).*$',/SUBEXPR,/EXTRACT))[1]
  end
  close, lun
  free_lun, lun

  return, var_list

end

