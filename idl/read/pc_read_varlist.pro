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

function pc_read_varlist, datadir=datadir, nvar=nvar, particles=particles
  compile_opt idl2

; Find the number snapshots.
  if n_elements(datadir) eq 0 then datadir = pc_get_datadir()
  if keyword_set(particles) then list = datadir + '/proc0/pvarN.list' else list = datadir + '/proc0/varN.list'
  nvar = file_lines(list)

; Read the file name of each snapshot.
  openr, lun, list, /get_lun
  varfile = ''
  readf, lun, varfile
  varlist = varfile
  for i = 2, nvar do begin
    readf, lun, varfile
    varlist = [varlist, varfile]
  endfor
  close, lun
  free_lun, lun

  return, varlist

end

