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
;	list = PC_READ_VARLIST([, DATADIR=string] [, NVAR=variable] [, TIME=variable] [, /ALLPROCS]
;		[, /PARTICLES] [, /POINTMASSSES] [, /DOWN])
;
; KEYWORDS:
;	DATADIR:	Set this keyword to a string containing the full
;		path of the data directory.  If omitted, './data' is
;		assumed.
;	NVAR:	Set this keyword to a variable that will contain the
;		total number of snapshots in the list.
;	TIME:   Set this keyword to a variable that will contain the
;		times when the snapshots are written
;       ALLPROCS:  Look for varlist in data/allprocs.
;	PARTICLES: Set this keyword to read the list of particle
;		   data instead of fluid.
;       POINTMASSSES:   Read list of pointmasses snapshots.
;       DOWN:   Read list of downsampled snapshots.
;
; MODIFICATION HISTORY:
;       Written by:     Chao-Chin Yang, February 18, 2013.
;	Modified by:	Johannes Tschernitz, March 22, 2024. Added keyword time
;	Modified by:	Matthias Rheinhardt, October, 3rd, 2026. Added keywords down, pointmasses
;-
function pc_read_varlist, datadir=datadir, nvar=nvar, particles=particles, allprocs=allprocs, time=time, down=down, pointmasses=pointmasses
  compile_opt idl2

  default, procdir, '/proc0/'
  if (keyword_set (down)) then $
    list_file = 'varN_down.list' $
  else $
    list_file = 'varN.list'
  if (keyword_set (particles)) then list_file = 'pvarN.list'
  if (keyword_set (pointmasses)) then list_file = 'qvarN.list'

; Find the list of snapshots.
  datadir = pc_get_datadir(datadir)
  if (size (allprocs, /type) ne 0) then begin
    if (keyword_set (allprocs)) then procdir = '/allprocs/'
  end else if (file_test (datadir+'/allprocs/'+list_file)) then begin
    procdir = '/allprocs/'
  end
  list_file = datadir + procdir + list_file

; Read the file name of each snapshot.
  nvar = file_lines (list_file)
  var_list = strarr (nvar)
  if arg_present(time) then time = dblarr(nvar)
  openr, lun, list_file, /get_lun
  for i = 0, nvar-1 do begin
    entry = ''
    readf, lun, entry
    if arg_present(time) then time[i] = double((stregex(entry,'([0-9]\.[0-9]+).*$',/SUBEXPR,/EXTRACT))[0])
    var_list[i] = (stregex(entry,'([A-Z]+d*[0-9]+).*$',/SUBEXPR,/EXTRACT))[1]
  end
  close, lun
  free_lun, lun

  return, var_list

end

