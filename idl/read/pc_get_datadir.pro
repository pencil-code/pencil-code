;
; $Id$
;
; pc_get_datadir
;
; Returns: default data directory; or, if exists, directory from 'datadir.in'.
; 
; 03-Aug-2007 (Anders & Chao-Chin)
; 16-Sep-2015/PABourdin: rewritten
;
function pc_get_datadir

	; Default value of datadir.
	default, data, 'data'
	default, datadir_in, 'datadir.in'

	if (file_test (datadir_in)) then begin
		openr, lun, datadir_in, /get_lun
		readf, lun, data
		close, lun
		free_lun, lun
	endif

	return, data
end
