; Returns: default data directory; or, if exists, directory from 'datadir.in'.
; 
; 03-Aug-2007 (Anders & Chao-Chin)
; 16-Sep-2015/PABourdin: rewritten

function pc_get_datadir, datadir

COMPILE_OPT IDL2, HIDDEN

      if is_defined(datadir) then begin

        datadir=strtrim(datadir,2)
        if strpos(datadir,'~') eq 0 then datadir='$HOME'+strmid(datadir,1)
        return, datadir 

      endif else begin

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

      endelse
end
