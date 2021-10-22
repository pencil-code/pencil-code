; Returns: default data directory; or, if exists, directory from 'datadir.in'.
; 
; 03-Aug-2007 (Anders & Chao-Chin)
; 16-Sep-2015/PABourdin: rewritten
;  7-Sep-2021/MR: returns now absolute pathname

function pc_get_datadir, datadir

COMPILE_OPT IDL2, HIDDEN

      if is_defined(datadir) then $
        if strtrim(datadir,2) ne '' then data=datadir

      if not is_defined(data) then begin

	; Default value of datadir.
	data='data'
	datadir_in='datadir.in'

	if (file_test (datadir_in)) then begin
          on_ioerror, ioerr
  	  openr, lun, datadir_in, /get_lun
	  readf, lun, data
ioerr:
	  free_lun, lun
          on_ioerror, NULL
	endif
        if strtrim(data,2) eq '' then data='data'

      endif

      data=strtrim(data,2)
      if strpos(data,'~') eq 0 then data=getenv('HOME')+strmid(data,1) $
      else if strmid(data,0,1) ne '/' then begin
        cd, current=wd
        data=strtrim(wd,2)+'/'+data
      endif
;
;  If no slash at end of string, append it.
;
      if stregex(data,'.*\/$') eq -1 then data=data+'/'
      return, data 

end
