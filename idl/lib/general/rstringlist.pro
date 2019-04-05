  function rstringlist, file, silent=silent
;
; Reads a list of strings in file (one or more per line) into returned array.
; Returns one-empty-element array if file can't be read or is contains only empty strings.
;
    on_ioerror, openerr
    openr, lun, file, /get_lun
    on_ioerror, NULL
    
    str='' & str_chain=''

    while ~ eof(lun) do begin
      readf, lun, str
      str_chain += str+' '
    endwhile

    close, lun
    free_lun, lun

    if strtrim(str_chain,2) eq '' then begin
      if ~keyword_set(silent) then $
        print, 'rstringlist: No strings in file "'+strtrim(file,2)+'"!!!'
      return, ['']
    endif else $
      return, strsplit(str_chain,/extract)

openerr:
    if ~keyword_set(silent) then $
      print, 'rstringlist: Error when opening file "'+strtrim(file,2)+'"!!!'
    return, ['']

  end
