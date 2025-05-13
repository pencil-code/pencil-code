  function rstringlist, file, silent=silent, commchar=commchar
;
; Reads a list of strings in file (one or more per line) into returned array.
; If commchar is given, ignore string swhich begin with one of the characters in commchar.
; Returns one-empty-element array if file can't be read or is contains only empty strings.
;
    on_ioerror, openerr
    openr, lun, file, /get_lun
    on_ioerror, NULL
    
    str='' & str_chain=''

    while ~ eof(lun) do begin
      readf, lun, str
      if keyword_set(commchar) then begin
;
; Ignore string if first non-blank character is in list of comment characters commchar
;
        if stregex(str,'^ *['+strtrim(commchar,2)+']') eq -1 then $
          str_chain += str+' '
      endif else $
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
