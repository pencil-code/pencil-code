;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   extract_help.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   01-Apr-2005
;;;
;;;  Description:
;;;    Extract help text (starting `;;; Description:' til next noncomment line)
;;;    from the current IDL function or procedure.
;;;  Usage:
;;;    function toto, .., HELP=help / pro toto, ..., HELP=help
;;;      if (keyword_set(help)) then extract_help, 'toto'
;;;      [...]
;;;    end
;;;  To do:
;;;    - Manually scan the !path, because HELP, /SOURCE does not show
;;;      routines that contain `compile_opt idl2, hidden' to avoid
;;;      compilation message

pro extract_help, name, $
                  DEBUG=debug, HELP=help

  if (keyword_set(help)) then extract_help, 'extract_help'

  default, debug, 0

  ;; Traverse through !path and find file
  paths = str_sep(!path, ':')
  Np = n_elements(paths)
  for i=0L, Np-1 do begin
    fname = paths[i] + '/' + name + '.pro'
    if (file_test(fname)) then begin
      i = Np                    ; found => stop looping
    endif else begin
      fname = ''                ; otherwise, clear fname
    endelse
  endfor
  ;
  if (fname eq '') then message, "Couldn't find <" + name + '.pro> in !path' 

  if (debug) then print, 'EXTRACT_HELP: fname = <' + fname + '>'

  fnameq = '"'+fname+'"' ; quote fname to allow silly spaces in file names
  cmd = 'perl -ne "print if /^;;;\s*(Description|Usage):/../^[^;]/;" ' $
      + fnameq
  if (debug) then print, 'EXTRACT_HELP: cmd = <' + cmd + '>'
  spawn, cmd, help_text

  header = ';;;  ' + ['Name: ', '  '+name] 

  print, flatten_strings([header,help_text],/NEWLINES)
  retall

end
; End of file extract_help.pro
