;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   flatten_strings.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   05-Oct-2002
;;;
;;;  Description:
;;;    Collapse contents of string array into a flat string
;;;  Arguments:
;;;    STRARR  -- string array
;;;  Return value:
;;;    flattened string
;;;  Key words:
;;;    NEWLINE  -- if true, join array elements by newlines
;;;    FINAL_NL -- if true, append a final newline with /NEWLINES

function flatten_strings, strarr, $
                          NEWLINES=nl, FINAL_NL=final_nl, $
                          HELP=help

  compile_opt idl2, hidden

  default, help,     0
  default, nl,       0
  default, final_nl, 0

  if (help) then begin
    print, extract_help('flatten_strings')
    retall
  endif

  if (nl) then begin
    glue = string(10B)
  endif else begin
    glue = ''
  endelse

  strarr1 = strarr              ; copy to avoid overwriting of original
  res = strarr1[0]

  nn=n_elements(strarr)

  while (nn gt 1) do begin
    strarr1 = strarr1[1:*]
    res = res + glue + strarr1[0]
    nn=n_elements(strarr1)
  endwhile

  if (final_nl) then res = res + glue

  return, res

end
; End of file flatten_strings.pro
