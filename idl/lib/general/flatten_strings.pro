;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   flatten_strings.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   05-Oct-2002
;;;
;;;  Description:
;;;   Collapse contents of string array into a flat string

function flatten_strings, strarr

  strarr1 = strarr              ; copy to avoid overwriting of original
  res = strarr1[0]
  while ((size(strarr1))[1] gt 1) do begin
    strarr1 = strarr1[1:*]
    res = res + strarr1[0]
  endwhile

  return, res

end
; End of file flatten_strings.pro
