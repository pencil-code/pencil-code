;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   pc_convert_vars_struct.pro     ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  $Id$
;;;
;;;  Description:
;;;   Converts a given vars structure into an array and its corresponding tags.
;;;

; Conversion of vars structure into an array.
function pc_convert_vars_struct, vars, varcontent, tags

  if (size (vars, /type) ne 8) then message, "The given data is not a structure."

  if (n_elements (varcontent) eq 0) then varcontent = pc_varcontent

  ; Create array out of given structure and pass recursively computed results
  totalvars = (size (varcontent))[1]
  tag_str = ''
  arr_str = ''
  for iv=0L, totalvars-1L do begin
    if (varcontent[iv].idlvar eq "uu") then begin
      tag_str += ', uu:[' + strtrim (iv, 2) + ',' + strtrim (iv+1, 2) + ',' + strtrim (iv+2, 2) + ']'
      tag_str += ', ux:' + strtrim (iv, 2) + ', uy:' + strtrim (iv+1, 2) + ', uz:' + strtrim (iv+2, 2)
      arr_str += ' & array[*,*,*,tags.uu] = vars.uu'
    endif else if (varcontent[iv].idlvar eq "aa") then begin
      tag_str += ', aa:[' + strtrim (iv, 2) + ',' + strtrim (iv+1, 2) + ',' + strtrim (iv+2, 2) + ']'
      tag_str += ', ax:' + strtrim (iv, 2) + ', ay:' + strtrim (iv+1, 2) + ', az:' + strtrim (iv+2, 2)
      arr_str += ' & array[*,*,*,tags.aa] = vars.aa'
    endif else begin
      tag_str += ', ' + varcontent[iv].idlvar + ':' + strtrim (iv, 2)
      arr_str += ' & array[*,*,*,tags.' + varcontent[iv].idlvar + '] = vars.' + varcontent[iv].idlvar
    endelse
    ; For vector quantities skip the required number of elements of the f array.
    iv += varcontent[iv].skip
  endfor
  tag_str = 'tags = { ' + strmid (tag_str, 2) + ' }'
  if (execute (tag_str) ne 1) then message, 'Error executing: ' + tag_str

  index = where (strcmp (tag_names (vars[0]), varcontent[0].idlvar, /fold_case))
  s = size (vars[0].(index), /structure)
  if (s.type eq 4) then begin
    array = fltarr (s.dimensions[0], s.dimensions[1], s.dimensions[2], totalvars)
  end else begin
    array = dblarr (s.dimensions[0], s.dimensions[1], s.dimensions[2], totalvars)
  end
  arr_str = strmid (arr_str, 3)
  if (execute (arr_str) ne 1) then message, 'Error executing: >' + arr_str + '<'

  return, array

end

