;;;;;;;;;;;;;;;;;;;;;;;
;;;   default.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;  $Date: 2004-07-14 19:24:56 $
;;;
;;;  Description:
;;;   Set variable or structure slot to default value unless it
;;;   already exists.
;;;  Usage:
;;;   default, var, val
;;;   default, STRUCT=struc, 'slot', val
;;;   default, STRUCT=struc, ['slot1',...,'slotN'], val
;;;   default, STRUCT=struc, ['slot1',...,'slotN'], [val1,...,valN]
;;;  Examples:
;;;   default, STRUCT=par, 'lchiral', 0L
;;;   default, STRUCT=par, ['lchiral','lcolor'], [0L, 1L]
;;;   default, STRUCT=par, ['lchiral','lcolor'], 0L

pro default, var, val, STRUCT=struct, HELP=help, DEBUG=debug
  compile_opt idl2,hidden

  if (keyword_set (help)) then begin
    print, "DEFAULT:  usage:"
    print, "  default, var, val"
    print, "  default, STRUCT=struc, 'slot', val"
    print, "  default, STRUCT=struc, ['slot1',..,'slotN'], val"
    print, "  default, STRUCT=struc, ['slot1',..,'slotN'], [val1,...,valN]"
    return
  endif

  if (keyword_set (struct)) then begin
    ; structure:
    nvar = n_elements(var)
    nval = n_elements(val)
    if ((nvar gt 1) and (nval eq 1)) then begin
      ; expand val
      val = replicate(val,nvar)
      nval = n_elements(val)
    endif

    if (keyword_set (debug)) then begin
      ; debug output
      print, 'nvar = ', nvar, ', nval = ', nval
      print, 'var = <', var, '>'
      print, 'val = <', val, '>'
    endif

    ; sanity tests
    s = size(var)
    if (s[s[0]+1] ne 7) then $
        message, "First arg must be a string when called with STRUCT"
    if ((nvar gt 1) and (nvar ne nval)) then $
        message, $
        "slot and val arrays must have equal size or val or slot be scalar"

    cmd = "struct = create_struct( struct"
    nnewtags = 0
    for i=0,nvar-1 do begin
      idx = min(where(tag_names(struct) eq strupcase(var[i])))
      if (idx lt 0) then begin
        ; tag didn't exist
        if (nvar gt 1) then begin
          ; several slots --> assign one by one
          cmd += ", '" + var[i] + "'" + ", val[" + strtrim(i,2) + "]"
        endif else begin
          ; just one slot --> assign whole variable (could be array)
          cmd += ", '" + var[i] + "'" + ", val"
        endelse
        nnewtags++
      endif
    endfor
    cmd += " )"
    if (nnewtags ge 1) then begin
      if (not execute(cmd)) then message, 'Could not create structure: ', cmd
    endif

    return
  end

  ; simple variable
  if (size (var, /type) eq 0) then var = val

end
