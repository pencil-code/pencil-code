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

  if (n_elements(help) ne 0) then begin
    print, "DEFAULT:  usage:"
    print, "  default, var, val"
    print, "  default, STRUCT=struc, 'slot', val"
    print, "  default, STRUCT=struc, ['slot1',..,'slotN'], val"
    print, "  default, STRUCT=struc, ['slot1',..,'slotN'], [val1,...,valN]"
    return
  endif

  if (n_elements(debug) eq 0) then debug=0

  if (n_elements(struct) eq 0) then begin ;; simple variable
    if n_elements(var) eq 0 then var=val
  endif else begin ;; structure slot
    ; expand val if necessary:
    nvar = n_elements(var)
    if ((nvar gt 1) and (n_elements(val) eq 1)) then $
        val = replicate(val,nvar)
    nval = n_elements(val)

    if (debug) then print, 'var, val = <', var, '>, <', val, '>' 
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
      if (idx lt 0) then begin  ; tag didn't exist
        if (nvar gt 1) then begin
          ;; several slots --> assign one by one
          cmd = cmd + ", '" + var[i] + "'" + ", val[" + strtrim(i,2) + "]"
        endif else begin
          ;; just one slot --> assign whole variable (could be array)
          cmd = cmd + ", '" + var[i] + "'" + ", val"
        endelse
        nnewtags = nnewtags + 1
      endif
    endfor

if (debug) then print, 'nval,nvar = ', nval, nvar
 
    cmd = cmd + " )"
    if (nnewtags ne 0) then begin
      if (debug) then print, 'Now running: ',cmd
      res = execute(cmd)
    endif

;;    if (n_elements(struct.var) eq 0) then struct.var=val
  endelse

end
