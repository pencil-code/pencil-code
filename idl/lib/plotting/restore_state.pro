;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   restore_state.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   26-Sep-2000
;;;  Version: 0.21
;;;
;;;  Description:
;;;    Restore the graphics and layout state saved by `save_state()'
;;;  Usage:
;;;    save_state()
;;;    !p.multi = .. & !x.range = .. & [etc..]
;;;    plot,...
;;;    restore_state()
;;;  Limitations: Data do not persist after `.rnew'

pro restore_state, VERBOSE=verb, NOOP=noop, FULL=full

  common _state, pstack,dstack,rgbstack,xstack,ystack,zstack,depth,bottom

  if (n_elements(depth) eq 0) then begin
    if (keyword_set(verb)) then print,'Restore_state: No information stored'
    return                      ; save_state has never been called
  endif

  default, full, 0

  if (not keyword_set(noop)) then begin
    if (bottom) then begin
      message, 'Already at lowest level', /INFO
    endif else begin
      if (full) then depth=1    ; does not clear pstack etc,
                                ; but next save_state will
      ;;; !d variable
      dstruct = dstack[depth-1] ; (You can't write to !d)
      ;; set_plot will reset some things (like background colour), so
      ;; only call when needed:
      if (!d.name ne dstruct.name) then set_plot, dstruct.name
      if (!d.name eq 'X') then wset, dstruct.window
      ;;; !p variable
      !p = pstack(depth-1)
      ;;; colour table
      rgb = rgbstack[*,*,depth-1]
      tvlct, rgb
      ;;; !x, !y, !z variables
      !x = xstack[depth-1]
      !y = ystack[depth-1]
      !z = zstack[depth-1]
      ;; Some fields need extra adjustment
;; Why this one: ?
;    !p.multi[0]=0               ; Next plot overwrites
      ;; .. in particular all the fields of !d
      if (depth gt 1) then begin ; Don't remove entries at lowest level
        pstack = pstack[0:depth-2]
        dstack = dstack[0:depth-2]
        rgbstack = rgbstack[*,*,0:depth-2]
        xstack = xstack[0:depth-2]
        ystack = ystack[0:depth-2]
        zstack = zstack[0:depth-2]
      endif else begin
        bottom = 1
      endelse
      depth = depth-1
    endelse
  endif

  if (keyword_set(verb)) then begin
    print,'Restore_state: level = ' + strtrim(depth,2)
  endif

end
; End of file restore_state.pro
