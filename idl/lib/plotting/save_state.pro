;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   save_state.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   26-Sep-2000
;;;  Version: 0.2
;;;
;;;  Description:
;;;    Save the graphics and layout state (environment variables !p,
;;;    !d, !x, !y, !z, color table) onto a stack
;;;  Usage:
;;;    save_state()
;;;    !p.multi = .. & !x.range = .. & [etc..]
;;;    plot,...
;;;    restore_state()
;;;  Limitations: Data do not persist after `.rnew'


pro save_state, VERBOSE=verb, NOOP=noop

  common _state, pstack,dstack,rgbstack,xstack,ystack,zstack,depth,bottom

  default, depth, 0

  tvlct, r,g,b, /GET
  rgb = [[r],[g],[b]]
  if (not keyword_set(noop)) then begin
    if (depth le 0) then begin  ; Lowest level
      pstack = [!p]
      dstack = [!d]
      rgbstack = [rgb]
      xstack = [!x]
      ystack = [!y]
      zstack = [!z]
    endif else begin
      pstack = [pstack,!p]
      dstack = [dstack,!d]
      rgbstack = [[[rgbstack]], [[rgb]]]
      xstack = [xstack,!x]
      ystack = [ystack,!y]
      zstack = [zstack,!z]
    endelse
    depth = depth+1
    bottom = 0
  endif

  if (keyword_set(verb)) then begin
    print,'Save_state: level = ' + strtrim(depth,2)
  endif

end
; End of file save_state.pro
