;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   contourfill.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   21-Nov-2000
;;;  Version: 0.5
;;;  Description:
;;;   A modified contour subroutine.
;;;   - Does not require `reform' on data[*,0,*]
;;;   - Does coloured filling by default
;;;   - Can overplot the grid
;;;   - Can add colorbar
;;;
;;;   New keywords:
;;;     /GRID     --  overplot the grid points
;;;     COLORBAR  --  add color bar (1 or 'bottom': place at bottom;
;;;                  'right': place to the right
;;;
;;;   Modifications:
;;;   13-apr-16/MR: modifications for use of irregular/triangulated grid analogous to contour

pro contourfill, z, x, y, $
                 NLEVELS=nlevels, LEVELS=levels_, FILL=fill, GRID=grid, $
                 COLORBAR=colbar, DEBUG=debug, $
                 C_COLORS=c_colors, _EXTRA=_extra

  default, nlevels, 60
  default, fill, 1
  default, debug, 0

  irreg = is_in(tag_names(_extra), 'triangulate', /abbrev) or $
          is_in(tag_names(_extra), 'irregular', /abbrev)

; Contour does not reset !z.type to zero if called with /ZLOG. This is
; silly (although it allows subsequent calls with /over to work
; correctly [and in most cases unexpectedly, so we just do not care]). 

  oldztype = !z.type            ; Save this to restore later

  if (debug) then $
    print, FORMAT='(A, 10(I3))', 'CONTOURFILL: Initial !p.multi = ', !p.multi

  if (n_elements(colbar) gt 0) then begin
    ;; Determine data type
    s = size(colbar)
    type = s[s[0]+1]
    if (type eq 7) then begin   ; a string
      case (colbar) of
        'none':   colbar=0
        'bottom': colbar=1
        'right':  colbar=2
        else: message, 'COLBAR must be either numerical or one of "none", "bottom", "right"'
      endcase
    endif
  endif else begin
    colbar = 0
  endelse
;
  array = reform(z)
  if ((size(array))[0] gt 2) then begin
    message, /INFO, 'Argument is 3-d or higher. Plotting first section only'
    array = array[*,*]
  endif
;
  if ~irreg then begin
    s = size(array)
    if (n_elements(x) eq 0) then x = indgen(s[1])
    if (n_elements(y) eq 0) then y = indgen(s[2])
  endif
;
  xmarg = !x.margin             ; Store and reset later
  ymarg = !y.margin
  
  if (colbar eq 1) then !y.margin = ymarg + [3,0]
  if (colbar eq 2) then !x.margin = xmarg + [0,8]

  ;; non-trivial to keep IDL from using very first or very last color:
  ;; currently doesn't work for z-logarithmic plots
;
  if (n_elements(levels_) gt 0) then levels=levels_
  default, levels, linspace(minmax(array,/NAN),nlevels,GHOST=1,/UNIQUE)

  contour, array, x, y, LEVELS=levels, FILL=fill, _EXTRA=_extra
;
  if (keyword_set(grid)) then begin
    if ~irreg then begin
      if (debug ne 0) then print,'s = ', s
      plots, spread(x, 1, s[2]), spread(y, 0, s[1]), PSYM=3, NOCLIP=0
    endif else $
      oplot, x, y, PSYM=3, COLOR=0 
  endif

  if (colbar ne 0) then begin
    dmin = min(array)
    dmax = max(array)
    case (colbar) of
      1: begin                  ; bottom
        colorbar, DIVISIONS=6, MIN=dmin, MAX=dmax, POSITION=0.07
      end
      2: begin                  ; right
        colorbar, DIVISIONS=6, MIN=dmin, MAX=dmax, /VERTICAL, POSITION=0.93
      end
    endcase        
  endif

  !x.margin = xmarg
  !y.margin = ymarg

  if (debug) then begin
    print, FORMAT='(A, 10(I3))', $
           'CONTOURFILL: Final !p.multi   = ', !p.multi
    print, '------------------------------------------------'
  endif

  !z.type = oldztype            ; Restore saved.

end
; End of file contourfill.pro
