;;;;;;;;;;;;;;;;;;;;;;;
;;;   warp_ct.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   20-Apr-2005
;;;
;;;  Description:
;;;    Distort a colour table for enhanced contrast, b/w printing, etc.
;;;    Currently, cubic distortion relative to the midpoint of the table.
;;;
;;;  Arguments:
;;;    c0, c1 -- lowest and highest colour of new colour table (as
;;;              index values in original table)
;;;  Keywords:
;;;    SLOPE -- 1 for linear (default), large for high-contrast, low
;;;             for low-contrast warping.
;;;             NB: For  slope > 1.5 , colours at the edge are
;;;             truncated (as the cubic profile will overshoot)
;;;    BW    -- flag for making first + last colour !p.color and !p.background
;;;    RESET -- flag for restoring original colormap
;;;
;;;  Examples:
;;;    warp_ct, 129, 255            ; only use upper half of colours
;;;    warp_ct, 70,, 255            ; lighter colours, linear warping
;;;    warp_ct, 70,, 255, SLOPE=10. ; lighter colours, high contrasts
;;;    warp_ct, /RESTORE            ; back to original colour table

pro warp_ct, c0, c1, $
             SLOPE=slope, $
             BW=bw, RESTORE=restore, $
             HELP=help

  if (keyword_set(help)) then extract_help, 'warp_ct'

  common _warp_ct, oldct, level

  default, slope,   1.
  default, bw,      0
  default, restore, 0
  default, level,   0

  ;; Store current colour table
  tvlct, /GET, red, green, blue
  Ncol = !d.table_size
  color = indgen(Ncol)

  if (restore) then begin
    if (level gt 0) then begin
      tvlct, oldct.r, oldct.g, oldct.b 
      level = 0
    endif else begin
      print, 'Already at lowest level..'
    endelse
  endif else begin
    ;; Warp colour table
    if (level eq 0) then begin
      ;; Store current table
      oldct = {r: red, g: green, b: blue}
      level = level + 1
    endif
    Ncol_2 = (Ncol-1)/2.
    zeta = (color-Ncol_2)*1./Ncol_2 ; in [-1,1]
    zeta_new = slope*zeta + (1-slope)*zeta^3
    idx = (c0+c1)/2. + (c1-c0)/2.*zeta_new
    idx = ((idx > c0) < c1)

print, 'minmax(zeta)=', minmax(zeta)
print, 'minmax(zeta_new)=', minmax(zeta_new)
print, 'c0, c1     =', c0, c1
print, 'minmax(idx)=',minmax(idx)

    ; r1 = red  [idx]
    ; g1 = green[idx]
    ; b1 = blue [idx]
    r1 = interpol(red*1.  , indgen(Ncol), idx)
    g1 = interpol(green*1., indgen(Ncol), idx)
    b1 = interpol(blue*1. , indgen(Ncol), idx)

    ;; Fix lowest and highest level, if requested
    if (bw) then begin
      fg = 0
      bg = 255
      r1[     0]=fg & g1[     0]=fg & b1[     0]=fg
      r1[Ncol-1]=bg & g1[Ncol-1]=bg & b1[Ncol-1]=bg
    endif

    ;; Now load new color table
    tvlct, r1, g1, b1

  endelse


end
; End of file warp_ct.pro
