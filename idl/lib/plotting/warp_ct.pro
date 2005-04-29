;;;;;;;;;;;;;;;;;;;;;;;
;;;   warp_ct.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   20-Apr-2005
;;;
;;;  Description:
;;;    Distort a colour table for enhanced contrast, b/w printing, etc.
;;;    Currently supported cubic distortion (using SLOPE) and
;;;    power-law (using EXPONENT), relative to the midpoint of the
;;;    table.
;;;
;;;  Arguments:
;;;    c0, c1   -- lowest and highest colour of new colour table (as
;;;                index values in original table)
;;;  Keywords:
;;;    SLOPE    -- slope for cubic warping.
;;;                < 1 : low-contrast warping,
;;;                = 1 : linear (default),
;;;                > 1 : high-contrast warping.
;;;                NB: For  slope > 1.5 , colours at the edge are
;;;                truncated (as the cubic profile will overshoot)
;;;    EXPONENT -- exponent for power-law warping.
;;;                < 1 : high-contrast warping,
;;;                = 1 : linear (default),
;;;                > 1 : low-contrast warping.
;;;    BW       -- flag for making first + last colour !p.color and
;;;                !p.background
;;;    NSTEPS   -- set number of steps in colour table (effectively
;;;                limiting !d.table_size, but in reality there will
;;;                still be the same number of colours).
;;;                If positive, apply discretization before warping (->
;;;                equally spaced colours in the original table);
;;;                if negative, apply after warping (-> equally spaced
;;;                in new colour table)
;;;    RESET    -- flag for restoring original colormap
;;;
;;;  Examples:
;;;    warp_ct, 129, 255               ; only use upper half of colours
;;;    warp_ct, 70,, 255               ; lighter colours, linear warping
;;;    warp_ct, 70,, 255, SLOPE=10.    ; lighter colours, high contrasts
;;;    warp_ct, 70,, 255, EXPONENT=0.3 ; lighter colours, high contrasts
;;;    warp_ct, 70,, 255, NSTEPS=15    ; no warping, but only 15 colors
;;;    warp_ct, 70,, 255, NSTEPS=-15   ; colorbar will be equidistant
;;;    warp_ct, /RESTORE               ; back to original colour table

pro warp_ct, c0, c1, $
             SLOPE=slope, $
             EXPONENT=exponent, $
             BW=bw, NSTEPS=nsteps, $
             RESTORE=restore, $
             HELP=help

  if (keyword_set(help)) then extract_help, 'warp_ct'

  common _warp_ct, oldct, level

  default, slope,    1.
  default, exponent, 1.
  default, bw,       0
  default, nsteps,   0
  default, restore,  0
  default, level,    0

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
    ;
    if (slope ne 1.) then begin
      zeta = slope*zeta + (1-slope)*zeta^3
    endif else if (exponent ne 1.) then begin
      zeta = sign(zeta)*(abs(zeta)^exponent)
    endif
    idx = (c0+c1)/2. + (c1-c0)/2.*zeta

    ;; Apply nsteps before warping
    if (nsteps gt 0) then begin
      idx = floor(idx*nsteps/Ncol+0.5)*Ncol/nsteps
    endif

    idx = ((idx > c0) < c1)     ; sanitize

    ; r1 = red  [idx]
    ; g1 = green[idx]
    ; b1 = blue [idx]
    r1 = interpol(red*1.  , indgen(Ncol), idx)
    g1 = interpol(green*1., indgen(Ncol), idx)
    b1 = interpol(blue*1. , indgen(Ncol), idx)

    ;; Apply nsteps after warping
    if (nsteps lt 0) then begin
      idx1 = findgen(Ncol)
      idx1 = floor(idx1*abs(nsteps)/Ncol+0.5)*Ncol/abs(nsteps)
      r1 = r1[idx1]
      g1 = g1[idx1]
      b1 = b1[idx1]
    endif

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
