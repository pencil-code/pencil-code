;;;;;;;;;;;;;;;;;
;;;  psa.pro  ;;;
;;;;;;;;;;;;;;;;;

;;; Author:  wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;; $Date: 2007-05-17 23:33:39 $
;;; $Revision: 1.7 $

;;;   Switch output device to PostScript
;;;   Usage:
;;; IDL> psa [...]
;;; IDL> plot, ...
;;; IDL> pse
;;;
;;; Special options are
;;;    FILENAME=    for a file name other than 'idl.ps'
;;;   /LANDSCAPE    for landscape format (_not_ seascape as per IDL default,
;;;                 but this produces an inverted bounding box, which
;;;                 confuses some viewers and applications)
;;;   /SEASCAPE_OK  for seascape format (annoying, but has standard
;;;                 bounding box)
;;;   /FULLPAGE     for portrait format covering the full page
;;;   /NOPSFONTS or /NO_PS_FONTS   for using vector fonts instead of
;;;                                hardware (PostScript) fonts
;;;    THICKNESS=   for setting the line thickness (default is 3)
;;;    /QUIET       for suppressing warnings
;;;
;;; All other options (including /COLOR) are passed on to the DEVICE procedure
;;;

pro psa, $
         NOPSFONTS=nops, NO_PS_FONTS=no_ps, $
         FILENAME=filename, $
         LANDSCAPE=landscape, $
         FULLPAGE=fullpage, $
         THICKNESS=thick, $
         BITS_PER_PIXEL=bitdepth, $
;         COLOR=color, $
         SEASCAPE_OK=seascape, $
         QUIET=quiet, $
         _EXTRA=extra
;; Key word NOPSFONTS or NO_PS_FONTS activates vector fonts
  ON_ERROR, 2

  ; COMMON _ps1,_oldthick,_fname,_olddev
  ;; Why not just save all of !p, !d, !x, !y, !z?
  COMMON _ps1,_fname,_oldsysvars

  if (n_elements(quiet) eq 0) then quiet=0

  if (n_elements(seascape) eq 0) then seascape=0
  if (seascape ne 0) then landscape=1

  ;; Lengths (in cm) for LANDSCAPE (use overlap of A4 and letter to
  ;; make sure this prints correctly)
  paperwidth  = 21.00
  paperheight = 27.94
  margin = 1.5
  width  = paperwidth-2*margin
  height = paperheight-2*margin

  ; _olddev=!d.name
  _oldsysvars = { p: !p, d: !d, x: !x, y: !y, z: !z }
  SET_PLOT,'ps'

  IF (N_ELEMENTS(bitdepth) EQ 0) THEN bitdepth=8
  ;; We now unconditionally set bits_per_pixel to 8.
  ;; The /COLOR keyword is still required, because it is handed down to
  ;; DEVICE via _EXTRA
  if (running_gdl() eq 0) then device, BITS_PER_PIXEL=8

  ; ;; If /COLOR keyword is given, set BITS_PER_PIXEL=8 or we will get
  ; ;; only 64 colours (i.e. light yellow instead of white)
  ; ;; NB: We don't want to absorb the /COLOR keyword here since
  ; ;; otherwise later calls to psa would reset color, i.e.
  ; ;;   psa, /color & .. & psa
  ; ;; would make the second plot black an white (which might be a good
  ; ;; thing, but is not currently so). Hence, we peek into the _EXTRA
  ; ;; structure to check for /COLOR.
  ; if (n_elements(extra) gt 0) then begin
  ;   if (has_tag (extra, 'color')) then begin
  ;     if (extra.color) then device, BITS_PER_PIXEL=8
  ;   endif
  ; endif

  IF (N_ELEMENTS(filename) EQ 0) THEN filename='idl.ps'
  IF (N_ELEMENTS(landscape) EQ 0) THEN BEGIN ; portrait
    IF (N_ELEMENTS(fullpage) EQ 0) THEN BEGIN
      DEVICE, FILENAME=filename, _EXTRA=extra, /PORTRAIT
    ENDIF ELSE BEGIN
      DEVICE, FILENAME=filename, _EXTRA=extra, /PORTRAIT, $
          XSIZE=width, YSIZE=height, $
          XOFFSET=margin, YOFFSET=margin
    ENDELSE
  ENDIF ELSE BEGIN              ; landscape
    DEVICE, FILENAME=filename, _EXTRA=extra, LANDSCAPE=landscape
    if (landscape) then begin
      ;; we want landscape, not seascape -- but there are potential problems
      ;
      if (seascape eq 0) then begin
        scalefact = -1.
        xoff      = paperwidth-margin
        yoff      = margin
        if (quiet eq 0) then begin
          print, 'WARNING:', $
              ' Negative SCALE_FACTOR produces inverted bounding box, which'
          print, 'WARNING:', $
              ' some viewers/applications cannot handle.'
          print, 'WARNING:', $
              ' Use the /SEASCAPE_OK option to avoid this (and get seascape).'
        endif
      endif else begin
        scalefact = 1.
        xoff      = margin
        yoff      = paperheight-margin
      endelse
      ;
      device, $
          XSIZE=height, YSIZE=width, $
          XOFFSET=xoff, YOFFSET=yoff, SCALE_FACTOR=scalefact
    endif
  ENDELSE
  IF (NOT (KEYWORD_SET(nops) OR KEYWORD_SET(no_ps))) THEN ps_fonts
  IF (N_ELEMENTS(thick) EQ 0) THEN thick=3
;; Remember the original thicknesses
  IF (n_elements(_oldthick) le 1) THEN BEGIN
    _oldthick = [1, !P.CHARTHICK, !P.THICK, !X.THICK, !Y.THICK, !Z.THICK]
  ENDIF
  !P.CHARTHICK=thick
  !P.THICK=thick
  !X.THICK=!P.THICK
  !Y.THICK=!P.THICK
  !Z.THICK=!P.THICK
  _fname = filename

end
