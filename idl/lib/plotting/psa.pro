;;;;;;;;;;;;;;;;;
;;;  psa.pro  ;;;
;;;;;;;;;;;;;;;;;

;;; Author:  wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;; $Date: 2003-11-24 16:35:01 $
;;; $Revision: 1.1 $

;;;   Switch output device to PostScript
;;;   Usage:
;;; IDL> psa [...]
;;; IDL> plot, ...
;;; IDL> pse
;;;
;;; Special options are
;;;   /LANDSCAPE  for landscape format (_not_ seascape as per IDL default)
;;;    FILENAME=  for a file name other than 'idl.ps'
;;;   /FULLPAGE   for portrait format covering the full page
;;;   /NOPSFONTS or /NO_PS_FONTS   for using vector fonts instead of
;;;                                hardware fonts
;;;    THICKNESS= for setting the line thickness (default is 3)
;;; All other options are passed on to the DEVICE procedure
;;;

pro psa, $
         NOPSFONTS=nops, NO_PS_FONTS=no_ps, $
         FILENAME=filename, $
         LANDSCAPE=landscape, $
         FULLPAGE=fullpage, $
         THICKNESS=thick,$
         _EXTRA=extra
;; Key word NOPSFONTS or NO_PS_FONTS activates vector fonts
  ON_ERROR,2

  COMMON _ps1,_oldthick,_fname,_olddev

  ;; Lengths (in cm) for LANDSCAPE
  paperwidth  = 21.00
  paperheight = 29.70
  margin = 1.5
  width  = paperwidth-2*margin
  height = paperheight-2*margin

  _olddev=!d.name
  SET_PLOT,'ps'
  IF (N_ELEMENTS(filename) EQ 0) THEN filename='idl.ps'
  IF (N_ELEMENTS(landscape) EQ 0) THEN BEGIN ; portrait
    IF (N_ELEMENTS(fullpage) EQ 0) THEN BEGIN
      DEVICE, FILENAME=filename, _EXTRA=extra, /PORTRAIT
    ENDIF ELSE BEGIN
      DEVICE, FILENAME=filename, _EXTRA=extra, /PORTRAIT, $
          XSIZE=18, YSIZE=26, XOFFSET=1, YOFFSET=2.2
    ENDELSE
  ENDIF ELSE BEGIN              ; landscape
    DEVICE, FILENAME=filename, _EXTRA=extra, LANDSCAPE=landscape
    if (landscape) then begin   ; we want landscape, not seascape
      device, $
          XSIZE=height, YSIZE=width, $
          XOFF=paperwidth-margin, YOFF=margin, scale_factor=-1.
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
