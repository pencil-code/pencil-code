;;;;;;;;;;;;;;;;;
;;;  psl.pro  ;;;
;;;;;;;;;;;;;;;;;

;;;   Switch output device to PostScript in landscape format
;;;   Usage:
;;; IDL> psl [,/PSFONTS] [, FILENAME=..]
;;; IDL> plot, ...
;;; IDL> pse

pro psl, $
         _EXTRA=extra
  on_error,2
  psa, _EXTRA=extra, /LANDSCAPE, XSIZE=23, YSIZE=18
end
