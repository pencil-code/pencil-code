;;;;;;;;;;;;;;;;;;;;;;;
;;;   wdtvscl.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   16-Sep-1998
;;;
;;;  Description:
;;;   A variant of TVSCL that accepts the keyword BOTTOM in addition
;;;   to TOP.
;;;   With ABSOLUT set, scale symmetric around zero, so a value of
;;;   zero gets always the middlemost colour.

PRO wdtvscl, data, x, y, channel, $
              TOP=top, BOTTOM=bottom, $
              ABSOLUT=absolut, $
              _EXTRA=extra
  ON_ERROR, 2
  IF (N_ELEMENTS(top) EQ 0) THEN top = !D.TABLE_SIZE-1
  IF (N_ELEMENTS(bottom) EQ 0) THEN bottom = 0
  IF (N_ELEMENTS(absolut) EQ 0) THEN absolut = 0
  IF (absolut EQ 0) THEN BEGIN
    dat_min = MIN(data)
    dat_max = MAX(data)
  ENDIF ELSE BEGIN
    dat_min = -MAX(abs(data))
    dat_max = MAX(abs(data))
  ENDELSE

;;; Make a homogeneous field as bright as possible
;  IF ((dat_min EQ dat_max) AND (!D.NAME='PS')) THEN BEGIN
;    dat_min = dat_min - 100
;    dat_max = dat_max
;    print,'dat_interv = ', dat_min, dat_max
;  ENDIF
  IF (dat_min EQ dat_max) THEN BEGIN
    IF (dat_min EQ 0) THEN BEGIN
      dat_min = -1
      dat_max = 1
    ENDIF ELSE BEGIN
      d_dat = abs(dat_min)
      dat_min = dat_min - d_dat
      dat_max = dat_max + d_dat
    ENDELSE
  ENDIF
;
  IF (N_PARAMS() EQ 1) THEN BEGIN ; MY_TVSCL, IMAGE
    TV, bottom + (top-bottom)*(data-dat_min)/(dat_max-dat_min), $
        _EXTRA=extra
  ENDIF
;
  IF (N_PARAMS() EQ 2) THEN BEGIN ; MY_TVSCL, IMAGE, POSITION
    TV, bottom + (top-bottom)*(data-dat_min)/(dat_max-dat_min), $
        x, $
        _EXTRA=extra
  ENDIF
;
  IF (N_PARAMS() EQ 3) THEN BEGIN ; MY_TVSCL, IMAGE, X, Y
    TV, bottom + (top-bottom)*(data-dat_min)/(dat_max-dat_min), $
        x, y, $
        _EXTRA=extra
  ENDIF
;
  IF (N_PARAMS() EQ 4) THEN BEGIN ; MY_TVSCL, IMAGE, X, Y, CHANNEL
    TV, bottom + (top-bottom)*(data-dat_min)/(dat_max-dat_min), $
        x, y, channel, $
        _EXTRA=extra
  ENDIF

END
; End of file wdtvscl.pro
