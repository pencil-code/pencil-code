PRO LINE_KEY, X , Y , LABEL , LINESTY , LINELENGTH = LINELENGTH $
            , THICK = THICK , SYMBOL = SYMBOL , CHARSIZE = CHARSIZE , VSHIFT = VSHIFT, _extra=extra 
;+
; NAME:
;       LINE_KEY
;
; PURPOSE:
;       Put a key for a curve in a plot, together with its corresponding 
;       linestyle or plot symbol.
;
; CALLING SEQUENCE:
;       LINE_KEY, X , Y , LABEL , LINESTY [ , LINELENGTH = LINELENGTH
;               , /SYMBOL , CHARSIZE = CHARSIZE ]
;
; INPUTS:
;       X               X-position of label in normalized coordinates
;
;       Y               Y-position of label in normalized coordinates
;
;       LABEL           The label, a string scalar.
;
;       LINESTY         The linestyle of the curve for which to put a key.
;
; OPTIONAL KEYWORDS:
;       LINELENGTH      Lenght of line-key in fraction of plot window.
;                       Defaults to 0.1.
;
;       THICK           The line thickness.
;
;      SYMBOL           If set, the plot symbol corresponding to SYMBOL
;                       is plotted instead of a line for a positive value
;                       a line with linestyle LINESTY *and* symbols 
;                       according to SYMBOL otherwise
;
;       CHARSIZE        The charsize used for the label.
;
;       VSHIFT          Vertical shift of the label text, vshift*charsize downwards.
;
; OUTPUTS:
;       Result: real or double precision scalar or array
;
; RESTRICTIONS
;       A plot window must be defined in advance.
;
; MODIFICATION HISTORY:
;       WRITTEN, Michael Andersen, Nordic Optical Telescope, November, 1997
;-

  ON_ERROR, 2

; Check for right number of input parameters 

  IF( N_PARAMS() NE 3 AND N_PARAMS() NE 4 ) THEN $
  MESSAGE, 'Called with wrong number of parameters'
  IF( N_PARAMS() EQ 3 AND NOT KEYWORD_SET( SYMBOL ) ) THEN $
  MESSAGE, 'Linestyle or SYMBOL keyword not specified'

; Calculate where to put the line and label in device coordinates.

  x_win_size    = !P.CLIP(2) - !P.CLIP(0)
  y_win_size    = !P.CLIP(3) - !P.CLIP(1)

  IF( NOT KEYWORD_SET( LINELENGTH ) ) THEN linelength = 0.1
  IF( NOT KEYWORD_SET( THICK      ) ) THEN thick      = 1.0
  IF( NOT KEYWORD_SET( CHARSIZE   ) ) THEN charsize   = 1.0
  IF(     KEYWORD_SET( SYMBOL     ) ) THEN $
    IF( SYMBOL GT 0 ) THEN linelength = 0.0    
  IF( NOT KEYWORD_SET( VSHIFT     ) ) THEN vshift     = 0.02

  lx1           =   x                        * x_win_size + !P.CLIP(0)
  lx2           = ( x + linelength         ) * x_win_size + !P.CLIP(0)
  xl            = ( x + linelength + 0.015 ) * x_win_size + !P.CLIP(0)
  ly            =   y                        * y_win_size + !P.CLIP(1)
  yl            = ( y - vshift * charsize    ) * y_win_size + !P.CLIP(1)

  IF( NOT KEYWORD_SET( SYMBOL ) ) $
    THEN PLOTS, [ lx1 , lx2 ] , [ ly , ly ] $
              , LINESTY = linesty , THICK = thick , /DEVICE, _extra=extra $
    ELSE PLOTS, [ lx1 , lx2 ] , [ ly , ly ] $
              , LINESTY = linesty, PSYM = SYMBOL  , /DEVICE, _extra=extra
  XYOUTS, xl , yl , label , CHARSIZE = charsize , /DEVICE

END
