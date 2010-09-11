;;;;;;;;;;;;;;;;;;;;;;
;;;  ps_fonts.pro  ;;;
;;;;;;;;;;;;;;;;;;;;;;

;;; Modified by wd (Dobler@kis.uni-freiburg.de)
;;; Date:    20-Nov-2001
;;; Version: 1.2

;;;
;;; Substitute the IDL-vector-fonts by PostScript-fonts
;;; ---------------------------------------------------
;;;
;;; this is done in a way, that the numbers !XX 
;;; result in both cases in the same style
;;;
;;; Helvetica              !3        sans serif
;;; NewCentury Schoolbook  !5        normal text (variant)
;;; Times-Roman            !6        for normal text
;;; Symbol                 !7        greek alphabet
;;; Times-Italic           !8        for math symbols
;;; Zapf-Chancery          !12       calligraphic
;;; NCSchoolb.-Italic      !15       math symbols (variant)
;;; Times-Roman-Bold       !17       bold text
;;; Times-Bold-Italic      !18       bold italics
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro ps_fonts
  COMMON _ps_fonts, _old_P_FONT
  if (n_elements(_old_P_FONT) le 1) then begin
    _old_P_FONT = [1,!P.FONT]
  endif

  !P.FONT = 0
  if (running_gdl() eq 0) then begin
    device, /HELVETICA, ISOLATIN1=1,			font_index=3
    device, /SCHOOLBOOK, ISOLATIN1=1,			font_index=5
    device, /TIMES, ISOLATIN1=1,				font_index=6
    device, /SYMBOL,					font_index=7
    device, /TIMES, /ITALIC, ISOLATIN1=1,			font_index=8
    device, /ZAPFCHANCERY, /MEDIUM, /ITALIC,ISOLATIN1=1,	font_index=12
    device, /SCHOOLBOOK, /ITALIC, ISOLATIN1=1,		font_index=15
    device, /TIMES, /BOLD, ISOLATIN1=1,			font_index=17
    device, /TIMES, /BOLD, /ITALIC, ISOLATIN1=1,		font_index=18
  endif
return
end
