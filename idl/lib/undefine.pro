;;;;;;;;;;;;;;;;;;;;;;;;
;;;   undefine.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;  Author: Philippe Bourdin
;;;  Date:   22-Jul-2013
;;;  Description:
;;;    Frees up to sixteen variables at once and sets them as 'undefined'.
;;;  Usage:
;;;    undefine, var1, var2, var3[, ...]

pro undefine, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16
	; shrink a variable first, then undefine it
	p1 = 0
	p2 = temporary (p1)
	p3 = temporary (p2)
	p4 = temporary (p3)
	p5 = temporary (p4)
	p6 = temporary (p5)
	p7 = temporary (p6)
	p8 = temporary (p7)
	p9 = temporary (p8)
	p10 = temporary (p9)
	p11 = temporary (p10)
	p12 = temporary (p11)
	p13 = temporary (p12)
	p14 = temporary (p13)
	p15 = temporary (p14)
	p16 = temporary (p15)
	dummy = temporary (p16)
end
