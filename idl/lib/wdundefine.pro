;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   wdundefine.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   15-Sep-2004
;;;  Description:
;;;    Get rid of one or several variables like DELVAR, but also works
;;;    in program scripts and subroutines (instead of silently
;;;    aborting like DELVAR does).
;;;    Based on undefine.pro by D. Fanning.

pro wdundefine, var1, var2, var3, var4, var5, var6, var7, var8

   if n_params() eq 0 then message, 'One argument required in call to UNDEFINE'

   if (n_elements(var1) gt 0) then tmp = size(temporary(var1))
   if (n_elements(var2) gt 0) then tmp = size(temporary(var2))
   if (n_elements(var3) gt 0) then tmp = size(temporary(var3))
   if (n_elements(var4) gt 0) then tmp = size(temporary(var4))
   if (n_elements(var5) gt 0) then tmp = size(temporary(var5))
   if (n_elements(var6) gt 0) then tmp = size(temporary(var6))
   if (n_elements(var7) gt 0) then tmp = size(temporary(var7))
   if (n_elements(var8) gt 0) then tmp = size(temporary(var8))
end
