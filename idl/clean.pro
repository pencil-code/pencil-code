;;;;;;;;;;;;;;;;;;;;;
;;;   clean.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   11-Apr-2004
;;;
;;;  Description:
;;;    delete those variables that prevent `.r start' after a premature
;;;    `.r r' of `.r rall'.
;;;  Usage:
;;;    [.r r before .r start -> problems]
;;;    @clean
;;;    .r start
;;;    .r r
;;;  Discussion:
;;;    `.rnew start' does segfault (or something like `double free or
;;;    [whatever]') with most IDL versions I am using. No fun.

delvar, x,y,z, mx,my,mz

; End of file clean.pro
