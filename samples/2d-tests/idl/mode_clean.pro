;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   mode_clean.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   13-Jan-2004
;;;
;;;  Description:
;;;   Remove modes kx=0,1,3,5,7,9,.. removing average and combining
;;;   f() and shifted f().

;; Remove zero mode
weight = 0.*lnrho
weight[l1:l2,m1:m2,n1+1:n2-1] = 1.0
weight[l1:l2,m1:m2,[n1,n2]  ] = 0.5
uu[*,*,*,0] = uu[*,*,*,0] $
              - total(weight*exp(lnrho)*uu[*,*,*,0]) $
                / total(weight*exp(lnrho))

;; Remove odd modes
uu[l1:l2,*,*,*]    = 0.5*(uu[l1:l2,*,*,*] +shift(uu[l1:l2,*,*,*] ,nx/2,0,0,0))
ss[l1:l2,*,*]      = 0.5*(ss[l1:l2,*,*]   +shift(ss[l1:l2,*,*]   ,nx/2,0,0))
lnrho[l1:l2,*,*,*] = 0.5*(lnrho[l1:l2,*,*]+shift(lnrho[l1:l2,*,*],nx/2,0,0))



end
; End of file mode_clean.pro
