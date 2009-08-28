;;;;;;;;;;;;;;;;;;;;;;;
;;;   calcvar.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   03-Aug-2002
;;;  $Id$
;;;
;;;  Description:
;;;   Calculate derived variables that are often needed, but should
;;;   not be calculated by default (to save memory).
;;;   Assumes that some variables and parameters exist.

var_ = ''
print, 'Which quantity do you want? '
print, '  rr, rb, jj, cs2, pp, Temp'
read, var_, PROMPT='>  '


case var_ of

    'rr'  : rr   = sqrt(xx^2+yy^2+zz^2)
    'bb'  : bb   = curl(aa)
    'jj'  : jj   = curl2(aa)
    'cs2' : cs2  = cs20*exp(gamma*ss+gamma_m1*(lnrho-lnrho0))
    'pp'  : pp   = cs20*rho0/gamma*exp(gamma*(ss+lnrho-lnrho0))
    'Temp': Temp = cs20/gamma_m1*exp(gamma*ss+gamma_m1*(lnrho-lnrho0))

endcase

end
; End of file calcvar.pro
