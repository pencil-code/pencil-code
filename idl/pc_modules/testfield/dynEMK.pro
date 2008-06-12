
; Determination of the tensor Lambda = ( lam11  lam12 )  
;			               ( lam21  lam22 )
; in dE_0/dt = \Lambda E_0
; from the observed (measured) quantities defined below

;  11-jun-08/MR: coded

sigma = 0.151		; Growth rate of E_0
omega = 0.792		; oscillation frequency of E_0
phi   = 2.08		; phase shift of E_0,y vs. E_0,x
q     = 2.0		; ratio |E_0,y|/|E_0,x|

D = omega/tan(phi)

lam11 = sigma - D
lam22 = 2*D +lam11
lam12 = sqrt(sigma^2+omega^2)/q
lam21 = - (omega^2 + D^2)/lam12 

print, "lam11, lam22, lam12, lam21=", lam11, lam22, lam12, lam21

end
