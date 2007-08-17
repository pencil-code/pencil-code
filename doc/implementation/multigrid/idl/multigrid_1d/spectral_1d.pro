;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   spectral_1d.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   16-Nov-2005
;;;
;;;  Description:
;;;    Solve stationary heat conduction porblem
;;;      T'' = -q(x) , T(0)=T(1)=0
;;;    using spectral method.
;;;    This is (possibly up to boundary conditions) the same problem as in
;;;    multigrid.pro, so we can compare

function q, x
  ;; Heating profile
  w = 0.03
  x0 = 2./3.
  return, exp(-(x-x0)^2/(2.*w^2))/sqrt(2*!pi)/w
end
; ---------------------------------------------------------------------- ;

Nx = 256L
Lx = 1.
x = linspace([0,Lx],Nx,/PERIODIC)
dx = x[1]-x[0]
dk = 2*!pi/Lx
qtilde = FFT(q(x),-1)

;; For plotting spectra:
k = (findgen(Nx)-Nx/2)*dk
spect = shift(qtilde,Nx/2)

;; For spectral method:
k = shift(k,-Nx/2)              ; 0, dk, 2dk, .., (Nx/2-1)dk, -Nx/2 dk, .., -dk
spect = qtilde

k_2 = 1./(k^2 > 1.e-20)
k_2[0] = 0.

temp_tilde = -qtilde*k_2
temp = FFT(temp_tilde,1)

;plot, x, temp-temp[0], XRANGE=[0.6,0.8]

plot, x, temp
oplot, x, -0.00210580 + 0.0037*abs(sin(!pi*(x-2./3.))), COLOR=120

end
; End of file spectral_1d.pro
