;
; Plot rhs of Poisson equation vs. numerical Laplacian of solution
; (should be the same)
;

;
; Usage:
;   @toto1
;

.r start
.r rall

.compile xder2_2nd
.compile yder2_2nd
.compile zder2_2nd

rho = exp(lnrho)
rho1 = del2(potself)

rr = sqrt(xx^2+yy^2+zz^2)
psym = 4
yr = max(rho1[l1+3:l2-3,m1+3:m2-3,n1+3:n2-3]) * [0.3,1.1]*1
yr = [0.99,1.01]
plot_binned, rr, rho1, $
             xr=[-0.3 ,0.3], $
             YRANGE=yr, $
             PSYM=psym
plot_binned, -rr, rho1, /OVER, PSYM=psym

;x_ = linspace([-0.2,2],1000)
;oplot, x_, exp(1.*exp(-x_^2/0.05^2)), COLOR=120

; Same thing , but manually
d2x = (shift(potself,1,0,0) + shift(potself,-1,0,0) - 2*potself)/dx^2
d2y = (shift(potself,0,1,0) + shift(potself,0,-1,0) - 2*potself)/dy^2
d2z = (shift(potself,0,0,1) + shift(potself,0,0,-1) - 2*potself)/dz^2

rho2 = d2x + d2y + d2z

psym = 1
plot_binned,  rr, rho2, /OVER, PSYM=psym
plot_binned, -rr, rho2, /OVER, PSYM=psym

plot_binned,  rr, rho, /OVER, PSYM=4, COLOR=120
plot_binned, -rr, rho, /OVER, PSYM=4, COLOR=120

;c = (1+cos(5*rho))
;plot_binned,  rr, rho + c*potself, /OVER, PSYM=2, COLOR=120
;plot_binned, -rr, rho + c*potself, /OVER, PSYM=2, COLOR=120

esrg_legend, SPOS='tr', /BOX, $
             ['del2(pot)', 'd2x+d2y_d2z(pot)', 'exp(-r^2/0.05^2)', 'rho'], $
             PSYM=[psym,4,0,4], $
             COLOR=[0,0,120,120]
             
