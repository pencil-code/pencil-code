;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   stability_classical.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (wdobler [at] gmail [dot] com)
;;;  Date:   18-Aug-2008
;;;
;;;  Description:
;;;    Plot stability polynomials for classical Runge-Kutta methods

L = 4
x = linspace(0,L)
yrange = [-1,1]*1.05

loadct, 5
blue = 50
orange = 150
red = 120
orange = 140
yellow = 185
dashed = 2

if (!d.name eq 'PS') then begin
  DEVICE, XSIZE=18, YSIZE=9, $
      XOFFSET=3
  device, /COLOR
endif
aspect_ratio = aspect_pos(minmax(yrange, /RANGE) / L, $
                          MARGIN=[0.1, 0.03, 0.15, 0.03])

plot, x, exp(-x), $
      YRANGE=yrange, XSTYLE=1, YSTYLE=3, $
      POS=aspect_ratio, $
      XTITLE='!3Cou!X', $
      YTITLE='!8A!3(Cou)!X'
ophline, [-1, 0, 1]

oplot, COLOR=blue  , x, 1 - x
oplot, COLOR=red   , x, 1 - x + x^2/2.
oplot, COLOR=orange, x, 1 - x + x^2/2. - x^3/6.
oplot, COLOR=yellow, x, 1 - x + x^2/2. - x^3/6. + x^4/24.

oplot, COLOR=blue  , x, 1.0 / (1+x), LINE=dashed

esrg_legend, SPOS='br', /BOX, $
             ['Exact', $
              'Expl. Euler', 'Expl. 2nd-order', 'Expl.3rd-order', 'Expl.4th-order', $
              'Impl. Euler'], $
             COLOR=[0, $
                    blue, red, orange, yellow, $
                    blue], $
             LINESTYLE=[0, $
                        0, 0, 0, 0, $
                        dashed]
                        


end
; End of file stability_classical.pro
