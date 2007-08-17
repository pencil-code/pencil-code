;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   restrict_test.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   01-Mar-2007
;;;
;;;  Description:
;;;    Test our different restriction (coarse-graining) schemes.
;;;    See restrict.pro for a discussion.

; ---------------------------------------------------------------------- ;
function f, x
;;
;; The function to restrict
;;
;  return, 1*sin(!pi*x) + 0.2*cos(16*!pi*x)
;  return, 1*sin(!pi*x) + 0.2*cos(16*!pi*x)*sin(!pi*x)
;  return, (x-0.4)^3
;  return, cos(15*!pi*x)
;  return, tanh(10.*(x-0.5))
;  return, exp(-100.*(x-0.5)^2)
;
;  x1 = (x-0.5)+1.e-6
;  k = 2*!pi*48
;  return, sin(k*x1)/(k*x1)

  return, x^2*(x lt 0.5) + (0.375-0.5*x^2)*(x ge 0.5)

end

; ---------------------------------------------------------------------- ;

@restrict

N_coarse = 9
;N_fine   = 2*N_coarse - 1       ; Classical multigrid (vertex-centred,
                                ; i.e. every other point is on both grids)
N_fine   = 2*N_coarse          ; Cell-centred, coarse-grid points are not
                                ; on fine grid

red  = 120./255.*!d.table_size
blue =  40./255.*!d.table_size

x_coarse = linspace([0, 1.D], N_coarse)
x_fine   = linspace([0, 1.D], N_fine)

f_coarse = f(x_coarse)
f_fine   = f(x_fine)

f_rest1 = restrict1(f_fine, x_fine, x_coarse)
f_rest2 = restrict2(f_fine, x_fine, x_coarse)

x_ = linspace([0.,1.],200)
f_ = f(x_)
yr = minmax([f_,f_fine,f_rest1,f_rest2]) + [-1,1]*0.05

plot,  x_,   f_, $
    XRANGE=stretchrange(minmax(x_fine), 0.03), $
    YSTYLE=1, YRANGE=yr

; plot,  x_,   f_, $
;     XRANGE=[0.,1]+[-1,1]*0.03, $
;     YSTYLE=1, YRANGE=[0,1.1]+[-1,1]*0.15

oplot, x_fine,   f_fine,   PSYM=1
oplot, x_coarse, f_coarse, PSYM=-5, LINE=2

oplot, x_coarse, f_rest1, PSYM=2, COLOR=red
oplot, x_coarse, f_rest2, PSYM=4, COLOR=blue



if (running_gdl_p()) then begin
  plots, [0.12,0.2], [0.9, 0.9 ], /NORMAL
  plots, [0.12,0.2], [0.85,0.85], /NORMAL, PSYM=1
  plots, [0.12,0.2], [0.8 ,0.8 ], /NORMAL, PSYM=-5, LINESTYLE=2
  plots, [0.12,0.2], [0.75,0.75], /NORMAL, PSYM=2, COLOR=red
  plots, [0.12,0.2], [0.7 ,0.7 ], /NORMAL, PSYM=4, COLOR=blue 
  ;
  xyouts, 0.22, 0.9 , /NORMAL, 'f',        CHARSIZE=1.4
  xyouts, 0.22, 0.85, /NORMAL, 'f_fine',   CHARSIZE=1.4
  xyouts, 0.22, 0.8 , /NORMAL, 'f_coarse', CHARSIZE=1.4
  xyouts, 0.22, 0.75, /NORMAL, 'f_rest1',  CHARSIZE=1.4
  xyouts, 0.22, 0.7 , /NORMAL, 'f_rest2',  CHARSIZE=1.4
endif else begin
  esrg_legend, /BOX, SPOS='tl', $
                [ 'f', 'f_fine', 'f_rest1', 'f_rest2', 'f_coarse' ], $
      COLOR=    [   0,        0,       red,      blue,          0 ], $
      LINESTYLE=[   0,        0,         0,         0,          2 ], $
      PSYM=     [   0,        1,         2,         4,         -5 ]
endelse


print, $
    '      x_coarse     f_coarse     f_rest1      f_rest2     delta1     delta2'
print, $
    '     ---------------------------------------------------------------------'
print, $
    transpose([ $
                  [x_coarse                                      ], $
                  [f_coarse                                      ], $
                  [restrict1(f_fine, x_fine, x_coarse)           ], $
                  [restrict2(f_fine, x_fine, x_coarse)           ], $
                  [restrict1(f_fine, x_fine, x_coarse) - f_coarse], $
                  [restrict2(f_fine, x_fine, x_coarse) - f_coarse]  $
                  $
              ])


end
; End of file restrict.pro
