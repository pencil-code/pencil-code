;;;;;;;;;;;;;;;;;;;;
;;;   evol.pro   ;;;
;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   12-Nov-2001
;;;
;;;  Description:
;;;   Time evolution (from n.dat).

default, up_to_date, 0
default, nfile, datatopdir + '/n.dat'

if (not up_to_date) then begin
  data = input_table(nfile)
endif

up_to_date = 1

tt = data[1,*]
;
urms = data[2,*]
umax = data[3,*]
;
divurms = data[10,*]
divumax = data[11,*]
;
rrms = data[8,*]
rmax = data[9,*]

save_state

!p.multi=[0,2,2]

;; Velocity
yr = minmax([umax,urms])
plot, tt, urms, $
    YRANGE=yr, YSTYLE=3, $
    TITLE='!6Velocity!X', $
    XTITLE='!8t!X', YTITLE='!8u!X'
oplot, tt, umax, LINE=2
esrg_legend, $
    ['!6rms!X', '!6max!X'], $
    LINE=[0,2], SPOS='tl', /BOX


;; div u
yr = minmax([divumax,divurms])
plot, tt, divurms, $
    YRANGE=yr, YSTYLE=3, $
    TITLE='!6div !17u!X', $
    XTITLE='!8t!X', YTITLE='!17div !8u!X'
oplot, tt, divumax, LINE=2

;; ln rho
yr = minmax([rmax,rrms])
plot, tt, rrms, $
    YRANGE=yr, YSTYLE=3, $
    TITLE='!6log-density!X', $
    XTITLE='!8t!X', YTITLE='!6ln !7r!X'
oplot, tt, rmax, LINE=2

restore_state

end
; End of file evol.pro
