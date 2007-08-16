;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   convergence.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   29-May-2007
;;;
;;;  Description:
;;;    Test convergence

;; These should all have the same amount of data
;; [And we assume that at least one dimension is >= 64]
c64 = input_table('32.ascii') 
c32 = input_table('32.ascii') 
c16 = input_table('16.ascii') 
c8 = input_table('8.ascii') 
c4 = input_table('4.ascii') 
c2 = input_table('2.ascii') 

N = (size(c64))[2]

idx = linspace(1,N,N)

plot, idx, /NODATA, $
      XRANGE=[0,150], $
      YRANGE=[1.e-18,1200], YSTYLE=1, $
      /YLOG, $
      XTITLE='!6Number !8N!6 of iterations!X', $
      YTITLE='!3||!6Residual!3||!X'

col = linspace(80, 160, 6)

oplot, idx, c64[2,*], COLOR=col[0]
oplot, idx, c32[2,*], COLOR=col[1]
oplot, idx, c16[2,*], COLOR=col[2]
oplot, idx, c8[2,*],  COLOR=col[3]
oplot, idx, c4[2,*],  COLOR=col[4]
oplot, idx, c2[2,*],  COLOR=col[5]

oplot, idx, 1.e2*0.48^idx, LINE=2
oplot, idx, 1.e2*0.81^idx, LINE=2
oplot, idx, 1.e2*0.92^idx, LINE=2

expos = '~ ' + ['0.48', '0.81', '0.92'] + '!UN!N'

esrg_legend, SPOS='tr', /BOX, $
             '!6'+['64', '32', '16', '8', '4', '2', expos]+'!X', $
             COLOR=[col,0,0,0], $
             LINESTYLE=[0,0,0,0,0,0, 2,2,2]

end
; End of file convergence.pro
