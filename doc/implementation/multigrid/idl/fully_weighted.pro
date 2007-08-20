;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   fully_weighted.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ucalgary.ca)
;;;  Date:   20-Aug-2007
;;;
;;;  Description:
;;;    Print spectral transfer function of fully-wighted coarse-graining
;;;    filter

kappa = linspace(0,!pi)
syms = texsyms()

save_state

!p.charsize=1.4

plot, kappa, (1+cos(kappa)) / 2., $
      XSTYLE=1, $
      XTITLE='!8k!D !N'+syms.delta+'!8x!X', $
      YTITLE='!6Spectral transfer function!X', $
      YMARGIN=[4,4], XMARGIN=[8,6]

ypos = 1.05

k_Ny  = !pi
k_Ny2 = !pi/2
scaley= 0.05
scalex = scaley*3
arrowx = [0., 0., -0.15, 0., 0.15] * scalex
arrowy = [1., 0.,  0.4, 0., 0.4] * scaley

chsize = !p.charsize*1
;
plots, k_Ny + arrowx, 1.01 + arrowy
xyouts, k_Ny-0.2, 1.1, $
        '!8k!B!6Ny,fine!N!X', $
        CHARSIZE=chsize
;
plots, k_Ny2 + arrowx, 1.01 + arrowy
xyouts, k_Ny2-0.17, 1.1, $
        '!8k!B!6Ny,coarse!N!X', $
        CHARSIZE=chsize

restore_state       

end
; End of file fully_weighted.pro
