;  $Id$
;
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=14,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
;  plots the test figure that is checked in in this directory
;
!x.title='!8E!6!dth!n [!7l!6g/cm!u3!n (km/s)!u2!n Mm!u3!n]'
!y.title='!8T!6 [K]'
!p.charsize=1.8
plot,ts.eth,ts.ttm,xr=[0,1600],yr=[0,3.5e4],ps=-6
;
; restore,'../heating_noionize/1.sav' & oplot,Eth,TTm,col=122
END
