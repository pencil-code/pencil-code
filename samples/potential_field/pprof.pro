;$Id$
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=4 & !p.thick=4 & !x.thick=4 & !y.thick=4
end
;
;  compare vertical current profiles
;
;  mv idl.ps ~/tex/mhd/alpha_pot/fig/pprof.ps
;
!p.multi=0
!p.charsize=1.6
!x.title='!8z!6'
!y.title='!8J!dx!n!6'
restore,'prof_c1.sav'
plot,zzz,jjjz,xr=[0.,.6],yr=[-.01,.13],ps=-4
;
restore,'prof_pot.sav'
oplot,zzz,jjjz,col=122,li=0
;
restore,'prof_pwd.sav'
oplot,zzz,jjjz,col=188,li=0
;
oplot,.5*[1,1],[-.1,.2]
oplot,x,x-x
end
