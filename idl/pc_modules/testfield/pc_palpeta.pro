;$Id$
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=4 & !p.thick=4 & !x.thick=4 & !y.thick=4
end
;
;  mv idl.ps ~/tex/notes/Poiseuille/fig/1024c_alpij.ps
;  mv idl.ps ~/tex/notes/Poiseuille/fig/256c_alpij.ps
;
restore,'alpetaij.sav'
;
!p.multi=[0,2,4]
!p.charsize=2.3
;
plot,zzz,alpijm(*,0,0),ytit='!7a!d!8xx!6!n'
oplot,zzz,zzz*0
;
plot,zzz,alpijm(*,1,1),ytit='!7a!d!8yy!6!n'
oplot,zzz,zzz*0
;
plot,zzz,alpijm(*,0,1),ytit='!7a!d!8xy!6!n'
oplot,zzz,zzz*0
;
plot,zzz,alpijm(*,1,0),ytit='!7a!d!8yx!6!n'
oplot,zzz,zzz*0
;
plot,zzz,etaijm(*,0,0),ytit='!7g!d!8xx!6!n'
oplot,zzz,zzz*0
;
plot,zzz,etaijm(*,1,1),ytit='!7g!d!8yy!6!n'
oplot,zzz,zzz*0
;
plot,zzz,etaijm(*,0,1),ytit='!7g!d!8xy!6!n'
oplot,zzz,zzz*0
;
plot,zzz,etaijm(*,1,0),ytit='!7g!d!8yx!6!n'
oplot,zzz,zzz*0
;
END
