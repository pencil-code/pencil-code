;$Id: add3.tex,v 1.595 2023/04/26 08:25:00 brandenb Exp $
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
!p.charsize=1.7
!x.margin=[6.8,1.8]
!y.margin=[3.2,0.4]
;
urms=.1
kk=[1.,512.]
kNy=max(kk)
!x.title='!8k!6/!8k!6!dNy!n'
!y.title='!7m!8k!6!u2!8n!6-1!n/!8c!6!ds!n'
xr=[1e-3,1]
yr=[3e-11,3.]
xr=[1e-3,1]
yr=[3e-11,3.]
xr=[1e-1,1]
yr=[2e-2,2.]
plot_oo,kk/kNy,1e-4*kk/urms,xr=xr,yr=yr
oplot,kk/kNy,1e-9*kk^3/urms,col=122
oplot,kk/kNy,2e-14*kk^5/urms,col=55
oplot,[1.,1.]*.5,yr,li=3,col=155
loadct,6
oplot,xr,[1.,1.]*.2,li=3,col=122
loadct,5
;
kk=[1.,64.]
kNy=max(kk)
oplot,kk/kNy,2e-3*kk/urms,li=2
oplot,kk/kNy,5e-6*kk^3/urms,li=2,col=122
oplot,kk/kNy,2e-8*kk^5/urms,li=2,col=55
oplot,[1.,1.]*.28,yr,li=3,col=155
;
oplot,xr,.1/xr,li=1,col=155
print,"$mv idl.ps ../figs/phyper.eps"
;
END
