;$Id: pdardel_strong.pro,v 1.11 2021/10/31 04:29:15 brandenb Exp $
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
!p.charsize=1.5
!x.margin=[9.2,.5]
!y.margin=[3.2,.5]
;
siz=1.8
siz2=1.2
sym = texsyms()
!x.title='# of procs'
!y.title=sym.mu+'!6s/step/point'
yr2=[5e-2,1e-0]
yr3=[1e-2,1e-0]
yr=[1e-5,3e-1]
xr=[5e1,7e4]
;
plot_oo,xr,yr,/nodata
xx=[1e1,1.8e4] & oplot,xx,1.02/xx,li=1
xx=[1e1,4.8e4] & oplot,xx,0.70/xx,li=2
;
circ_sym,1.0,1
file='dardel_256.txt'
a=rtable(file,head=1,2)
oplot,a(0,*),a(1,*),ps=8,col=55
oplot,[1,1]*512,yr2,col=55
xyouts,150,.001,'!6256!u3!n',siz=siz,col=55
xyouts,400,.015,'!6512',siz=siz2,col=55
;
file='dardel_512_gnu.txt'
a=rtable(file,head=1,2)
oplot,a(0,*),a(1,*),ps=8,col=155
oplot,[1,1]*2048,yr2,col=155
xyouts,500,.0005,'!6512!u3!n',siz=siz,col=155
xyouts,1300,.015,'!62048',siz=siz2,col=155
;
loadct,6
circ_sym,1.3,0
file='dardel_1024_gnu.txt'
a=rtable(file,head=1,2)
oplot,a(0,*),a(1,*),ps=8,col=122
oplot,[1,1]*4096,yr2,col=122
xyouts,1300,.0001,'!61024!u3!n',siz=siz,col=122
xyouts,2600,.015,'!64096',siz=siz2,col=122
loadct,5
;
circ_sym,1.3,1
file='dardel_2048.txt'
a=rtable(file,head=1,2)
oplot,a(0,*),a(1,*),ps=8,col=122
oplot,[1,1]*8192,yr2,col=122;,li=2
xyouts,4000,.00003,'!62048!u3!n',siz=siz,col=122
xyouts,5500,.015,'!616384',siz=siz2,col=122
;
circ_sym,1.3,1
file='dardel_4096.txt'
a=rtable(file,head=1,2)
oplot,a(0,*),a(1,*),ps=8
oplot,[1,1]*16384,yr2
xyouts,9000,.000016,'!64096!u3!n',siz=siz
xyouts,13000,.015,'!632384',siz=siz2
;
print,"$mv idl.ps ../fig/pdardel_strong.ps"
print,"$mv idl.ps ../figs/pdardel_strong.ps"
END
