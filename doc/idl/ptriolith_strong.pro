;$Id: ptriolith_strong.pro,v 1.7 2021/10/31 07:26:26 brandenb Exp $
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=13,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
;  mv idl.ps ../fig/ptriolith_strong.ps
;  mv idl.ps ~/tex/prop/computing/boulder/janus/fig/ptriolith_strong.ps
;  mv idl.ps ~/tex/prop/computing/boulder/summit/fig/ptriolith_strong.ps
;
!p.charsize=1.5
!x.margin=[9.2,.5]
!y.margin=[3.2,.5]
;
siz=1.8
siz2=1.4
sym = texsyms()
!x.title='# of procs'
!y.title=sym.mu+'!6s/step/point'
yr2=[3e-2,1e-0]
yr3=[1e-2,1e-0]
yr=[3e-5,1e-0]
;
file='triolith.dat'
a=rtable(file,head=1,2)
plot_oo,a(0,*),a(1,*),ps=1,xr=[1,5e4],yr=yr,/nodata
xx=[1e1,1.8e4] & oplot,xx,1.02/xx,li=1
;
circ_sym,1.3,0
file='triolith_256.dat'
a=rtable(file,head=1,2)
oplot,a(0,*),a(1,*),ps=8,col=55
oplot,[1,1]*256,yr2,col=55
xyouts,70,.002,'!6256!u3!n',siz=siz,col=55
xyouts,200,.015,'!6256',siz=siz2,col=55
;
file='triolith_512.dat'
a=rtable(file,head=1,2)
oplot,a(0,*),a(1,*),ps=4,col=155
oplot,[1,1]*1024,yr2,col=155
xyouts,260,.0004,'!6512!u3!n',siz=siz,col=155
xyouts,650,.015,'!61024',siz=siz2,col=155
;
circ_sym,1.3,1
file='triolith_2304.dat'
a=rtable(file,head=1,2)
oplot,a(0,*),a(1,*),ps=8,col=122
oplot,[1,1]*9216,yr2,col=122,li=2
oplot,[1,1]*18432,yr3,col=122
xyouts,1000,.00008,'!62304!u3!n',siz=siz,col=122
xyouts,5000,.015,'!69216',siz=siz2,col=122
xyouts,10000,.005,'!618432',siz=siz2,col=122
;
loadct,6
circ_sym,1.8,1
a=rtable('summit.dat',2)
oplot,a[0,*],a[1,*],col=122,ps=-8
xyouts,24,.07,'!6Summit!c     576!u3!n',col=122
loadct,5
;
print,"$mv idl.ps ../figs/ptriolith_strong.eps"
END
