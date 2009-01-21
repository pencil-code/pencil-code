;$Id$
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=4 & !p.thick=4 & !x.thick=4 & !y.thick=4 
end
;
;  plot plot total (kinetic and thermal energy)
;  obtain plot range from local parameters.pro file 
;
default,yr1,[0,7.2]
default,yr2,[6.87,6.88]
@parameters
;
;  mv idl.ps
;
;  read time-series file
;
pc_read_ts,obj=obj
tt=obj.t
ekin=obj.ekin
eth=obj.eth
;
!p.charsize=1.6
plot,tt,eth+ekin,yr=yr1,xtit='!8t!6',ytit='!6energy'
oplot,tt,eth,li=1
oplot,tt,ekin,li=2
;
;  inset with finer range
;
!p.position=[.35,.35,.85,.75]
plot,tt,eth+ekin,yr=yr2,/noerase
!p.position=0
;
END
