; $Id: ptimings.pro,v 1.3 2003-05-13 04:54:53 brandenb Exp $
;
red = 122
blue = 55
organge = 166
fg = !p.color
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
  ;red=0 & blue=0 & fg=0
end
;
;  mv idl.ps ../figs/ptimings.ps
;
fact=1.
a=rtable('timings.dat',2,head=1)
b=rtable('horseshoe.dat',2,head=1)
c=rtable('kabul.dat',2,head=1)
d=rtable('horseshoe_mega.dat',2,head=1)
n1=reform(a(0,*)) & t1=reform(a(1,*))
n2=reform(b(0,*)) & t2=reform(b(1,*))
n3=reform(c(0,*)) & t3=reform(c(1,*))
n4=reform(d(0,*)) & t4=reform(d(1,*))
;
save_state
sym = texsyms()
!p.multi=0
!p.charsize=1.2
!x.title='# of procs'
!y.title=sym.mu+'!6s/step/point'
!x.range=[.8,160]
!y.range=fact*[.07,15]
!x.style=3
!y.style=3

plot_oo, n1, fact*t1, PSYM=-1, LINE=fg
oplot,   n2, fact*t2, PSYM=-5,li=0, COLOR=red
oplot,   n3, fact*t3, PSYM=-6,li=2, COLOR=blue
oplot,   n4, fact*t4, PSYM=-7,li=3, COLOR=organge
;
esrg_legend, ['!6Origin3000!X', '!6Horseshoe!X', '!6KIS cluster!X', 'GigaBits'], $
    LINE=[1,0,2,3], $
    COLOR=[fg,red,blue,organge], $
    PSYM=[-1,-5,-6,-7], $
    SPOS='tr', /BOX

restore_state
;
; xx=[1,120] & oplot,xx,4./xx^.7
print,'import ptimings.jpg'
print,'scp2 ptimings.jpg $scr/ccp2001'

end
