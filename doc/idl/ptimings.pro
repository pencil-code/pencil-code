; $Id: ptimings.pro,v 1.2 2002-10-14 18:31:26 dobler Exp $
;
red = 122
blue = 55
fg = !p.color
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
  red=0 & blue=0 & fg=0
end
;
;  mv idl.ps ../figs/ptimings.ps
;
fact=1.
a=rtable('timings.dat',2,head=1)
b=rtable('horseshoe.dat',2,head=1)
c=rtable('kabul.dat',2,head=1)
n1=reform(a(0,*)) & t1=reform(a(1,*))
n2=reform(b(0,*)) & t2=reform(b(1,*))
n3=reform(c(0,*)) & t3=reform(c(1,*))
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

plot_oo, n1, fact*t1, PSYM=-1, LINE=1
oplot,   n2, fact*t2, PSYM=-5,li=0, COLOR=122
oplot,   n3, fact*t3, PSYM=-6,li=2, COLOR=55
;
esrg_legend, ['!6Origin3000!X', '!6Horseshoe!X', '!6KIS cluster!X'], $
    LINE=[1,0,2], $
    COLOR=[fg,red,blue], $
    PSYM=[-1,-5,-6], $
    SPOS='tr', /BOX

restore_state
;
; xx=[1,120] & oplot,xx,4./xx^.7
print,'import ptimings.jpg'
print,'scp2 ptimings.jpg $scr/ccp2001'

end
