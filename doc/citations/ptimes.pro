;$Id: add3.tex,v 1.502 2015/02/19 06:35:19 brandenb Exp $
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
;  mv idl.ps fig/ptimes.ps
;  convert idl.ps ~/Figures/Pencil/2015/ptimes.png
;
!p.charsize=1.7
!x.margin=[6.1,.5]
!y.margin=[3.2,.5]
!x.title='!6year'
!y.title='!6number of papers'
;
a=rtable('times.txt',2)
n=a(0,*)
y=a(1,*)
print,n
print
print,total(n)
plot,y,n,ps=10,yr=[0,62]
END
