;$Id: add3.tex,v 1.570 2019/10/02 21:06:15 brandenb Exp $
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
end
;
!p.charsize=1.7
!x.margin=[6.8,.5]
!y.margin=[3.2,.5]
!x.title='!6check-in rank'
!y.title='!6# of line changes'
;
file='line_additions.txt'
a=rtable(file,2)
add=a[0,*]
sub=a[1,*]
dif=add-sub
n=n_elements(add)
i=findgen(n)
xr=[0,42]
yr=[5e1,1e6]
plot_io,xr,yr,/nodata
circ_sym,1.3,1 & oplot,i,add,ps=-8,col=122
circ_sym,1.3,0 & oplot,i,sub,ps=-8,col=55 & loadct,6
circ_sym,1.3,0 & oplot,i,dif,ps=-8,col=122 & loadct,5
p1=linfit(i,alog(add))
p2=linfit(i,alog(sub))
p3=linfit(i,alog(dif))
oplot,xr,exp(xr*p1[1]+p1[0]),l=2,col=122
oplot,xr,exp(xr*p2[1]+p2[0]),l=2,col=55 & loadct,6
oplot,xr,exp(xr*p3[1]+p3[0]),l=2,col=122 & loadct,5
oplot,xr*0+34,yr,li=3
;
xx=2. & dx=3.
legend,xx,dx,10^3.2,0,col=122,'!6additions'
legend,xx,dx,10^2.7,0,col= 55,'!6removals' & loadct,6
legend,xx,dx,10^2.2,0,col=122,'!6difference' & loadct,5
;
print,"$mv idl.ps ~/public_html/pencil-code/newsletter/2020/2/fig/pline_additions.ps"
END
