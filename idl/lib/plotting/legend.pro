pro legend,x1,dx,y,lin,text,size=size,ps=ps,symsize=symsize,col=col,thick=thick,log=log,charthick=charthick,vshift=vshift
if n_params(0) eq 0 then begin
  print,'legend,x1,dx,y,lin,text,size=size'
  return
endif
;
default,ps,0
default,vshift,0
default,symsize,1.6
if !d.name eq 'PS' then default,col,0 else default,col,255
;
x=[x1,x1+dx/2.,x1+dx]
if keyword_set(log) then x(1)=sqrt(x(0)*x(2))
x2=x1+1.2*dx
oplot,x,[y,y,y],lin=lin,ps=ps,symsize=symsize,col=col,thick=thick
if n_elements(size) eq 0 then xyouts,x2,y+vshift,text,size=1.8,col=col,charthick=charthick
if n_elements(size) ne 0 then xyouts,x2,y+vshift,text,size=size,col=col, charthick=charthick
;
end
