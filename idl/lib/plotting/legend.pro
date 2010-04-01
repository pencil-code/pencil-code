;$Id$
pro legend,x1,dx,y,lin,text,size=size,ps=ps,symsize=symsize,col=col,thick=thick,log=log
if n_params(0) eq 0 then begin
  print,'legend,x1,dx,y,lin,text,size=size'
  return
endif
;
;  plot legend
;
default,ps,0
default,symsize,1.6
if !d.name eq 'PS' then default,col,0 else default,col,255
;
;  the line starts at x1 and its length is dx
;
x=[x1,x1+dx/2.,x1+dx]
if keyword_set(log) then x(1)=sqrt(x(0)*x(2))
x2=x1+1.2*dx
oplot,x,[y,y,y],lin=lin,col=col,thick=thick
;
;  place symbol in the middle of the line
;
if ps ne 0 then oplot,[1,1]*(x1+dx/2.),[y,y],ps=ps,symsize=symsize,col=col
;
;  add legend symbols
;
if n_elements(size) eq 0 then xyouts,x2,y,text,size=1.8,col=col
if n_elements(size) ne 0 then xyouts,x2,y,text,size=size,col=col
;
end
