 function integral,x,y,xrange=xrange,accumulate=accumulate
;+
; ROUTINE:        integral
;
; USEAGE:         result=integral(x,y,xrange=xrange,accumulate=accumulate)
;
; PURPOSE:        integration by trapezoid rule
;
; INPUT:
;   x             vector of x axis points. If XRANGE is not set,
;                 limits of integration will be from x(0) to x(n_elements(x)-1)
;                 If XRANGE is not set, X need not be monitonic.
;
;   y             vector of corresponding y axis points
;
; KEYWORD_INPUT:
;
;   xrange        2 element vector, lower and upper limit of integration
;                 The use of XRANGE with a non-monitonic variation in x may
;                 produce incorrect results.
;
;   accumulate    if set, return value is a vector giving the accumulated
;                 integral at each value of x. XRANGE can't be used when
;                 this option is set.
;
; OUTPUT:         result of integration
;
; EXAMPLE:
;;                 /4
;; find            |  x dx
;;                 /0 
;  x=findgen(5) & y=x
;  print,integral(x,y)     
;  8.00000                             ; answer for linear integrand is exact
;
;
;;                 /5.5
;; find            |  x^2 dx
;;                 /0 
;
;  x=findgen(7) & y=x^2               &$
;  print,"    numerical     analytic" &$
;  print,integral(x,y,xrange=[0,5.5]), 5.5^3/3
;     56.3750      55.4583             
;
;; use more resolution in x to improve answer
;
;  x=findgen(551)/100. & y=x^2        &$
;  print,"    numerical     analytic" &$
;  print,integral(x,y), 5.5^3/3
;     55.4584      55.4583             ; much better
;
;  author:  Paul Ricchiazzi                            3NOV92
;           Institute for Computational Earth System Science
;           University of California, Santa Barbara
;
; REVISIONS:
; sep93: fixed error in treatment of xrange, added examples
; apr96: allow use of xrange with monitonically decreasing x
;-

if n_params() eq 0 then begin
  xhelp,'integral'
  return,0
endif  

nn=n_elements(x)
if nn ne n_elements(y) then message,'x and y vectors must be same size'
dx=shift(x,-1)-x
yy=.5*(shift(y,-1)+y)

if keyword_set(accumulate) then begin
  sum=fltarr(nn)
  for i=1,nn-1 do sum(i)=sum(i-1)+dx(i-1)*yy(i-1)
  return,sum
endif    
      
if n_elements(xrange) eq 0 then return,total(dx(0:nn-2)*yy(0:nn-2))

; rest of code is to treat end points when xrange is set...

;
xmin=min(xrange)
xmax=max(xrange)

ii=where(x ge xmin and x le xmax,nc)

if nc eq nn then return,total(dx(0:nn-2)*yy(0:nn-2))

n1=ii(0)

if n1 eq -1 then return,0.

n2=ii(nc-1)
sum=0.

; sum up points fully inside range

if n2 gt n1 then sum=total(dx(n1:n2-1)*yy(n1:n2-1))

; now add in rest of area inside integration limits


if n1 gt 0 then begin
  if x(n1) lt x(n2) then x1=xmin else x1=xmax
  y1=y(n1)+(x1-x(n1))*(y(n1+1)-y(n1))/(x(n1+1)-x(n1))
  sum1=.5*(y1+y(n1))*(x(n1)-x1)
  sum=sum+sum1
endif
if n2 lt nn-1 then begin
  if x(n1) lt x(n2) then x2=xmax else x2=xmin
  y2=y(n2)+(x2-x(n2))*(y(n2+1)-y(n2))/(x(n2+1)-x(n2))
  sum2=.5*(y2+y(n2))*(x2-x(n2))
  sum=sum+sum2
endif

return,sum

end
