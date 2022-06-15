function sine_step,x,x0,width,shift=shift

default, shift, 0

if (width eq 0) then width = 1e-37

xi = (x-x0)/width - shift
xi = xi > (-1.0)
xi = xi < 1.0
;;print,minmax(xi)
return, 0.5*(1+sin(0.5*!dpi*xi))

END














