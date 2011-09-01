;$Id$
FUNCTION pc_round,x,d
;
;  rounds off a number to the accuracy of d
;
;  compute sign of x
;
if x gt 0. then xsgn=1. else xsgn=-1.
;
;  compute default
;
xabs=abs(x)
if xabs gt 1. then begin
  d0=xsgn*10.^(fix(alog10(xabs)))
endif else begin
  d0=xsgn*10.^(fix(alog10(.1*xabs)))
endelse
print,'d0=',d0
;
;  set default
;
default,d,d0
return,fix(x/d+.5)*d
END
