FUNCTION grange,x1,x2,n
;
;  generates a uniform sequence of values in the range [x1,x2]
;  $Id: grange.pro,v 1.1 2004-03-24 11:06:14 brandenb Exp $
;
return,findgen(n)/(n-1)*(x2-x1)+x1
END
