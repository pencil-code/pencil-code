FUNCTION grange,x1,x2,n,log=log
;
;  generates a uniform sequence of values in the range [x1,x2]
;  $Id$
;
if keyword_set(log) then begin
  lx1=alog(x1)
  lx2=alog(x2)
  return,exp(findgen(n)/(n-1)*(lx2-lx1)+lx1)
endif else begin
  return,findgen(n)/(n-1)*(x2-x1)+x1
endelse
;
END
