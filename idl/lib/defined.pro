; $Id$
;
;  If argument is undefined, it is initialized to zero
;  15-jun-02/axel: coded
;
;  06-jun-2004/wolf: Added following remark
;    Does anybody use this? Beware:
;    undefined(x) will not always do what one expects: if x is defined,
;    but zero undefined will return 0, i.e. logical false. 
;

FUNCTION defined,var
  if n_elements(var) eq 0 then var=0 else var=var
  return,var
END
