; $Id: defined.pro,v 1.1 2002-06-15 09:30:51 brandenb Exp $
;
;  If argument is undefined, it is initialized to zero
;  15-jun-02/axel: coded
;
FUNCTION defined,var
  if n_elements(var) eq 0 then var=0 else var=var
  return,var
END
