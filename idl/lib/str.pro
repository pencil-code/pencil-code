FUNCTION str,x,_extra=e
;+
;  $Id: str.pro,v 1.1 2002-07-22 22:59:40 brandenb Exp $
;-
  return,strtrim(string(x,_extra=e),2)
END

