FUNCTION str,x,_extra=e
COMPILE_OPT IDL2,HIDDEN
;+
;  $Id: str.pro,v 1.2 2004-05-05 16:42:57 mee Exp $
;-
  return,strtrim(string(x,_extra=e),2)
END

