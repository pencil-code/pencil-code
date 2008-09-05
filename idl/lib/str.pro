FUNCTION str,x,_extra=e
COMPILE_OPT IDL2,HIDDEN
;+
;  $Id$
;-
  return,strtrim(string(x,_extra=e),2)
END

