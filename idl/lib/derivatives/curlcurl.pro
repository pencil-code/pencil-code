;
; $Id: curlcurl.pro,v 1.1 2008-01-18 10:33:19 ajohan Exp $
;
; Calculate two consecutive curls of a 3-D vector field.
;
; 18-jan-08/anders: coded
;
function curlcurl, f
  COMPILE_OPT IDL2, HIDDEN
;
  w=make_array(size=size(f),/nozero)
  w=graddiv(f)-del2(f)
;
  return, w
;
end
