;
;  $Id: yderxder_6th_ghost.pro,v 1.1 2008-01-18 10:29:39 ajohan Exp $
;
;  Second derivative d2f/dydx
;  - 6th-order
;  - with ghost cells
;
;***********************************************************************
function yderxder,f
  COMPILE_OPT IDL2,HIDDEN
;
  return, xderyder(f)
;
end
