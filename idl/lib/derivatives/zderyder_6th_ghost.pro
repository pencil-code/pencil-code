;
;  $Id: zderyder_6th_ghost.pro,v 1.1 2008-01-18 10:29:39 ajohan Exp $
;
;  Second derivative d2f/dzdy
;  - 6th-order
;  - with ghost cells
;
;***********************************************************************
function zderyder,f
  COMPILE_OPT IDL2,HIDDEN
;
  return, yderzder(f)
;
end
