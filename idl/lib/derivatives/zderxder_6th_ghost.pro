;
;  $Id: zderxder_6th_ghost.pro,v 1.1 2008-01-18 10:29:39 ajohan Exp $
;
;  Second derivative d2f/dzdx
;  - 6th-order
;  - with ghost cells
;
;***********************************************************************
function zderxder,f
  COMPILE_OPT IDL2,HIDDEN
;
  return, xderzder(f)
;
end
