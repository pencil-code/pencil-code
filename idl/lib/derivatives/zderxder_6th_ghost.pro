;;
;;  $Id: zderxder_6th_ghost.pro,v 1.2 2008-06-10 13:07:41 ajohan Exp $
;;
;;  Second derivative d2f/dzdx
;;  - 6th-order
;;  - with ghost cells
;;
function zderxder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  return, xderzder(f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
end
