;;
;;  $Id: yderxder_6th_ghost.pro,v 1.2 2008-06-10 13:07:41 ajohan Exp $
;;
;;  Second derivative d2f/dydx
;;  - 6th-order
;;  - with ghost cells
;;
function yderxder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  return, xderyder(f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
end
