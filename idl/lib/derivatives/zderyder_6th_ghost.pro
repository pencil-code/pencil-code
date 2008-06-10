;;
;;  $Id: zderyder_6th_ghost.pro,v 1.2 2008-06-10 13:07:41 ajohan Exp $
;;
;;  Second derivative d2f/dzdy
;;  - 6th-order
;;  - with ghost cells
;;
function zderyder,f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t
  COMPILE_OPT IDL2,HIDDEN
;
  return, yderzder(f,ghost=ghost,bcx=bcx,bcy=bcy,bcz=bcz,param=param,t=t)
;
end
