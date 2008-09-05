;;
;;  $Id$
;;
;;  Simple program that evaluates whether a Pencil Code variable is periodic.
;;
;;  Usage: pc_is_periodic, var
;;
;;  Here var is a three-dimensional physical variable with ghost zones.
;;
pro pc_is_periodic, f

  sizef=size(f)

  mx=sizef[1]
  my=sizef[2]
  mz=sizef[3]

  l1=3
  l2=mx-1-3
  m1=3
  m2=my-1-3
  n1=3
  n2=mz-1-3
  
  lperiodicx_bot=0
  lperiodicx_top=0
  lperiodicy_bot=0
  lperiodicy_top=0
  lperiodicz_bot=0
  lperiodicz_top=0

  if (min(f(l2-2:l2,m1:m2,n1:n2) eq f(0:l1-1, m1:m2       ,n1:n2    )) eq 1) $
      then lperiodicx_bot=1
  if (min(f(l1:l1+2,m1:m2,n1:n2) eq f(l2+1:l2+3,m1:m2     ,n1:n2    )) eq 1) $
      then lperiodicx_top=1
  if (min(f(l1:l2,m2-2:m2,n1:n2) eq f(l1:l2  ,0:m1-1      ,n1:n2    )) eq 1) $
     then lperiodicy_bot=1
  if (min(f(l1:l2,m1:m1+2,n1:n2) eq f(l1:l2  ,m2+1:m2+3   ,n1:n2    )) eq 1) $
     then lperiodicy_top=1
  if (min(f(l1:l2,m1:m2,n2-2:n2) eq f(l1:l2  ,m1:m2       ,0:n1-1   )) eq 1) $
     then lperiodicz_bot=1
  if (min(f(l1:l2,m1:m2,n1:n1+2) eq f(l1:l2  ,m1:m2       ,n2+1:n2+3)) eq 1) $
     then lperiodicz_top=1

  print, 'lperiodicx=', lperiodicx_bot, lperiodicx_top
  print, 'lperiodicy=', lperiodicy_bot, lperiodicy_top
  print, 'lperiodicz=', lperiodicz_bot, lperiodicz_top

end
