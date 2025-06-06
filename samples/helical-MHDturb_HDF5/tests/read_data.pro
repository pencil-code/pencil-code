;;;
;;; Read time series and data cube, write a few values
;;;

datadir = '../data'

pc_read_ts, datadir=datadir, obj=ts

pc_read_var, datadir=datadir, obj=var, /trimall

fo="(a,5g15.7)"

openw, 1, 'read_data.out'
printf, 1, 'ts.times : ', ts.t(0:4)
printf, 1, 'aa(5,5,0:4,1) : ', reform(var.aa(5,5,0:4,1)), fo=fo
close, 1
;
;  exit, to prevent idl from waiting for further input
;
exit
