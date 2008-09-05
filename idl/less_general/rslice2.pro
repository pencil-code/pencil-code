; $Id$
;
;  show two different slices
;
file_slice1='tmp/proc0/uz.xz'
file_slice2='tmp/proc0/lnrho.xz'
;
t=0.
xz_slice1=fltarr(nx,nz)
xz_slice2=fltarr(nx,nz)
;
close,1 & openr,1,file_slice1,/f77
close,2 & openr,2,file_slice2,/f77
;
!p.multi=[0,1,2]
while not eof(1) do begin
  readu,1,xz_slice1,t
  readu,2,xz_slice2,t
  print,t
  plot,xz_slice1(11,*),yr=[-1,1]*1e+1
  plot_io,exp(xz_slice2(11,*));,yr=[-1,1]*1e+1
  tvscl,xz_slice1,0
  tvscl,xz_slice2,2
  wait,.1
end
;
!p.multi=0
close,1
close,2
END
