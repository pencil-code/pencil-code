; $Id: rslice_xy.pro,v 1.1 2002-08-17 10:37:49 brandenb Exp $
;
;  reads xy slices
;  this routine is not very general yet and needs to be adjusted
;  before it can be general purpose.
;
file_slice='tmp/proc0/ux.xy'
file_slice='tmp/proc0/uz.xy'
file_slice='tmp/proc0/divu.xy'
file_slice='tmp/proc0/bx.xy'
;
t=0.
xy_slice=fltarr(nx,ny)
;
close,1
openr,1,file_slice,/f77
while not eof(1) do begin
  readu,1,xy_slice,t
  print,t,min(xy_slice),max(xy_slice)
  ;plot,xy_slice(11,*),yr=[-1,1]*1e-3
  tvscl,xy_slice
  wait,.1
end
;
close,1
END
