; $Id: rslice.pro,v 1.4 2002-08-17 08:40:37 brandenb Exp $
;
;  reads xz slices
;  this routine is not very general yet and needs to be adjusted
;  before it can be general purpose.
;
file_slice='tmp/proc0/uz.xz'
;
t=0.
xz_slice=fltarr(nx,nz)
;
close,1
openr,1,file_slice,/f77
while not eof(1) do begin
  readu,1,xz_slice,t
  print,t
  ;plot,xz_slice(11,*),yr=[-1,1]*1e-3
  tvscl,xz_slice
  wait,.1
end
;
close,1
END
