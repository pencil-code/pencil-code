; $Id: rslice.pro,v 1.2 2002-06-02 07:51:39 brandenb Exp $
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
  plot,xz_slice(11,*),yr=[-1,1]*1e-3
  tvscl,xz_slice
  wait,.1
end
;
close,1
END
