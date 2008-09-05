; $Id$
;
;  reads xz slices
;  this routine is not very general yet and needs to be adjusted
;  before it can be general purpose.
;
file_slice=datatopdir+'/proc0/uz.xz'
file_slice=datatopdir+'/proc0/ux.xz'
file_slice=datatopdir+'/proc0/lnrho.xz'
file_slice=datatopdir+'/proc0/bx.xz'
;
t=0.
xz_slice=fltarr(nx,nz)
slice_ypos=0.
;
close,1
openr,1,file_slice,/f77
while not eof(1) do begin
  readu,1,xz_slice,t,slice_ypos
  print,t
  ;plot,xz_slice(11,*),yr=[-1,1]*1e-3
  tvscl,xz_slice
  wait,.1
end
;
close,1
END
