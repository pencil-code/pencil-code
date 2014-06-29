;
; DOCUMENT ME!!! What do these lines of code do, and why?
;
function rebinbox, a, zoom, zdir=zdir, sample=sample

sa=size(a)
nx0=sa(1) & ny0=sa(2)
if sa(0) eq 1 then ny0=1

if not keyword_set(zdir) then begin
  b=fltarr(nx0+1,ny0+1)
  b(0:nx0-1,0:ny0-1)=a
  b(nx0,0:ny0-1)=a(0,*)
  b(0:nx0-1,ny0)=a(*,0)
  b(nx0,ny0)=a(0,0)
  temp = rebin(b, (nx0+1)*zoom, (ny0+1)*zoom, sample=sample)
  return,temp(0:nx0*zoom-1,0:ny0*zoom-1)
endif else begin
  b=fltarr(nx0+1,ny0)
  b(0:nx0-1,*)=a
  b(nx0,*)=a(0,*)
  temp = rebin(b, (nx0+1)*zoom, ny0*zoom, sample=sample)
  return,temp(0:nx0*zoom-1,0:(ny0-1)*zoom)
endelse

end
