t=0.
bxmz=fltarr(nz*nprocz)
bymz=fltarr(nz*nprocz)
;
close,1
openr,1,'tmp/zaverages.dat'
;
fo='(8e10.3)'
while not eof(1) do begin
  readf,1,t & print,t
  readf,1,bxmz,bymz,fo=fo
  plot,bxmz
  oplot,bymz,li=1
  wait,.1
end
close,1
;
END
