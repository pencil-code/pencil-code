; $Id: rxyaver.pro,v 1.6 2003-10-16 17:18:39 brandenb Exp $
;
;  reads the xyaver.dat file
;
t=0.
bxmz=fltarr(nz*nprocz)
bymz=fltarr(nz*nprocz)
nprocz=2
;
close,1
openr,1,datatopdir+'/xyaverages.dat'
;
it=0
fo='(8e10.3)'
default,w,.1
while not eof(1) do begin
  readf,1,t & print,t
  readf,1,bxmz,bymz,fo=fo
  ;
  if it eq 0 then begin
    bxmzt=bxmz
    bymzt=bymz
    tt=t
  endif else begin
    bxmzt=[bxmzt,bxmz]
    bymzt=[bymzt,bymz]
    tt=[tt,t]
  endelse
  ;plot,bxmz
  ;oplot,bymz,li=1
  ;wait,w
  it=it+1
end
close,1
;
nt=n_elements(tt)
bxmzt=reform(bxmzt,nz*nprocz,nt)
bymzt=reform(bymzt,nz*nprocz,nt)
END
