; $Id: rxyaver.pro,v 1.9 2005-06-08 12:57:20 brandenb Exp $
;
;  read global sizes
;
default, datatopdir, 'data'
default, dimfile, 'dim.dat'
;
mx=0L & my=0L & mz=0L & nvar=0L & naux=0L
prec=''
nghostx=0L & nghosty=0L & nghostz=0L
;
nprocx=0L & nprocy=0L & nprocz=0L
close,1
openr,1,datatopdir+'/'+dimfile
readf,1,mx,my,mz,nvar,naux
readf,1,prec
readf,1,nghostx,nghosty,nghostz
readf,1,nprocx,nprocy,nprocz
close,1
;
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
;
;  reads the xyaver.dat file
;
t=0.
;nprocz=1
bxmz=fltarr(nz*nprocz)
bymz=fltarr(nz*nprocz)
alpijz=fltarr(nz*nprocz,3,3)
etaijkz=fltarr(nz*nprocz,3,3)
;
close,1
openr,1,datatopdir+'/xyaverages.dat'
print,'read from file: ',datatopdir+'/xyaverages.dat'
;
it=0
fo='(8e12.4)'
default,w,.1
while not eof(1) do begin
  readf,1,t & print,t
  readf,1,bxmz,bymz,alpijz,etaijkz,fo=fo
stop
  ;readf,1,bxmz,bymz,alpijz,fo=fo
  ;readf,1,bxmz,bymz,fo=fo
  ;
  if it eq 0 then begin
    bxmzt=bxmz
    bymzt=bymz
    alpijzt=alpijz
    etaijkzt=etaijkz
    tt=t
  endif else begin
    bxmzt=[bxmzt,bxmz]
    bymzt=[bymzt,bymz]
    alpijzt=[alpijzt,alpijz]
    etaijkzt=[etaijkzt,etaijkz]
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
alpijzt=reform(alpijzt,nz*nprocz,nt,3,3)
etaijkzt=reform(etaijkzt,nz*nprocz,nt,3,3)
END
