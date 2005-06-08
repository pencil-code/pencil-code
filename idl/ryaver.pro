; $Id: ryaver.pro,v 1.1 2005-06-08 04:30:03 brandenb Exp $
;
;  reads the yaver.dat file
;  need to supply nprocy by hand...
;
t=0.
@data/index
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
bmxz=fltarr(nx,nz*nprocz,nnamexz)
;
close,1
openr,1,datatopdir+'/yaverages.dat'
;
fo='(8e12.4)'
default,w,.1
while not eof(1) do begin
  readf,1,t
  readf,1,bmxz,fo=fo
  ;oplot,bymxz,li=1
  print,t,max(bmxz)
  wait,w
end
close,1
;
END
