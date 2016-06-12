; $Id$
;
;  reads the yaver.dat file, puts the result into fmxz array
;  this routine keeps only the last time
;
t=0.
@data/index
datatopdir = pc_get_datadir(datatopdir)
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
fmxz=fltarr(nx,nz,nnamexz)
;
close,1
openr,1,datatopdir+'/yaverages.dat'
;
fo='(8e12.4)'
default,w,.01
while not eof(1) do begin
  readf,1,t
  readf,1,fmxz,fo=fo
  print,t,max(fmxz)
  wait,w
end
close,1
;
END
