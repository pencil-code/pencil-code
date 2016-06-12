; $Id$
;
;  reads the zaver.dat file
;  need to supply nprocy by hand...
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
bmxy=fltarr(nx,ny*nprocy,nnamexy)
;
close,1
openr,1,datatopdir+'/zaverages.dat'
;
fo='(8e12.4)'
default,w,.1
while not eof(1) do begin
  readf,1,t
  readf,1,bmxy,fo=fo
  ;oplot,bymxy,li=1
  print,t,max(bmxy)
  wait,w
end
close,1
;
END
