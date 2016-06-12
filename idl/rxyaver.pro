; $Id$
;
;  reads xyaver.dat file
;  read first global array sizes
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
;  read mesh
;
pc_read_grid,obj=obj
z=obj.z
;
;  usuable array boundaries
;
n1=3 & n2=mz-4 & n12=n1+indgen(nz)
zzz=z(n1:n2)
;
;  set input array
;
fmz=fltarr(nz*nprocz,nnamez)
;
;  read data
;
close,1
openr,1,datatopdir+'/xyaverages.dat'
print,'read from file: ',datatopdir+'/xyaverages.dat'
;
it=0
fo='(8e12.5)'
default,w,.01
while not eof(1) do begin
  readf,1,t & print,t
  readf,1,fmz,fo=fo
  if it eq 0 then begin
    fmzt=fmz
    tt=t
  endif else begin
    fmzt=[fmzt,fmz]
    tt=[tt,t]
  endelse
  it=it+1
end
ending:
close,1
;
nt=n_elements(tt)
fmzt=reform(fmzt,nz,nt,nnamez)
;
END
