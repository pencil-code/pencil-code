;
; $Id$
;
;  read zprofile file
;
function pc_read_zprof, varname , datadir=datadir, debug=debug
;
IF (not keyword_set(datadir)) THEN datadir=pc_get_datadir()
;
;  read expected dimension and processor number
;
pc_read_dim,nzgrid=nzgrid,nprocz=nprocz,nprocy=nprocy,nprocx=nprocx,datadir=datadir
nz=nzgrid/nprocz
if keyword_set(debug) then print,'nzgrid,nprocz,nprocy,nprocx=',nzgrid,nprocz,nprocy,nprocx
;
;  define full arrays
;
zfull=fltarr(nzgrid)
afull=fltarr(nzgrid)
;
;  loop over all processors and records in file
;
izcount=0
for iprocz=0,nprocz-1 do begin
  procname='/proc'+str(iprocz*nprocy*nprocx)
  filename=datadir+procname+'/zprof_'+varname+'.dat'
  openr,1,filename
    for iz=0,nz-1 do begin
      readf,1,z,a
      zfull(izcount)=z
      afull(izcount)=a
      izcount=izcount+1
    endfor
  close,1
endfor
return,afull
;
end
