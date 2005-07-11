; $Id: pc_read_zprof.pro,v 1.1 2005-07-11 12:33:12 brandenb Exp $
;
;  read zprofile file
;
function pc_read_zprof, varname , datadir=datadir, debug=debug
;
default, datadir, 'data'
;
;  read expected dimension and processor number
;
pc_read_dim,nzgrid=nzgrid,nprocz=nprocz
nz=nzgrid/nprocz
if keyword_set(debug) then print,'nzgrid,nprocz=',nzgrid,nprocz
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
  procname='/proc'+str(iprocz)
  filename=datadir+procname+'/zprof_'+varname+'.dat'
  openr,1,filename
    for iz=0,nz-1 do begin
      readf,1,z,a
      zfull(izcount)=z
      afull(izcount)=a
      izcount=izcount+1
    endfor
  close,1
  return,afull
endfor
;
end
