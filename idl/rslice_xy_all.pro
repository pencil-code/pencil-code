pro rslice_xy_all,file,plane
;
; $Id: rslice_xy_all.pro,v 1.3 2002-10-02 20:11:14 dobler Exp $
;
; This program reads video snapshots from all the processors
; in the xy or xz plane.
;
;
dummy=''
;datatopdir='data'
close,1
openr,1,datatopdir+'/'+'dim.dat'
readf,1,dummy
readf,1,dummy
readf,1,nghostx,nghosty,nghostz
readf,1,nprocx,nprocy,nprocz
close,1
;
; Assuming same size on every processor
;
datadir=datatopdir+'/proc0'
close,1
openr,1,datadir+'/dim.dat'
readf,1, nnx,nny,nnz,nna
readf,1, dummy
readf,1, dummy
readf,1, ipx,ipy,ipz
close,1
;
nx=nnx-nghostx*2
ny=nny-nghosty*2
nz=nnz-nghostz*2
;
; Setting the right plane
;
if (plane EQ 'xy') then begin
  grid=ny
  nprocgrid=nprocy
  deltaproc=1
  print,'Showing the xy-plane'
endif else begin
  grid=nz
  nprocgrid=nprocz
  deltaproc=nprocy
  print,'Showing the xz-plane'
endelse 

t=0.
loc_slice=fltarr(nx,grid)
slice_glob=fltarr(nx,grid*nprocgrid)
;
for i=1,nprocgrid do begin
  j=fix((i-1)*deltaproc)
  close,i
  print,'Reading data/proc'+str(j)+'/'+file
  openr,i,datatopdir+'/proc'+str(j)+'/'+file,/f77
end
;
while not eof(1) do begin
  for i=1,nprocgrid do begin
    ;
    start=(i-1)*grid
    stop =i*grid-1
    readu,i,loc_slice,t
      slice_glob(*,start:stop)=loc_slice
  end
  print,t,min(slice_glob),max(slice_glob)
  tvscl,slice_glob
;  surface,slice_glob
  wait,.1
end
;
for i=1,nprocgrid do begin
  close,i
end
;
END


