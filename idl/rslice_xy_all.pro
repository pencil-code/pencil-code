pro rslice_xy_all,file
;
; $Id: rslice_xy_all.pro,v 1.1 2002-08-22 13:49:09 nilshau Exp $
;
; This program reads video snapshots from all the processors
; in the xy plane.
;
;
dummy=''
datatopdir='tmp'
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
datadir='tmp/proc0'
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
;
t=0.
xy_slice=fltarr(nx,ny)
slice_glob=fltarr(nx,ny*nprocy)
;
for i=1,nprocy do begin
  j=i-1
  close,i
  openr,i,'tmp/proc'+str(j)+'/'+file,/f77
end
;
while not eof(1) do begin
  for i=1,nprocy do begin
    ;
    ystart=(i-1)*ny
    ystop =i*ny-1
    readu,i,xy_slice,t
      slice_glob(*,ystart:ystop)=xy_slice
  end
  print,t,min(xy_slice),max(xy_slice)
  tvscl,slice_glob
;  surface,slice_glob
  wait,.1
end
;
for i=1,nprocy do begin
  close,i
end
;
END


