pro rslice_xy_all,file,plane,mpeg=mpeg,itmax=itmax,fmax=fmax,fmin=fmin,nrepeat=nrepeat,OLDFILE=OLDFILE
;
; $Id$
;
; This program reads video snapshots from all the processors
; in the xy or xz plane, depending on whether plane='xy' or plane='xz'.
;
; if the keyword /mpeg is given, the file movie.mpg is written.
; itmax is the maximum number of time slices
; nrepeat is the number of repeated images (to slow down movie)
;
; if /OLDFILE is given, uses slice files wich donot contain the
; position variable
;
;
;
default,itmax,10
default,nrepeat,0
;
dummy=''
datatopdir='data'
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
slice_pos=0.
;
for i=1,nprocgrid do begin
  j=fix((i-1)*deltaproc)
  close,i
  print,'Reading data/proc'+str(j)+'/'+file
  openr,i,datatopdir+'/proc'+str(j)+'/'+file,/f77
end
;
;  open MPEG file, if keyword is set
;
if keyword_set(mpeg) then begin
  Nwx=512 & Nwy=512
  if (!d.name eq 'X') then window,2,xs=Nwx,ys=Nwy
  mpeg_name = 'movie.mpg'
  mpegID = mpeg_open([Nwx,Nwy],FILENAME=mpeg_name)
  itmpeg=0 ;(image counter)
end
;
it=0 ;(image counter)
while not eof(1) and it le itmax do begin
  for i=1,nprocgrid do begin
    ;
    start=(i-1)*grid
    stop =i*grid-1
if keyword_set(OLDFILE) then begin
    readu,i,loc_slice,t
end else begin
    readu,i,loc_slice,t,slice_pos
end
      slice_glob(*,start:stop)=loc_slice
  end
  ffmin=min(slice_glob)
  ffmax=max(slice_glob)
  default,fmin,ffmin
  default,fmax,ffmax
  image=bytscl(slice_glob,min=fmin,max=fmax)
  tv,image
  if keyword_set(mpeg) then begin
    image = tvrd(true=1)
    for irepeat=0,nrepeat do begin
      mpeg_put, mpegID, window=2, FRAME=itmpeg, /ORDER
      itmpeg=itmpeg+1 ;(counter)
    end
    print,itmpeg,t,ffmin,ffmax
  end else begin
    print,it,t,ffmin,ffmax
  end
  it=it+1 ;(counter)
  wait,.1
end
;
;  write & close mpeg file
;
if keyword_set(mpeg) then begin
  print,'Writing MPEG file..'
  mpeg_save, mpegID, FILENAME=mpeg_name
  mpeg_close, mpegID
end
;
;  close all files
;
for i=1,nprocgrid do close,i
;
END
