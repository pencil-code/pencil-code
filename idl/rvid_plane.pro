pro rvid_plane,field,mpeg=mpeg,png=png,tmin=tmin,tmax=tmax,amax=amax,amin=amin,$
  nrepeat=nrepeat,wait=wait,njump=njump,datadir=datadir,OLDFILE=OLDFILE,$
  test=test,proc=proc,ix=ix,iy=iy,ps=ps
;
; $Id: rvid_plane.pro,v 1.5 2003-10-18 20:44:13 brandenb Exp $
;
;  reads and displays data in a plane (currently with tvscl)
;  and plots a curve as well (cross-section through iy)
;
;  if the keyword /mpeg is given, the file movie.mpg is written.
;  tmin is the time after which data are written
;  nrepeat is the number of repeated images (to slow down movie)
;
;  Typical calling sequence
;  rvid_plane,'uz',amin=-1e-1,amax=1e-1,/proc
;
default,proc,0
default,ix,-1
default,iy,-1
default,ps,0
default,extension,'xz'
default,amax,.05
default,amin,-amax
default,field,'lnrho'
default,datadir,'data'
default,nrepeat,0
default,njump,0
default,tmin,0.
default,tmax,1e38
default,wait,.03
;
if keyword_set(proc) then begin
  file_slice=datadir+'/proc0/slice_'+field+'.'+extension
endif else begin
  file_slice=datadir+'/slice_'+field+'.'+extension
endelse
print,!d.y_size
;
;  Read the dimensions and precision (single or double) from dim.dat
;
mx=0L & my=0L & mz=0L & nvar=0L & prec=''
nghostx=0L & nghosty=0L & nghostz=0L
;
close,1
openr,1,datadir+'/'+'dim.dat'
readf,1,mx,my,mz,nvar
readf,1,prec
readf,1,nghostx,nghosty,nghostz
close,1
;
;  double precision?
;
if prec eq 'D' then unit=1d0 else unit=1e0
;
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
;
t=0.*unit & islice=0
;
if extension eq 'xy' then plane=fltarr(nx,ny)*unit
if extension eq 'xz' then plane=fltarr(nx,nz)*unit
help,plane
slice_xpos=0.*unit
slice_ypos=0.*unit
slice_zpos=0.*unit
slice_z2pos=0.*unit
;
;  open MPEG file, if keyword is set
;
dev='x' ;(default)
if keyword_set(png) then begin
  resolution=[!d.x_size,!d.y_size] ; set window size
  print,'z-buffer resolution (in pixels)=',resolution
  set_plot, 'z'                   ; switch to Z buffer
  device, SET_RESOLUTION=resolution ; set window size
  itpng=0 ;(image counter)
  dev='z'
end else if keyword_set(mpeg) then begin
  ;Nwx=400 & Nwy=320
  Nwx=!d.x_size & Nwy=!d.y_size
  if (!d.name eq 'X') then window,2,xs=Nwx,ys=Nwy
  mpeg_name = 'movie.mpg'
  print,'write mpeg movie: ',mpeg_name
  mpegID = mpeg_open([Nwx,Nwy],FILENAME=mpeg_name)
  itmpeg=0 ;(image counter)
end
;
;  allow for jumping over njump time slices
;  initialize counter
;
ijump=njump ;(make sure the first one is written)
;
close,1 & openr,1,file_slice,/f77
while not eof(1) do begin
if keyword_set(OLDFILE) then begin ; For files without position
  readu,1,plane,t
end else begin
  readu,1,plane,t,slice_z2pos
end
;
if keyword_set(test) then begin
  print,t,min([plane,xy,xz,yz]),max([plane,xy,xz,yz])
end else begin
  if t ge tmin and t le tmax then begin
    if ijump eq njump then begin
      ;if iy ne -1 then plot,plane(*,iy),yr=[amin,amax],ps=ps
      ;if ix ne -1 then plot,plane(ix,*),yr=[amin,amax],ps=ps
      ;
      ;  show image scaled between amin and amax and filling whole screen
      ;
      tv,congrid(bytscl(plane,min=amin,max=amax),!d.x_size,!d.y_size)
      xyouts,.93,1.13,'!8t!6='+string(t,fo="(f6.1)"),col=1,siz=2
      if keyword_set(png) then begin
        istr2 = strtrim(string(itpng,'(I20.4)'),2) ;(only up to 9999 frames)
        image = tvrd()
        ;
        ;  make background white, and write png file
        ;
        ;bad=where(image eq 0) & image(bad)=255
        tvlct, red, green, blue, /GET
        imgname = 'img_'+istr2+'.png'
        write_png, imgname, image, red, green, blue
        itpng=itpng+1 ;(counter)
        ;
      end else if keyword_set(mpeg) then begin
        ;
        ;  write directly mpeg file
        ;  for idl_5.5 and later this requires the mpeg license
        ;
        image = tvrd(true=1)
        for irepeat=0,nrepeat do begin
          mpeg_put, mpegID, window=2, FRAME=itmpeg, /ORDER
          itmpeg=itmpeg+1 ;(counter)
        end
        print,islice,itmpeg,t,min([plane]),max([plane])
      end else begin
        ;
        ; default: output on the screen
        ;
        print,islice,t,min([plane]),max([plane])
      end
      ijump=0
      wait,wait
      ;
      ; check whether file has been written
      ;
      if keyword_set(png) then spawn,'ls -l '+imgname
      ;
    end else begin
      ijump=ijump+1
    end
  end
  islice=islice+1
end
end
close,1
;
;  write & close mpeg file
;
if keyword_set(mpeg) then begin
  print,'Writing MPEG file..'
  mpeg_save, mpegID, FILENAME=mpeg_name
  mpeg_close, mpegID
end
;
END
