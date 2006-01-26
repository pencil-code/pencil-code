pro rvid_plane,field,mpeg=mpeg,png=png,tmin=tmin,tmax=tmax,max=amax,$
               min=amin,extension=extension,nrepeat=nrepeat,wait=wait,$
               njump=njump,datadir=datadir,OLDFILE=OLDFILE,test=test,$
               proc=proc,ix=ix,iy=iy,ps=ps,iplane=iplane,imgdir=imgdir,$
               global_scaling=global_scaling,shell=shell,r_int=r_int,$
               r_ext=r_ext,zoom=zoom,colmpeg=colmpeg,exponential=exponential, $
               contourplot=contourplot,color=color,sqroot=sqroot
;
; $Id: rvid_plane.pro,v 1.19 2006-01-26 12:24:39 ajohan Exp $
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
default,iplane,0
default,wait,.03
default,r_int,0.5
default,r_ext,1.0
default,zoom,1.0
default,dimfile,'dim.dat'
default,varfile,'var.dat'
default,imgdir,'.'
default,color,1
default,pixelsize,1
default,ximg,1
default,yimg,1
;
;  Read the dimensions and precision (single or double) from dim.dat
;
mx=0L & my=0L & mz=0L & nvar=0L & prec=''
nghostx=0L & nghosty=0L & nghostz=0L
nprocx=0L & nprocy=0L & nprocz=0L
;
close,1
openr,1,datadir+'/'+'dim.dat'
readf,1,mx,my,mz,nvar
readf,1,prec
readf,1,nghostx,nghosty,nghostz
readf,1,nprocx,nprocy,nprocz
close,1
;
;  Set reasonable extension for 2-D runs.
;
if ( (mx ne 7) and (my ne 7) and (mz eq 7) ) then extension='xy'
if ( (mx ne 7) and (my eq 7) and (mz ne 7) ) then extension='xz'
if ( (mx eq 7) and (my ne 7) and (mz ne 7) ) then extension='yz'
;
if n_elements(proc) ne 0 then begin
  file_slice=datadir+'/proc'+str(proc)+'/slice_'+field+'.'+extension
endif else begin
  file_slice=datadir+'/slice_'+field+'.'+extension
endelse
;
;  double precision?
;
if prec eq 'D' then unit=1d0 else unit=1e0
;
nx=mx-2*nghostx
ny=my-2*nghosty
nz=mz-2*nghostz
ncpus = nprocx*nprocy*nprocz
;
if keyword_set(shell) then begin
  ;
  ; to mask outside shell, need full grid;  read from varfiles, as in rall.pro
  ;
  datalocdir=datadir+'/proc0'
  mxloc=0L & myloc=0L & mzloc=0L
  ;
  close,1
  openr,1,datalocdir+'/'+dimfile
  readf,1,mxloc,myloc,mzloc
  close,1
  ;
  nxloc=mxloc-2*nghostx
  nyloc=myloc-2*nghosty
  nzloc=mzloc-2*nghostz
  ;
  one=unit & zero=0.0*unit
  x=fltarr(mx)*one & y=fltarr(my)*one & z=fltarr(mz)*one
  xloc=fltarr(mxloc)*one & yloc=fltarr(myloc)*one & zloc=fltarr(mzloc)*one
  readstring=''
  ;
  for i=0,ncpus-1 do begin        ; read data from individual files
    if n_elements(proc) ne 0 then begin
      datalocdir=datadir+'/proc'+str(proc)
    endif else begin
      datalocdir=datadir+'/proc'+strtrim(i,2)
    endelse
    ; read processor position
    dummy=''
    ipx=0L &ipy=0L &ipz=0L
    close,1
    openr,1,datalocdir+'/'+dimfile
    readf,1, dummy
    readf,1, dummy
    readf,1, dummy
    readf,1, ipx,ipy,ipz
    close,1
    openr,1, datalocdir+'/'+varfile, /F77
    if (execute('readu,1'+readstring) ne 1) then $
          message, 'Error reading: ' + 'readu,1'+readstring
    readu,1, t, xloc, yloc, zloc
    close,1
    ;
    ;  Don't overwrite ghost zones of processor to the left (and
    ;  accordingly in y and z direction makes a difference on the
    ;  diagonals)
    ;
    if (ipx eq 0) then begin
      i0x=ipx*nxloc & i1x=i0x+mxloc-1
      i0xloc=0 & i1xloc=mxloc-1
    endif else begin
      i0x=ipx*nxloc+nghostx & i1x=i0x+mxloc-1-nghostx
      i0xloc=nghostx & i1xloc=mxloc-1
    endelse
    ;
    if (ipy eq 0) then begin
      i0y=ipy*nyloc & i1y=i0y+myloc-1
      i0yloc=0 & i1yloc=myloc-1
    endif else begin
      i0y=ipy*nyloc+nghosty & i1y=i0y+myloc-1-nghosty
      i0yloc=nghosty & i1yloc=myloc-1
    endelse
    ;
    if (ipz eq 0) then begin
      i0z=ipz*nzloc & i1z=i0z+mzloc-1
      i0zloc=0 & i1zloc=mzloc-1
    endif else begin
      i0z=ipz*nzloc+nghostz & i1z=i0z+mzloc-1-nghostz
      i0zloc=nghostz & i1zloc=mzloc-1
    endelse
    ;
    x[i0x:i1x] = xloc[i0xloc:i1xloc]
    y[i0y:i1y] = yloc[i0yloc:i1yloc]
    z[i0z:i1z] = zloc[i0zloc:i1zloc]
    ;
  endfor
  ; 
  xx = spread(x, [1,2], [my,mz])
  yy = spread(y, [0,2], [mx,mz])
  zz = spread(z, [0,1], [mx,my])
  rr = sqrt(xx^2+yy^2+zz^2)
  
  ; assume slices are all central for now -- perhaps generalize later
  ; nb: need pass these into boxbotex_scl for use after scaling of image;
  ;     otherwise pixelisation can be severe...
  ; nb: at present using the same z-value for both horizontal slices;
  ;     hardwired into boxbotex_scl, also.
  ix=mx/2 & iy=my/2 & iz=mz/2 & iz2=iz
  if extension eq 'xy' then rrxy =rr(nghostx:mx-nghostx-1,nghosty:my-nghosty-1,iz)
  if extension eq 'Xy' then rrxy2=rr(nghostx:mx-nghostx-1,nghosty:my-nghosty-1,iz2)
  if extension eq 'xz' then rrxz =rr(nghostx:mx-nghostx-1,iy,nghostz:mz-nghostz-1)
  if extension eq 'yz' then rryz =rr(ix,nghosty:my-nghosty-1,nghostz:mz-nghostz-1)
  ;
endif
;
t=0.*unit & islice=0
;
if (extension eq 'xy') then plane=fltarr(nx,ny)*unit
if (extension eq 'Xy') then plane=fltarr(nx,ny)*unit
if (extension eq 'xz') then plane=fltarr(nx,nz)*unit
if (extension eq 'yz') then plane=fltarr(ny,nz)*unit
size_plane=size(plane)
print, 'Array size: ', size_plane[0:size_plane[0]]
;
slice_xpos=0.*unit
slice_ypos=0.*unit
slice_zpos=0.*unit
slice_z2pos=0.*unit
;
;  open MPEG file, if keyword is set
;
dev='x' ;(default)
if keyword_set(png) then begin
  Nwx=zoom*size_plane[1] & Nwy=zoom*size_plane[2]
;  resolution=[!d.x_size,!d.y_size] ; set window size
  resolution=[Nwx,Nwy] ; set window size
  print, 'z-buffer resolution in pixels '+ $
      '(set with zoom=', strtrim(zoom,2), ') =', strtrim(resolution,2)
  set_plot, 'z'                   ; switch to Z buffer
  device, SET_RESOLUTION=resolution ; set window size
  itpng=0 ;(image counter)
  dev='z'
end else if keyword_set(mpeg) then begin
  ;Nwx=400 & Nwy=320
  ;Nwx=!d.x_size & Nwy=!d.y_size
  Nwx=zoom*size_plane[1] & Nwy=zoom*size_plane[2]
  resolution=[Nwx,Nwy] ; set window size
  print,'z-buffer resolution (in pixels)=',resolution
  set_plot, 'z'                   ; switch to Z buffer
  device, SET_RESOLUTION=resolution ; set window size
  dev='z'
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
if keyword_set(global_scaling) then begin
  first=1L
  close,1 & openr,1,file_slice,/f77
  while not eof(1) do begin
    if keyword_set(OLDFILE) then begin ; For files without position
      readu,1,plane,t
    endif else begin
      readu,1,plane,t,slice_z2pos
    endelse
    if keyword_set(exponential) then begin
      if (first) then begin
        amax=exp(max(plane))
        amin=exp(min(plane))
        first=0L
      endif else begin
        amax=max([amax,exp(max(plane))])
        amin=min([amin,exp(min(plane))])
      endelse
    endif else if keyword_set(sqroot) then begin
      if (first) then begin
        amax=sqrt(max(plane))
        amin=sqrt(min(plane))
        first=0L
      endif else begin
        amax=max([amax,sqrt(max(plane))])
        amin=min([amin,sqrt(min(plane))])
      endelse
    endif else begin
      if (first) then begin
        amax=max(plane)
        amin=min(plane)
        first=0L
      endif else begin
        amax=max([amax,max(plane)])
        amin=min([amin,min(plane)])
      endelse
    endelse
  end
  close,1
  print,'Scale using global min, max: ', amin, amax
endif
;
close,1 & openr,1,file_slice,/f77
ifirst=1
while not eof(1) do begin
  if keyword_set(OLDFILE) then begin ; For files without position
    readu,1,plane,t
  end else begin
    readu,1,plane,t,slice_z2pos
  end
;
; rescale data with optional parameter zoom
;
  if keyword_set(exponential) then begin
    plane2=rebinbox(exp(plane),zoom)
  endif else if keyword_set(sqroot) then begin
    plane2=rebinbox(sqrt(plane),zoom)
  endif else begin
    plane2=rebinbox(plane,zoom)
  endelse
;
; do masking, if shell set
;
  if keyword_set(shell) then begin
    white=255
    if extension eq 'xy' then begin
      zrr = rebinbox(reform(rrxy,nx,ny),zoom)
      indxy=where(zrr lt r_int or zrr gt r_ext)
      plane2(indxy)=white
    endif
    if extension eq 'Xy' then begin
      zrr2 = rebinbox(reform(rrxy2,nx,ny),zoom)
      indxy2=where(zrr2 lt r_int or zrr2 gt r_ext)
      plane2(indxy2)=white
    endif
    if extension eq 'xz' then begin
      yrr = rebinbox(reform(rrxz,nx,nz),zoom,/zdir)
      indxz=where(yrr lt r_int or yrr gt r_ext)
      plane2(indxz)=white
    endif
    if extension eq 'yz' then begin
      xrr = rebinbox(reform(rryz,ny,nz),zoom,/zdir)
      indyz=where(xrr lt r_int or xrr gt r_ext)
      plane2(indyz)=white
    endif
  endif
;
  if keyword_set(test) then begin
    print,t,min([plane2,xy,xz,yz]),max([plane2,xy,xz,yz])
  end else begin
    if t ge tmin and t le tmax then begin
      if ijump eq njump then begin
        ;if iy ne -1 then plot,plane2(*,iy),yr=[amin,amax],ps=ps
        ;if ix ne -1 then plot,plane2(ix,*),yr=[amin,amax],ps=ps
;
;  show image scaled between amin and amax and filling whole screen
;
        if keyword_set(contourplot) then begin
          contourfill,plane2,lev=grange(amin,amax,60)
        endif else begin
          tv, bytscl(plane2,min=amin,max=amax), iplane
        endelse
        ;tv,congrid(bytscl(plane2,min=amin,max=amax),!d.x_size,!d.y_size)
        xyouts, 0.05, 0.9, /normal, $
            '!8t!6='+string(t,fo="(f6.1)"), color=color, size=0.5*zoom
        if keyword_set(png) then begin
          istr2 = strtrim(string(itpng,'(I20.4)'),2) ;(only up to 9999 frames)
          image = tvrd()
;
;  make background white, and write png file
;
          ;bad=where(image eq 0) & image(bad)=255
          tvlct, red, green, blue, /GET
          imgname = imgdir+'/img_'+istr2+'.png'
          write_png, imgname, image, red, green, blue
          itpng=itpng+1 ;(counter)
          ;
        end else if keyword_set(mpeg) then begin
;
;  write directly mpeg file
;  for idl_5.5 and later this requires the mpeg license
;
          image = tvrd(true=1)
          if keyword_set(colmpeg) then begin
; ngrs seem to need to work explictly with 24-bit color to get 
; color mpegs to come out on my local machines...
            image24 = bytarr(3,Nwx,Nwy)
            tvlct, red, green, blue, /GET
          endif
          for irepeat=0,nrepeat do begin
            if keyword_set(colmpeg) then begin
              image24[0,*,*]=red(image[0,*,*])
              image24[1,*,*]=green(image[0,*,*])
              image24[2,*,*]=blue(image[0,*,*])
              mpeg_put, mpegID, image=image24, FRAME=itmpeg, /ORDER
            endif else begin
              mpeg_put, mpegID, window=2, FRAME=itmpeg, /ORDER
            endelse
            itmpeg=itmpeg+1 ;(counter)
          end
          print,islice,itmpeg,t,min([plane2]),max([plane2])
        end else begin
;
; default: output on the screen
;
          if (ifirst) then $
              print, '----islice--------t----------min------------max--------'
          print,islice,t,min([plane2]),max([plane2])
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
  ifirst=0
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
