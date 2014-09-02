;
pro rvid_box_cyl,field,$
             bottom=bottom,no_orig_aspect=no_orig_aspect,ctable=ctable,$
             xrot=xrot,yrot=yrot,tunit=tunit,timetitle=timetitle,bar=bar,$
             xsize=xsize,ysize=ysize,xpos=xpos,ypos=ypos,zpos=zpos,$
             amax=amax,amin=amin,tmax=tmax,tmin=tmin,skip=skip,dev=dev,$
             norm=norm,itpng=itpng,png=png,imgprefix=imgprefix,imgdir=imgdir
;
common pc_precision, zero, one
;
default,dev,'x'
;
if keyword_set(dev) then begin
   dev1=dev
end else begin
   dev1=''
endelse
;
pc_read_grid,obj=grid,/trim,/q
pc_set_precision, dim=dim, /quiet
pc_read_dim, obj=dim, datadir=datadir, /quiet
nx=dim.nx & ny=dim.ny & nz=dim.nz
;
default,field,'rho'
default,bottom,0.0
default,bar,0.0
default,no_orig_aspect,0.0
default,png,0.0
default,xrot,30.0
default,zrot,30.0
default,ctable,5
default,tunit,2*!pi
default,timetitle,'T!d0!n'
default,xsize,600
default,ysize,600
default,xpos,0.2
default,ypos,0.5
default,zpos,0.6
default,amax,0.05
default,amin,-amax
default,tmin,0.0
default,tmax,1e38
default,skip,0
default,norm,1.
default,itpng,0
default,imgprefix,'img_'
default,imgdir,'.'
;
loadct,ctable
;
if (not keyword_set(datatopdir)) then datatopdir=pc_get_datadir()
datadir=datatopdir

file_slice1=datadir+'/slice_'+field+'.xy'
file_slice2=datadir+'/slice_'+field+'.xz'
file_slice3=datadir+'/slice_'+field+'.yz'
file_slice4=datadir+'/slice_'+field+'.xy2'

close, 1 & openr, 1, file_slice1, /f77, swap_endian=swap_endian
close, 2 & openr, 2, file_slice2, /f77, swap_endian=swap_endian
close, 3 & openr, 3, file_slice3, /f77, swap_endian=swap_endian
close, 4 & openr, 4, file_slice4, /f77, swap_endian=swap_endian
islice=0L
;
t=zero
tmp_xy2=fltarr(nx,ny)*one
tmp_xy=fltarr(nx,ny)*one
tmp_xz=fltarr(nx,nz)*one
tmp_yz=fltarr(ny,nz)*one
slice_xpos=zero
slice_ypos=zero
slice_zpos=zero
slice_z2pos=zero
;
rad=grid.x & phi=grid.y & zed=grid.z
nrad=n_elements(rad) 
nphi=n_elements(phi) 
nzed=n_elements(zed) 
;
tmp=(rad-max(rad)*cos(max(phi)))/max(rad)
scale=1.
rr=scale*tmp
;
;correct aspect ratio 
;
tmp=(zed-min(zed))/grid.Lz & zz=tmp-mean(tmp)
aspect_ratio = grid.Lz/grid.Lx
if (keyword_set(no_orig_aspect)) then aspect_ratio=1.
zz=aspect_ratio*zz
;
scale3,xrange=[0,1],yrange=[0,1],zrange=[0,1],$
       ax=(360+xrot) mod 360,az=zrot
;
if (!d.name eq 'X') then wdwset,xs=xsize,ys=ysize
;
while ( (not eof(1)) and (t le tmax) ) do begin
;
  if ( (t ge tmin) and (t le tmax) and (islice mod (skip+1) eq 0) ) then begin
    readu, 1, tmp_xy, t, slice_zpos
    readu, 2, tmp_xz, t, slice_ypos
    readu, 3, tmp_yz, t, slice_xpos
    readu, 4, tmp_xy2, t, slice_z2pos
;
    xy_int = round((tmp_xy-amin)/(amax-amin) * 255)
    xz_int = round((tmp_xz-amin)/(amax-amin) * 255)
    yz_int = round((tmp_yz-amin)/(amax-amin) * 255)
    xy2_int = round((tmp_xy2-amin)/(amax-amin) * 255)
;
;for iphi=0,nphi-2 do begin; nphi-2 do begin;nphi-2 do begin
;
;r-phi
;
    if dev1 ne 'z' then begin
       set_plot,'z'
       device, set_resolution=[xsize,ysize]
    end
;
    for iphi=0,nphi-2 do begin  ; nphi-2 do begin;nphi-2 do begin
;
;  3 _____2
;   |     |
;   |     |
;   |_____|
;   0     1
;
; xy
;
       x0=rr[0     ]*cos(phi[iphi])   & y0=rr[0     ]*sin(phi[iphi])
       x1=rr[nrad-1]*cos(phi[iphi])   & y1=rr[nrad-1]*sin(phi[iphi])
;   
       x2=rr[nrad-1]*cos(phi[iphi+1]) & y2=rr[nrad-1]*sin(phi[iphi+1])
       x3=rr[0     ]*cos(phi[iphi+1]) & y3=rr[0     ]*sin(phi[iphi+1])
;
       X = [x0,x1,x2,x3]           + xpos 
       Y = [y0,y1,y2,y3]           + ypos    
       z = replicate(zz[nzed/2],4) + zpos
;
       polyfill,x,y,z,pattern=xy_int[*,iphi],$
                image_coord=[[0,0],[nrad-1,0],[nrad-1,0],[0,0]],/t3d
;
; xy2
;
       if (keyword_set(bottom)) then begin
          z = replicate(0.,4) 
          polyfill,x,y,z,pattern=xy2_int[*,iphi],$
                   image_coord=[[0,0],[nrad-1,0],[nrad-1,0],[0,0]],/t3d
       endif
;
    endfor
;
; r-z
;
    for ized=0,nzed-2 do begin  ; nphi-2 do begin;nphi-2 do begin
;
;  3 _____2
;   |     |
;   |     |
;   |_____|
;   0     1
;
       x0=rr[0     ]*cos(phi[nphi/2])   & z0=zz[ized]
       x1=rr[nrad-1]*cos(phi[nphi/2])   & z1=zz[ized]
;   
       x2=rr[nrad-1]*cos(phi[nphi/2]) & z2=zz[ized+1]
       x3=rr[0     ]*cos(phi[nphi/2]) & z3=zz[ized+1]
;
       X = [x0,x1,x2,x3]   + xpos 
       Y = replicate(0.,4) + ypos    
       z = [z0,z1,z2,z3]   + zpos
;
       polyfill,x,y,z,pattern=xz_int[*,ized],$
                image_coord=[[0,0],[nrad-1,0],[nrad-1,0],[0,0]],/t3d
    endfor
;
    for iphi=0,nphi-2 do begin
;
;  3 _____2
;   |     |
;   |     |
;   |_____|
;   0     1
;
       x0=rr[nrad/2]*cos(phi[iphi])   & y0=rr[nrad/2]*sin(phi[iphi])
       x1=rr[nrad/2]*cos(phi[iphi+1]) & y1=rr[nrad/2]*sin(phi[iphi+1])
;                                                                 
       x2=rr[nrad/2]*cos(phi[iphi+1]) & y2=rr[nrad/2]*sin(phi[iphi+1])
       x3=rr[nrad/2]*cos(phi[iphi])   & y3=rr[nrad/2]*sin(phi[iphi])
;
       z0=zz[0]      & z1=zz[0]
       z2=zz[nzed-1] & z3=zz[nzed-1]
;
       X = [x0,x1,x2,x3] + xpos 
       Y = [y0,y1,y2,y3] + ypos    
       z = [z0,z1,z2,z3] + zpos
;
       polyfill,x,y,z,pattern=reform(yz_int[iphi,*]),$
                image_coord=[[0,0],[nzed-1,0],[nzed-1,0],[0,0]],/t3d
    endfor
;
; Time label
;
    time=t/tunit 
    time_int = fix(time)
    time_dec = fix(100*(time-time_int))
    if (time_dec ge 10) then begin
       time_str=strtrim(time_int,2)+'.'+strtrim(time_dec,2)
    endif else begin
       time_str=strtrim(time_int,2)+'.0'+strtrim(time_dec,2)
    endelse
;
    xyouts,.01,0.95,'!8t!x='+time_str+' '+timetitle,col=1,siz=1.
;
; Draw color bar.
;
    if (keyword_set(bar)) then begin
       default, bsize, 1.
       default, bformat, '(f5.2)'
       default, bnorm, 1.
       default, divbar, 2
       default, blabel, ''
       !p.title=blabel
;   colorbar, pos=[.89,.15,.90,.85], color=1, div=divbar,$
       colorbar_co, pos=[.1,.05,.9,.06], color=1, div=divbar,$
                 range=[amin,amax]/bnorm, $ ;,/right, /vertical, $
                 format=bformat, charsize=bsize, title=title
       !p.title=''
    endif
;
    a=tvrd() 
    bad=where(a eq 0) & a(bad)=255
    
    if (keyword_set(png)) then begin
       istr2 = strtrim(string(itpng,'(I20.4)'),2) ;(only up to 9999 frames)
       tvlct, red, green, blue, /GET
       imgname = imgprefix+istr2+'.png'
       write_png, imgdir+'/'+imgname, a, red, green, blue
       itpng=itpng+1 
    endif

    if dev1 ne 'z' then begin
       set_plot,'x'
       if (islice eq 1) then WINDOW, retain=2,XSIZE=xsize,YSIZE=ysize
       tv,a       ;,/true
    endif
;
 endif else begin       ; Read only time.                                                 
    dummy=zero
    readu, 1, xy, t
    readu, 2, dummy & readu, 3, dummy & readu, 4, dummy
 endelse
;
;  Ready for next video slice.
;
 islice=islice+1
;
; output
;
 if (islice eq 1) then $
    print, '   islice        t        min/norm     max/norm        amin     amax'
 print, islice, t, $
        min([min(tmp_xy2),min(tmp_xy),min(tmp_xz),min(tmp_yz)])/norm, $
        max([max(tmp_xy2),max(tmp_xy),max(tmp_xz),max(tmp_yz)])/norm, $
        amin, amax, format='(i9,e12.4,4f13.7)'
     
endwhile
;
end
