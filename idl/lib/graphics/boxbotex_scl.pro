pro boxbotex_scl,inxy,imxy,imxz,imyz,$
         xmax,ymax,xval=xval,yval=yval,zval=zval,zof=zof,$
         zoom=zoom,lmax=lmax,xrot=xrot,zrot=zrot,dev=dev,npx=npx,npy=npy,$
         amax=amax,amin=amin,thick=thick,zpos=zpos,scale=scale,title=title,$
         length=length,xpos=xpos,ip=ip,box=box,$
         centred=centred,shell=shell,r_int=r_int,r_ext=r_ext,$
         zrr1=zrr1,zrr2=zrr2,yrr=yrr,xrr=xrr,magnify=magnify,zmagnify=zmagnify,$
         nobottom=nobottom, norm=norm, sample=sample
;
; n=15
; inxy=reform(uuu(*,*,n,0))
; imxy=reform(uuu(*,*,0,0))
; imxz=reform(uuu(*,0,*,0))
; imyz=reform(uuu(0,*,*,0))
; boxbotex_scl,inxy,imxy,imxz,imyz,1.,1.,zof=.7,zpos=.34,ip=3,amin=-.1,amax=.1,/box
;
;  the keyword box is used to overplot a bounding box
;
if ip gt 4 then print,minmax(inxy)
if ip gt 4 then print,minmax(imxy)
if ip gt 4 then print,minmax(imxz)
if ip gt 4 then print,minmax(imyz)
;
;  special version of boxbot which inputs only data on the walls
;
if n_params(0) eq 0 then begin
  print,'pro boxbot, ?????,'
  print,'  xval=xval,yval=yval,zval=zval,zof=zof,$'
  print,'  zoom=zoom,lmax=lmax,xrot=xrot,zrot=zrot,dev=dev,npx=npx,npy=npy,$'
  print,'  amax=amax,amin=amin,thick=thick,zpos=zpos,scale=scale,title=title,$'
  print,'  length=length'
  return
endif
;
; THIS VERSION SHOWS THE BOTTOM OF THE BOX ALSO. PCM 10/9/93
; NAME:         TVBOX
; PURPOSE:
;       Show a 3d array as a box of chosen size with each face showing a
;       greyscale image of the array, overlaid with the velocity fields on
;       each face.
; CATEGORY:
;       Display, graphics.
; CALLING SEQUENCE:
;       BOX, array, xvel, yvel, zvel
; INPUTS:
;       Image = the 3 dimensional array to display.
;       xvel, yvel, zvel = the 3 dimensional velocity arrays.
;       xmax, ymax = size of cube (zmax=1.0)
; OPTIONAL INPUT PARAMETERS:
;       Zoom = Multiplication factor by which to interpolate and smooth
;               out the original arrays.
;       xval, yval, zval = proportion of the array to display. Defaults
;               are zval=1, xval=yval=0. ie top, front and left side.
;       xrot = rotation angle of image about x axis ( -90 > xrot < 90 ).
;              Default 30.
;       zrot = rotation angle of image about z axis ( 0 > zrot < 90 ).
;              Default 30.
;       dev : if set the image is written to the specified device (ps or cgm).
;       lmax = scaling factor for vector field arrows. Default 0.1
;       npx, npy : no. of time box is replicated in x and y directions.
;       amax, amin : range over which colour table is scaled, for displayed
;                     array. Default max & min(array)
;       thick : thickness of arrows. Default 1 pixel.
;       xpos : position of bottom of box image in window. Default 0.0
;       zpos : position of bottom of box image in window. Default 0.0
;
; OUTPUTS:
;       No explicit outputs.
; COMMON BLOCKS:
;       None.
; SIDE EFFECTS:
;       Display is changed.
; RESTRICTIONS:
;
; PROCEDURE:
;      Use Z-buffering and the POLYFILL procedure to warp x,y and z planes
;      of a box onto the appropriate faces. Read the image from the Z-buffer
;      and copy it to the window (or Postscript device).
;      Overlay the velocity fields on each face.
; MODIFICATION HISTORY:
;      dpb, DAMTP  Apr, 1992.
;      axel, Nov 2002, adapted for plotting only scalar field; using
;        as input only the 4 bounding surfaces.
;-

on_error,2            ;Return to caller if an error occurs
s = size(imxy)         ;Get size of image
nx = s[1]
ny = s[2]
s = size(imxz)         ;Get size of image
nz = s[2]
;
if n_elements(amax) eq 0 then begin
  print,'need to specify amin and amax'
  return
endif
;
if n_elements(xrot) eq 0 then xrot=30
if n_elements(zrot) eq 0 then zrot=30
if n_elements(zoom) eq 0 then zoom=5
if n_elements(xval) eq 0 then xval=0
if n_elements(yval) eq 0 then yval=0
if n_elements(zval) eq 0 then zval=1.0
if n_elements(xpos) eq 0 then xpos=0.0
if n_elements(zpos) eq 0 then zpos=1.1
if n_elements(zof) eq 0 then zof=1.40
if n_elements(scale) eq 0 then scale=1.0
if n_elements(length) eq 0 then length=1.0
if n_elements(norm) eq 0 then norm=1.0
if n_elements(magnify) eq 0 then magnify=1.0
if n_elements(zmagnify) eq 0 then zmagnify=1.0
if n_elements(ip) eq 0 then ip=0
if keyword_set(shell) then begin            
  if n_elements(zrr1) eq 0 or n_elements(zrr2) eq 0 or $
     n_elements(xrr) eq 0 or n_elements(yrr) eq 0 then $
     print,'must pass in arrays of rr in the 4 slices, to use shell in boxbotex_scl'
endif
;
ncols=245
mincol=1
rev=1
if keyword_set(dev) then begin
  if dev eq 'ps' or dev eq 'eps' then begin
    ncols=256
    mincol=32
    rev=-1
    atop=-amin
    amin=-amax
    amax=atop
  endif else if dev eq 'cgm' then ncols=245
  dev1=dev
end else begin
  dev1=''
end
if n_elements(npx) eq 0 then npx=1
if n_elements(npy) eq 0 then npy=1
nxi=fix((1.0-xval)*(nx-1))+1
nyi=fix((1.0-yval)*(ny-1))+1
if xrot gt 0 then nzi=fix(zval*nz) else nzi=fix((1.0-zval)*(nz-1))+1
if ip gt 3 then print,'nxi,nyi,nzi=',nxi,nyi,nzi
if ip gt 3 then print,'amax,amin,xmax,ymax=',amax,amin,xmax,ymax
;
zim=bytarr(nxi*npx,nyi*npy)
zimbot=zim
if xrot gt 0 then begin
  xim=bytarr(nyi*npy,nzi)
  yim=bytarr(nxi*npx,nzi)
  ;
  ;  xz-plane
  ;
  for i = 0,npx-1 do $
    yim(i*nxi:(i+1)*nxi-1,*) = bytscl(rev*imxz(nx-nxi:nx-1,0:nzi-1)/norm $
                               ,max=amax,min=amin,top=ncols-2-mincol)+mincol
  ;
  ;  yz-plane
  ;
  for j = 0,npy-1 do $
    xim(j*nyi:(j+1)*nyi-1,*) = bytscl(rev*imyz(ny-nyi:ny-1,0:nzi-1)/norm $
                               ,max=amax,min=amin,top=ncols-2-mincol)+mincol
  ;
  ;  top and bottom xy-planes
  ;
  for i = 0,npy-1 do begin
    for j = 0,npx-1 do begin
      zim(j*nxi:(j+1)*nxi-1,i*nyi:(i+1)*nyi-1) = bytscl(rev* $
                                inxy(nx-nxi:nx-1,ny-nyi:ny-1)/norm $
                              ,max=amax,min=amin,top=ncols-2-mincol)+mincol
      zimbot(j*nxi:(j+1)*nxi-1,i*nyi:(i+1)*nyi-1) = bytscl(rev* $
                                imxy(nx-nxi:nx-1,ny-nyi:ny-1)/norm $
                              ,max=amax,min=amin,top=ncols-2-mincol)+mincol
    endfor
  endfor
  zimg = rebinbox(reform(zim,nxi*npx,nyi*npy), zoom, sample=sample)
  zimgbot = rebinbox(reform(zimbot,nxi*npx,nyi*npy), zoom, sample=sample)
  yimg = rebinbox(reform(yim,nxi*npx,nzi), zoom, /zdir, sample=sample)
  ximg = rebinbox(reform(xim,nyi*npy,nzi),zoom, /zdir, sample=sample)
  ;
  ; set up masking for transparency, if using shell
  ;
  if keyword_set(shell) then begin
    zrrg = rebinbox(reform(zrr1,nxi*npx,nyi*npy), zoom, sample=sample)
    zrr2g = rebinbox(reform(zrr2,nxi*npx,nyi*npy), zoom, sample=sample)
    yrrg = rebinbox(reform(yrr,nxi*npx,nzi), zoom, /zdir, sample=sample)
    xrrg = rebinbox(reform(xrr,nyi*npy,nzi), zoom, /zdir, sample=sample)
    indz=where(zrrg lt r_int or zrrg gt r_ext,nindz)
    indz2=where(zrr2g lt r_int or zrr2g gt r_ext,nindz2)
    indx=where(xrrg lt r_int or xrrg gt r_ext,nindx)
    indy=where(yrrg lt r_int or yrrg gt r_ext,nindy)
    if nindz ne 0 then zimg(indz)=mincol-1
    if nindz2 ne 0 then zimgbot(indz2)=mincol-1
    if nindx ne 0 then ximg(indx)=mincol-1
    if nindy ne 0 then yimg(indy)=mincol-1
  endif
  ;
  ;  horizontal velocity in the two xy-planes
  ;
  ;uz    = rebinbox(xwxy,zoom)
  ;vz    = rebinbox(ywxy,zoom)
  ;uzbot = rebinbox(xvxy,zoom)
  ;vzbot = rebinbox(yvxy,zoom)
endif else begin
  print,"xrot < 0: don't have data for this"
  return
endelse
;
;vy = rebinbox(xvxz,zoom,/zdir)
;wy = rebinbox(zvxz,zoom,/zdir)
;ux = rebinbox(yvyz,zoom,/zdir)
;wx = rebinbox(zvyz,zoom,/zdir)
;
;ux1=fltarr(ny*npy*zoom,(nz-1)*zoom+1)
;wx1=ux1
;vy1=fltarr(nx*npx*zoom,(nz-1)*zoom+1)
;wy1=vy1
;uz1=fltarr(nx*npx*zoom,ny*npy*zoom)
;vz1=uz1
;uz1bot=uz1
;vz1bot=uz1
;for i = 0,npy-1 do begin
;  ux1(i*ny*zoom:(i+1)*ny*zoom-1,*)=ux(*,*)
;  wx1(i*ny*zoom:(i+1)*ny*zoom-1,*)=wx(*,*)
;endfor
;for j = 0,npx-1 do begin
;  vy1(j*nx*zoom:(j+1)*nx*zoom-1,*)=vy(*,*)
;  wy1(j*nx*zoom:(j+1)*nx*zoom-1,*)=wy(*,*)
;endfor
;for i = 0,npx-1 do begin
;  for j = 0,npy-1 do begin
;    uz1(i*nx*zoom:(i+1)*nx*zoom-1,j*ny*zoom:(j+1)*ny*zoom-1)=uz(*,*)
;    vz1(i*nx*zoom:(i+1)*nx*zoom-1,j*ny*zoom:(j+1)*ny*zoom-1)=vz(*,*)
;    uz1bot(i*nx*zoom:(i+1)*nx*zoom-1,j*ny*zoom:(j+1)*ny*zoom-1)=uzbot(*,*)
;    vz1bot(i*nx*zoom:(i+1)*nx*zoom-1,j*ny*zoom:(j+1)*ny*zoom-1)=vzbot(*,*)
;  endfor
;endfor

;help,zimg,ximg,yimg,ux,vy,uz
if ip gt 4 then print,max(ximg),max(yimg),max(zimg),min(zimg)
zs=size(zimg)
xs=size(ximg)
ys=size(yimg)

;if keyword_set(dev) then begin
;  x_size=640 & y_size=480
;  ;  x_size=200 & y_size=200
;endif else begin
  x_size=!d.x_size & y_size=!d.y_size
;endelse
;
if dev1 ne 'z' then begin
  set_plot,'z'
  device,set_res=[x_size,y_size]
end
erase
;
;Set up scaling
sf=1.25/scale
if (xmax*npx ge ymax*npy) then maxscale=xmax*npx*sf else maxscale=ymax*npy*sf
;AB: maxscale adjusted via inverse magnify parameter
maxscale=maxscale/magnify
range=[0,maxscale]
scale3,xrange=range,yrange=range,zrange=range/zmagnify,$
        ax=(360+xrot) mod 360,az=zrot
z0=zpos+(sf-1)/(2*sf)*maxscale & z1=z0+zval
x0=xpos+(sf-1)/(2*sf)*maxscale & x1=x0+xmax*npx+xval
;x0=xval+(sf-1)/(2*sf)*maxscale & x1=x0+xmax*npx
y0=yval+(sf-1)/(2*sf)*maxscale & y1=y0+ymax*npy
;
; set up verts for planes on edges of box, or through centre
;
if not keyword_set(centred) then begin
  ; traditional edge-planes of box
  verts=[[x0,y0,z0],[x1,y0,z0],[x1,y1,z0],[x0,y1,z0],$
         [x0,y0,z1],[x1,y0,z1],[x1,y1,z1],[x0,y1,z1],$
         [x0,y0,z0-zof],[x1,y0,z0-zof],[x1,y1,z0-zof],[x0,y1,z0-zof]]

  polyfill,verts(*,[3,0,4,7]),/t3d,pattern=ximg, $
        image_coord=[[xs(1)-1,0],[0,0],[0,xs(2)-1],[xs(1)-1,xs(2)-1]]
  polyfill,verts(*,[0,1,5,4]),/t3d,pattern=yimg, $
        image_coord=[[0,0],[ys(1)-1,0],[ys(1)-1,ys(2)-1],[0,ys(2)-1]]
  if xrot gt 0 then begin
    polyfill,verts(*,[4,5,6,7]),/t3d,pattern=zimg, $
        image_coord=[[0,0],[zs(1)-1,0],[zs(1)-1,zs(2)-1],[0,zs(2)-1]] 
    if not keyword_set(nobottom) then begin
      polyfill,verts(*,[8,9,10,11]),/t3d,pattern=zimgbot, $
          image_coord=[[0,0],[zs(1)-1,0],[zs(1)-1,zs(2)-1],[0,zs(2)-1]] 
    endif
  endif else begin
    polyfill,verts(*,[0,1,2,3]),/t3d,pattern=zimg, $
        image_coord=[[0,0],[zs(1)-1,0],[zs(1)-1,zs(2)-1],[0,zs(2)-1]]
  endelse
endif else begin
  ; centre-slices through box; with transparency if shell options set
  xm=x0+0.5*xmax*npx+xval & ym=y0+0.5*ymax*npy & zm=z0+0.5*zval
  verts=[[xm,y0,z0],[x1,ym,z0],[xm,y1,z0],[x0,ym,z0],$
         [x0,y0,zm],[x1,y0,zm],[x1,y1,zm],[x0,y1,zm],$
         [xm,y0,z1],[x1,ym,z1],[xm,y1,z1],[x0,ym,z1],$
         [x0,y0,z0-zof],[x1,y0,z0-zof],[x1,y1,z0-zof],[x0,y1,z0-zof]]
  polyfill,verts(*,[2,0,8,10]),/t3d,pattern=ximg, transparent=mincol, $
        image_coord=[[xs(1)-1,0],[0,0],[0,xs(2)-1],[xs(1)-1,xs(2)-1]]
  polyfill,verts(*,[3,1,9,11]),/t3d,pattern=yimg, transparent=mincol, $
        image_coord=[[0,0],[ys(1)-1,0],[ys(1)-1,ys(2)-1],[0,ys(2)-1]]
  if xrot gt 0 then begin
    polyfill,verts(*,[4,5,6,7]),/t3d,pattern=zimg, transparent=mincol, $
        image_coord=[[0,0],[zs(1)-1,0],[zs(1)-1,zs(2)-1],[0,zs(2)-1]] 
    if not keyword_set(nobottom) then begin
      polyfill,verts(*,[12,13,14,15]),/t3d,pattern=zimgbot, transparent=mincol, $
          image_coord=[[0,0],[zs(1)-1,0],[zs(1)-1,zs(2)-1],[0,zs(2)-1]] 
    endif
  endif else begin
    polyfill,verts(*,[0,1,2,3]),/t3d,pattern=zimgbot, $
        image_coord=[[0,0],[zs(1)-1,0],[zs(1)-1,zs(2)-1],[0,zs(2)-1]]
  endelse

endelse

a=tvrd()
;
;  make background white
;
bad=where(a eq 0) & a(bad)=255
;
;  plot
;
;if keyword_set(dev) then begin
; if dev eq 'ps' then begin
;   set_plot,dev
;   device,xsize=18,ysize=12,xoff=1.5,yoff=10,bits=8,encap=0
;   a(where(a eq 0B))=255B
;   tv,a,xsize=1.0,ysize=1.0,/norm
;   if keyword_set(title) then begin 
;      ;device,/TIMES
;      ;xyouts,0.18,0.82,title,font=0,size=2.5,/norm
;      xyouts,0.97,1.34,title,siz=2.
;   endif
; endif else begin
;   if dev eq 'eps' then begin
;     set_plot,'ps'
;     device,file='idl.eps',/encap
;     device,xsize=18,ysize=18,xoff=1.5,yoff=10,bits=8
;     a(where(a eq 0B))=255B
;     tv,a,xsize=1.0,ysize=1.0,/norm
;     if keyword_set(title) then begin 
;        device,/TIMES
;        xyouts,0.18,0.82,title,font=0,size=2.5,/norm
;     endif
;   endif else begin
;     if dev eq 'cgm' then begin
;       set_plot,dev
;       device,colors=ncols
;;       restore,'~dpb/p3d/col.tab3'
;;       tvlct,r,g,b
;;       loadct,file='/local/d/pcm/idl/colors1.tbl'
;;       mloadct
;;       help,/dev
;       tv,a
;     endif
;   endelse
; endelse
;endif else begin
  ;
  ;  default
  ;
  if dev1 ne 'z' then begin
    set_plot,'x'
    tv,a
  end
   if keyword_set(title) then begin 
      ;xyouts,0.18,0.82,title,size=2.5,color=0,/norm
      xyouts,1.02,1.34,title,siz=2.,col=122
   endif
;endelse
;
;  overplot streak lines of the vector field
;
if n_elements(lmax) eq 0 then lmax=0.1
if xrot gt 0 then begin
   ;
   ;  top surface
   ;
   if keyword_set(box) then velfld_box,verts,face=0,thick=thick
   ;
   ;  plot vectors at the bottom
   ;
  if keyword_set(box) then velfld_box,verts,face=4,thick=thick,zof=zof
endif else begin
  if keyword_set(box) then velfld_box,verts,face=4,thick=thick
endelse
;
;  plot vectors on the side walls
;  right hand side: xz-plane
;
if keyword_set(box) then velfld_box,verts,face=1,thick=thick
;
;  left hand side: yz-plane
;
if keyword_set(box) then velfld_box,verts,face=2,thick=thick
;
;  reset device
;
if keyword_set(dev) then begin
  if dev1 ne 'z' and dev1 ne 'x' then begin
    device,/close
    set_plot,'x'
    print,'image in idl.',dev,', device X reset'
  end
endif

end
