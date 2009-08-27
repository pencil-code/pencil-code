pro boxbotex_scl_sph,imrp,imrt1,imrt2,imrt3,imrt4,$
;         xmax,ymax,rad=rad,tht=tht,phi=phi,xc=xc,yc=yc,ip=ip,$
;         zof=zof,zpos=zpos,$
;         amin=amin,amax=amax,dev=dev,$
;         shell=shell,centred=centred,scale=scale,$
;         r_int=r_int,r_ext=r_ext,trr=trr,prr=prr,$
;         nobottom=nobottom,norm=norm,xrot=xrot,zrot=zrot,zoomz=zoomz
         xmax,ymax,xval=xval,yval=yval,zval=zval,zof=zof,$
         zoom=zoom,lmax=lmax,xrot=xrot,zrot=zrot,dev=dev,npx=npx,npy=npy,$
         amax=amax,amin=amin,thick=thick,zpos=zpos,scale=scale,title=title,$
         length=length,xpos=xpos,ip=ip,box=box,$
         r_int=r_int,r_ext=r_ext,orig_aspect=orig_aspect,$
         trr=trr,prr=prr,ptt=ptt,magnify=magnify,$
         nobottom=nobottom,norm=norm,rad=rad,tht=tht,phi=phi,$
         xc=xc,yc=yc,xm=xm,zm=zm,zoomz=zoomz,nointerpz=nointerpz

;
nr=n_elements(xm)
nt=n_elements(tht)
r0=xm[0]  & rn=xm[nr-1] & Lr=rn-r0
tht0=tht[0] & thtn=tht[nt-1]

;swap y and z
nx=n_elements(xc) 
ny=n_elements(yc)
nz=n_elements(zm)

Lx=xc[nx-1]-xc[0]
Ly=yc[ny-1]-yc[0]

;
;  the keyword box is used to overplot a bounding box
;
if ip gt 4 then print,minmax(imrp)
if ip gt 4 then print,minmax(imrt1)
if ip gt 4 then print,minmax(imrt2)
if ip gt 4 then print,minmax(imrt3)
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
;
;
;swap y-z here, and let's thing cartesian - phi=y, tht=z
;
on_error,2            ;Return to caller if an error occurs
s = size(imrp)        ;Get size of image
nx = s[1]
ny = s[2]
s = size(imrt1)       ;Get size of image
nr = s[1]
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
if n_elements(ip) eq 0 then ip=0
if ((n_elements(trr) eq 0) or (n_elements(prr) eq 0)) then $
  print,'must pass in arrays of rr in the 5 slices, to use shell in boxbotex_scl'
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
npr=npx
nxi=fix((1.0-xval)*(nx-1))+1
nyi=fix((1.0-yval)*(ny-1))+1
nri=fix((1.0-xval)*(nr-1))+1
if xrot gt 0 then nzi=fix(zval*nz) else nzi=fix((1.0-zval)*(nz-1))+1
if ip gt 3 then print,'nxi,nyi,nri,nzi=',nxi,nyi,nri,nzi
if ip gt 3 then print,'amax,amin,xmax,ymax=',amax,amin,xmax,ymax
;
thtim=bytarr(nxi*npx,nyi*npy)
if xrot gt 0 then begin
  phi1im=bytarr(nri*npr,nzi)
  phi2im=bytarr(nri*npr,nzi)
  phi3im=bytarr(nri*npr,nzi)
  phi4im=bytarr(nri*npr,nzi)
;
;  rtheta wedges
;
  for i=0,npx-1 do begin
   phi1im(i*nri:(i+1)*nri-1,*) = bytscl(rev*imrt1(nr-nri:nr-1,0:nzi-1)/norm $
                                 ,max=amax,min=amin,top=ncols-2-mincol)+mincol
;
   phi2im(i*nri:(i+1)*nri-1,*) = bytscl(rev*imrt2(nr-nri:nr-1,0:nzi-1)/norm $
                                 ,max=amax,min=amin,top=ncols-2-mincol)+mincol
;
   phi3im(i*nri:(i+1)*nri-1,*) = bytscl(rev*imrt3(nr-nri:nr-1,0:nzi-1)/norm $
                                 ,max=amax,min=amin,top=ncols-2-mincol)+mincol
;
   phi4im(i*nri:(i+1)*nri-1,*) = bytscl(rev*imrt4(nr-nri:nr-1,0:nzi-1)/norm $
                                 ,max=amax,min=amin,top=ncols-2-mincol)+mincol
  endfor
;
;  top and bottom xy-planes
;
  for i = 0,npy-1 do begin
    for j = 0,npx-1 do begin
      thtim(j*nxi:(j+1)*nxi-1,i*nyi:(i+1)*nyi-1) = bytscl(rev* $
                                imrp(nx-nxi:nx-1,ny-nyi:ny-1)/norm $
                              ,max=amax,min=amin,top=ncols-2-mincol)+mincol
    endfor
  endfor
;
; Use simple rebin instead of the cryptic rebinbox
;
  thtimg =rebin(reform(thtim ),nxi*npx*zoom,nyi*npy*zoom)
  phi1img=rebin(reform(phi1im),nri*npr*zoom,nzi*zoom)
  phi2img=rebin(reform(phi2im),nri*npr*zoom,nzi*zoom)
  phi3img=rebin(reform(phi3im),nri*npr*zoom,nzi*zoom)
  phi4img=rebin(reform(phi4im),nri*npr*zoom,nzi*zoom)
;
; set up masking for transparency, if using shell
;
  trrg = rebin(reform(trr),nxi*npx*zoom,nyi*npy*zoom)
  prrg = rebin(reform(prr),nri*npr*zoom,nzi*zoom)
  indt=where(trrg lt r_int or trrg gt r_ext,nindt)
  indp=where(prrg lt r_int or prrg gt r_ext,nindp)
  if nindt ne 0 then thtimg(indt)=mincol-1
  if nindp ne 0 then begin
    phi1img(indp)=mincol-1
    phi2img(indp)=mincol-1
    phi3img(indp)=mincol-1
    phi4img(indp)=mincol-1
  endif
;
; set up theta masking
; 
  pttg=rebin(reform(ptt),nri*npr*zoom,nzi*zoom)
  indpt=where(pttg lt tht0 or pttg gt thtn,nindpt)
  if nindpt ne 0 then begin
    phi1img(indpt)=mincol-1
    phi2img(indpt)=mincol-1
    phi3img(indpt)=mincol-1
    phi4img(indpt)=mincol-1
  endif
  ;
endif else begin
  print,"xrot < 0: don't have data for this"
  return
endelse
;
if ip gt 4 then print,max(thtimg),max(phi1img),max(phi2img),min(phi3img),min(phi4img)
zs=size(thtimg)
xs=size(phi1img)
ys=xs
;
x_size=!d.x_size & y_size=!d.y_size
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
scale3,xrange=range,yrange=range,zrange=range,$
        ax=(360+xrot) mod 360,az=zrot
z0=zpos+(sf-1)/(2*sf)*maxscale & z1=z0+zval
x0=xpos+(sf-1)/(2*sf)*maxscale & x1=x0+xmax*npx+xval
;x0=xval+(sf-1)/(2*sf)*maxscale & x1=x0+xmax*npx
y0=yval+(sf-1)/(2*sf)*maxscale & y1=y0+ymax*npy
;
; vertices of the center of the box
;
xm=x0+0.5*xmax*npx+xval & ym=y0+0.5*ymax*npy & zm=z0+0.5*zval
;
; z of the spherical shell
;
  Lxbox=x1-x0 & xratio=Lxbox/Lx
  if (keyword_set(nointerpz)) then begin
    spc=r_int*xratio
   endif else begin
    spc=r_int*sin(tht0)*xratio
  endelse
  spc2=(xc[nx-1]-r_ext)*xratio
  yf0=ym+spc & yf1=y1-spc2 & yfm=yf0+(yf1-yf0)/2.

  Lrbox=yf1-yf0 & rratio=Lrbox/Lr

  if (keyword_set(nointerpz)) then begin
    ; define four different vertices, obeying the
    ; opening angle of the spherical wedges, and the
    ; aspect ratio of the axes. Polyfill will attach 
    ; the image inside the polygon defined by them
    zu1=zm+r0*cos(tht[0])*rratio*zoomz 
    zd1=zm+r0*cos(tht[nt-1])*rratio*zoomz
    zu2=zm+rn*cos(tht[0])*rratio*zoomz 
    zd2=zm+rn*cos(tht[nt-1])*rratio*zoomz
  endif else begin
    if (not keyword_set(orig_aspect)) then begin
      ; define only up and down, obeying the aspect ratio
      ; of the axes. The rectangle will be masked later
      zu1=zm+rn*cos(tht[0])*rratio*zoomz    & zu2=zu1
      zd1=zm+rn*cos(tht[nt-1])*rratio*zoomz & zd2=zd1
    endif else begin
      ; don't care about the aspect ratio, use z0 and z1
      ; this is more like what rvid_box does
      zu1=z1 & zu2=zu1
      zd1=z0 & zd2=zd1
    endelse
  endelse
;
; set up verts for planes on edges of box, or through centre
;
;x=0 plot - phi=pi/2
  yf0=ym+spc & yf1=y1-spc2 & yfm=yf0+(yf1-yf0)/2.
  vertxs=[[xm,yf1,zd2],[xm,yf0,zd1],$
          [xm,yf0,zu1],[xm,yf1,zu2]]
  polyfill,vertxs,/t3d,pattern=phi4img, transparent=mincol, $
    image_coord=[[xs(1)-1,0],[0,0],[0,xs(2)-1],[xs(1)-1,xs(2)-1]]

;second x plot - phi=3pi/2
  yf0=ym-spc & yf1=y0+spc2 & yfm=yf0+(yf1-yf0)/2.
  vertxs=[[xm,yf1,zd2],[xm,yf0,zd1],$
          [xm,yf0,zu1],[xm,yf1,zu2]]
  polyfill,vertxs,/t3d,pattern=phi2img, transparent=mincol, $
    image_coord=[[xs(1)-1,0],[0,0],[0,xs(2)-1],[xs(1)-1,xs(2)-1]]

;y=0 plot - phi=0
  xf0=xm+spc & xf1=x1-spc2 & xfm=xf0+(xf1-xf0)/2.
  vertys=[[xf0,ym,zd1],[xf1,ym,zd2],$
          [xf1,ym,zu2],[xf0,ym,zu1]]
  polyfill,vertys,/t3d,pattern=phi3img, transparent=mincol, $
    image_coord=[[0,0],[ys(1)-1,0],[ys(1)-1,ys(2)-1],[0,ys(2)-1]]

;second y plot - phi=-pi
  xf0=xm-spc & xf1=x0+spc2 & xfm=xf0+(xf1-xf0)/2.
  vertys=[[xf0,ym,zd1],[xf1,ym,zd2],$
          [xf1,ym,zu2],[xf0,ym,zu1]]
  polyfill,vertys,/t3d,pattern=phi1img, transparent=mincol, $
    image_coord=[[0,0],[ys(1)-1,0],[ys(1)-1,ys(2)-1],[0,ys(2)-1]]

;midplane plots
verts=[[xm,y0,z0],[x1,ym,z0],[xm,y1,z0],[x0,ym,z0],$
         [x0,y0,zm],[x1,y0,zm],[x1,y1,zm],[x0,y1,zm],$
         [xm,y0,z1],[x1,ym,z1],[xm,y1,z1],[x0,ym,z1],$
         [x0,y0,z0-zof],[x1,y0,z0-zof],[x1,y1,z0-zof],[x0,y1,z0-zof]]
if xrot gt 0 then begin
    ;center plot
    polyfill,verts(*,[4,5,6,7]),/t3d,pattern=thtimg, transparent=mincol, $
      image_coord=[[0,0],[zs(1)-1,0],[zs(1)-1,zs(2)-1],[0,zs(2)-1]] 
    if not keyword_set(nobottom) then begin
        ;bottom plot
        polyfill,verts(*,[12,13,14,15]),/t3d,pattern=thtimg, transparent=mincol, $
          image_coord=[[0,0],[zs(1)-1,0],[zs(1)-1,zs(2)-1],[0,zs(2)-1]] 
    endif
endif else begin
    polyfill,verts(*,[0,1,2,3]),/t3d,pattern=thtimg, $
      image_coord=[[0,0],[zs(1)-1,0],[zs(1)-1,zs(2)-1],[0,zs(2)-1]]
endelse

a=tvrd()
;
;  make background white
;
bad=where(a eq 0) & a(bad)=255
;
;  plot
;
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
