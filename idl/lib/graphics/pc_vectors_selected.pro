;  $Id$
;
PRO pc_vectors_selected,l,m,n,U,V,W,X,Y,Z,length=length,ax=ax,az=az,nbox=nbox,$
	color=color,symsize=symsize,scale=scale,black=black,field=field,$
	per=per,center=center,nobox=nobox,noarrow=noarrow,back=back
;
if n_params(0) eq 0 then begin
  print,'PRO VECGD,UVW,X,Y,Z,length=length,ax=ax,az=az,nbox=nbox,color=color,symsize=symsize,scale=scale'
  return
end

;
;+
; NAME:
;     VECGD
; PURPOSE:
;     Produce a three dimensional velocity field plot.
;           Draw a directed arrow at each point showing
;           the direction and magnitude of the field.
;           
; CATEGORY:
;     Plotting
; CALLING SEQUENCE:
;     VECGD,U,V,W  ...  or:
;     VECGD,U,V,W,X,Y,Z
; INPUTS:
;     U = The X component of the 3 dimensional field.  Both
;           U, V and W must be 3 dimensional and have the same
;           dimensions.
;     V = The Y component of the 3 dimensional field.
;     W = The Z component of the 3 dimensional field.
;           The vector at point (i,j) has the magnitude of
;           (U(i,j)^2 + V(i,j)^2) + W(i,j)^2)^0.5
; OPTIONAL INPUT PARAMETERS:
;     X = Optional abcissae values, Size of X must = 1st dimension
;           of U, V and W.  Must be a vector. [0,nx]
;     Y = Optional ordinate values.  Same rules as X.
;     Z = Optional ordinate values.  Same rules as X.
;     length    = length (in grid units) of longest vector [1.0]
;     ax, az = view from angle ax, az rotation about x, z axis [30.,30.]
;     nbox   = 0, erase, new box and new perspective
;              1, noerase, same perspective
;              2, erase, new perspective + ?
;              other, noerase, new perspective, no box
;     color  = color of vectors [1]
; OUTPUTS:
;     None.
; COMMON BLOCKS:
;     None.
; SIDE EFFECTS:
;     Plotting on the selected device is performed.  System
;     variables concerning Plotting are changed.
; RESTRICTIONS:
;     None.
; PROCEDURE:
;     Straightforward.  The system variables !Xtitle, !Ytitle and
;     !Mtitle may be set to title the axes.
; MODIFICATION HISTORY:
;     Dms, Rsi, Oct, 1983.
;     For Sun, DMS, RSI, April, 1989.
;     RFS, Oct 1991, converted to 3D, eliminated arrowheads
;     AB, Oct 1992, symsize as free parameter
;     AB, May 1993, scale to focus on the center region, different fields
;-
;  colors
;
      if n_elements(color) eq 0 then begin
        autocol=1
        color=!d.n_colors
        color=255
      endif else begin
        autocol=0
      endelse
      if n_elements(length) eq 0 then length=3.0
      if n_elements(thresh) eq 0 then thresh=0.0
      if n_elements(ax)     eq 0 then ax=30
      if n_elements(az)     eq 0 then az=30
      if n_elements(nbox)   eq 0 then nbox=0
      if n_elements(symsize) eq 0 then symsize=1.0
      if n_elements(scale)  eq 0 then scale=1.
      if n_elements(field)  eq 0 then field=0
      blck=0 & if n_elements(black) eq 1 then blck=1
      arrw=1 & if n_elements(noarrow) eq 1 then arrw=0
;
;  set color
;
;if (field eq 0) then icol0=[0,63]
;if (field eq 1) then icol0=[1,30]
;if (field eq 2) then icol0=[32,30]
;
;if (field eq -1) then icol0=[0,256]
;if (field eq 0) then icol0=[0,!d.n_colors]
;if (field eq 1) then icol0=[1,!d.n_colors/2-1]
;if (field eq 2) then icol0=[!d.n_colors/2+1,!d.n_colors/2-1]
;
if (field eq -1) then icol0=[0,256]
if (field eq 0) then icol0=[0,255]
if (field eq 1) then icol0=[1,127]
if (field eq 2) then icol0=[128,254]
;
      s = size(u)
      t = size(v)
      q = size(w)
      ;
      mag = sqrt(u^2+v^2+w^2)       ;magnitude.
      ;
      ;print,minmax(mag)
      ugood = U
      vgood = V
      wgood = W
      x0 = min(x)                   ;get scaling
      x1 = max(x)
      y0 = min(y)
      y1 = max(y)
      z0 = min(z)
      z1 = max(z)
      rmin=[x0,y0,z0]
      rmax=[x1,y1,z1]
      ;
      maxmag0 = max(mag)
      minmag0 = min(mag)
;     if (length le 0) then begin
;       maxmag=(mag/abs(length) > 1e-22)
;     endif else begin
;       maxmag = max(mag)/length
;     endelse
      ;dx = (x1-x0)/maxmag*ugood      ;components.
      ;dy = (y1-y0)/maxmag*vgood
      ;dz = (z1-z0)/maxmag*wgood

maxmag=length/maxmag0

      dx = maxmag*ugood      ;components.
      dy = maxmag*vgood
      dz = maxmag*wgood

      xstylesav=!x.style
      ystylesav=!y.style
      zstylesav=!z.style
      !x.style=1                    ;exact box dimensions
      !y.style=1
      !z.style=1
;
; make box, set perspective using t3dset
; specht
;     box,rmin,rmax,ax,az,nosurf=nbox,scale=scale,per=per
if nbox eq 0 then begin
     ;       erase
     ;       scal3d,max(x),max(y),max(z),ax=ax,az=az,per=per,scl=scale,$
     ;         center=center,nobox=nobox
  xyzr=[min([x,y,z]),max([x,y,z])]
     ;surface,fltarr(2,2),/nodata,/save,xr=minmax(x),yr=minmax(y),zr=minmax(z),ax=ax,az=az
  if keyword_set(back) then begin
    col0=1
    col0=255-back
    surface,fltarr(2,2),/nodata,/save,xr=xyzr,yr=xyzr,zr=xyzr,ax=ax,az=az,back=back,col=col0,xtit='x',ytit='y',ztit='z'
  end else begin
    col0=255
    surface,fltarr(2,2),/nodata,/save,xr=xyzr,yr=xyzr,zr=xyzr,ax=ax,az=az,xtit='x',ytit='y',ztit='z'
  end
x0=min(x) & x1=max(x)
y0=min(y) & y1=max(y)
z0=min(z) & z1=max(z)
        plots,[x0,x1,x1],[y1,y1,y0],[z0,z0,z0],/t3d,/data,col=col0
        plots,[x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],[z1,z1,z1,z1,z1],/t3d,/data,col=col0
        plots,[x0,x0],[y0,y0],[z0,z1],/t3d,/data,col=col0
        plots,[x1,x1],[y0,y0],[z0,z1],/t3d,/data,col=col0
        plots,[x1,x1],[y1,y1],[z0,z1],/t3d,/data,col=col0
;
      endif
;
; draw vectors
; without arrowheads (for now)
;
      ;print,'SPECHT1, arrw=',arrw
      for i=0L,n_elements(l)-1L do begin ;Each point
;print,l(i),m(i),n(i),dx(i),dy(i),dz(i)
        if arrw eq 0 then begin
          x0 = l(i)
          x1 = x0 + dx(i)
          y0 = m(i)
          y1 = y0 + dy(i)
          z0 = n(i)
          z1 = z0 + dz(i)
        endif else begin
          ct=.10*symsize      ;(opening angle of arrow)
          st=1.-.15*symsize   ;(length of arrow)
          x0 = l(i)
          x1 = x0 + dx(i)
          x2 = x0 + st*dx(i) + ct*dz(i)
          x3 = x0 + st*dx(i) + ct*dy(i)
          y0 = m(i)
          y1 = y0 + dy(i)
          y2 = y0 + st*dy(i)
          y3 = y0 + st*dy(i) - ct*dx(i)
          z0 = n(i)
          z1 = z0 + dz(i)
          z2 = z0 + st*dz(i) - ct*dx(i)
          z3 = z0 + st*dz(i)
        endelse
	  ;
	  ;  set color and vector length
	  ;
          if autocol eq 0 then col=color else $
                               col=fix(icol0(0)+(icol0(1)-icol0(0))*(mag(i)-minmag0)/(maxmag0-minmag0))
          ;
          ;sym=symsize*mag(i)/maxmag0
;
;base point (instead of arrows, at least for now)
;          plots,[x0,x0],[y0,y0],[z0,z0],/noclip,/t3d,/data,color=col,$
;            psym=5,symsize=sym
;
;  symsize should be less than 1.5
;
                             ; plot vectors, base to head
;         plots,[x0,x1],[y0,y1],[z0,z1],/noclip,/t3d,/data,color=col
;                            ; 2D arrows
;         plots,[x0,x1,x1-(ct*dx-st*dy),x1,x1-(ct*dx+st*dy)], $
;               [y0,y1,y1-(ct*dy+st*dx),y1,y1-(ct*dy-st*dx)], $
;               [z0,z1,z1-(ct*dz+st*dx),z1,z1-(ct*dz-st*dx)], $
;               /noclip,/t3d,/data,color=col
;
;         plots,[x0,x1,x2],[y0,y1,y2],[z0,z1,z2],/noclip,/t3d,/data,color=col
;                            ; 2D arrows
;
        if arrw eq 0 then begin
if blck eq 0 then plots,[x0,x1],[y0,y1],[z0,z1],/noclip,/t3d,/data,color=col
if blck eq 0 then plots,[x1,x1],[y1,y1],[z1,z1],/noclip,/t3d,/data,color=col,ps=4
   ;if blck eq 0 then plots,[x0,x0],[y0,y0],[z0,z0],/noclip,/t3d,/data,color=col,ps=4
if blck eq 1 then plots,[x0,x1],[y0,y1],[z0,z1],/noclip,/t3d,/data
        endif else begin
if blck eq 0 then plots,[x0,x1,x2,x1,x3],[y0,y1,y2,y1,y3],[z0,z1,z2,z1,z3],/noclip,/t3d,/data,color=col
if blck eq 1 then plots,[x0,x1,x2,x1,x3],[y0,y1,y2,y1,y3],[z0,z1,z2,z1,z3],/noclip,/t3d,/data
        endelse
                             ; 2D arrows
      endfor

      !x.style=xstylesav
      !y.style=ystylesav
      !z.style=zstylesav
      return
end
