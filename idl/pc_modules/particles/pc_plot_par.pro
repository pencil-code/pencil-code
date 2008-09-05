;
;  $Id$
;
pro pc_plot_par, xx, $
    x=x, y=y, z=z, com=com, shiftx=shiftx, shifty=shifty, shiftz=shiftz, $
    position=position, ps=ps, color=color, pcolor=pcolor, fcolor=fcolor, $
    noerase=noerase, drawgrid=drawgrid, thick=thick, sym=sym, $
    ax=ax, az=az, $
    filename=filename, imgdir=imgdir, datadir=datadir, $
    lxy=lxy, lxz=lxz, lyz=lyz, quiet=quiet

default, x, 0 & default, y, 0 & default, z, 0
default, shiftx, 0 & default, shifty, 0 & default, shiftz, 0
default, com, 0
default, position, [0.1,0.1,0.9,0.9]
default, ps, 0
default, color, 0
default, noerase, 0
default, drawgrid, 0
default, thick, 3.0
default, sym, 3
default, filename, 'particles.eps'
default, imgdir, '.'
default, datadir, './data/'
default, lxy, 0
default, lxz, 0
default, lyz, 0
default, ax, 45
default, az, 45
default, quiet, 0

if (n_elements(x) ne 1 or n_elements(y) ne 1 or n_elements(z) ne 1) then begin
  nx=n_elements(x) & ny=n_elements(y) & nz=n_elements(z)
  x0=x[0] & x1=x[nx-1]
  y0=y[0] & y1=y[ny-1]
  z0=z[0] & z1=z[nz-1]
endif else begin
  pc_read_dim, obj=dim, datadir=datadir, /quiet
  pc_read_param, obj=par, datadir=datadir, /quiet
  nx=dim.nxgrid & ny=dim.nygrid & nz=dim.nzgrid
  x0=par.xyz0[0] & x1=x0+par.Lxyz[0]
  y0=par.xyz0[1] & y1=y0+par.Lxyz[1]
  z0=par.xyz0[2] & z1=z0+par.Lxyz[2]
endelse
;;
;;  Force 2-D plane.
;;
if (lxy) then nz=1
if (lxz) then ny=1
if (lyz) then nx=1

if (shifty ne 0.0) then begin
  for k=0L,n_elements(xx[*,1])-1 do begin
    xx[k,1]=xx[k,1]+shifty
    if (xx[k,1] gt y1) then xx[k,1]=xx[k,1]-(y1-y0)
    if (xx[k,1] lt y0) then xx[k,1]=xx[k,1]+(y1-y0)
  endfor
endif

if (com eq 1) then begin
  x0=x0+mean(xx[*,0]) & x1=x1+mean(xx[*,0])
  y0=y0+mean(xx[*,1]) & y1=y1+mean(xx[*,1])
  z0=z0+mean(xx[*,2]) & z1=z1+mean(xx[*,2])
endif

if ( (nx ne 1) and (ny ne 1) and (nz ne 1) ) then begin
  xsize=14.0
  ysize=11.0
endif else begin
  xsize=12.0
  ysize=12.0
endelse

if (ps) then begin
  set_plot, 'ps'
  device, /encapsulated, color=color, xsize=xsize, ysize=ysize, $
      font_size=11, filename=imgdir+'/'+filename
endif

pcharthick=!p.charthick & pthick=!p.thick
xthick=!x.thick & ythick=!y.thick & zthick=!z.thick

!p.charthick=thick & !p.thick=thick
!x.thick=thick & !y.thick=thick & !z.thick=thick

if (color) then loadct, 12
default, fcolor, 100
default, pcolor, 200

if ( (nx ne 1) and (ny ne 1) and (nz ne 1) ) then begin

  surface,  [[0.0,0.0,0.0],[0.0,0.0,0.0]], col=fcolor, $
      xrange=[x0,x1], yrange=[y0,y1], zrange=[z0,z1], $
      xstyle=1, ystyle=1, zstyle=1, /save, /nodata, position=position, $
      xtickformat='(A1)', ytickformat='(A1)', ztickformat='(A1)', $
      ax=ax, az=az, noerase=noerase
  axis, xaxis=1, x0, y1, z1, /t3d, xtickformat='(A1)',col=fcolor
  axis, xaxis=1, x0, y1, z0, /t3d, xtickformat='(A1)',col=fcolor
  axis, yaxis=1, x1, y0, z0, /t3d, ytickformat='(A1)',col=fcolor
  axis, yaxis=1, x1, y0, z1, /t3d, ytickformat='(A1)',col=fcolor
  axis, zaxis=0, x1, y1, z0, /t3d, ztickformat='(A1)',col=fcolor
  axis, zaxis=0, x1, y0, z0, /t3d, ztickformat='(A1)',col=fcolor

  plots, xx[*,0], xx[*,1], xx[*,2], psym=sym, col=pcolor, /t3d

  axis, zaxis=1, x0, y0, z0, /t3d, ztickformat='(A1)',col=fcolor
  axis, yaxis=0, x0, y0, z1, /t3d, ytickformat='(A1)',col=fcolor
  axis, xaxis=0, x0, y0, z1, /t3d, xtickformat='(A1)',col=fcolor
 
endif else if ( (nx ne 1) and (ny ne 1) and (nz eq 1) ) then begin
  
  plot, xx[*,0], xrange=[x0,x1], yrange=[y0,y1], /nodata, xstyle=1, ystyle=1
  plots, xx[*,0], xx[*,1], psym=sym
 
endif else if ( (nx ne 1) and (ny eq 1) and (nz ne 1) ) then begin
  
  plot, xx[*,0], xrange=[x0,x1], yrange=[z0,z1], /nodata, xstyle=1, ystyle=1
  plots, xx[*,0], xx[*,2], psym=sym
  if (drawgrid) then begin
    oplot, [x[nx/2]  ,x[nx/2]]  , [z[0]     ,z[nz-1]]
    oplot, [x[nx/2+1],x[nx/2+1]], [z[0]     ,z[nz-1]]
    oplot, [x[0]     ,x[nx-1]]  , [z[nz/2]  ,z[nz/2]]
    oplot, [x[0]     ,x[nx-1]]  , [z[nz/2+1],z[nz/2+1]]
  endif

endif else if ( (nx eq 1) and (ny ne 1) and (nz ne 1) ) then begin
  
  plot, xx[*,0], xrange=[y0,y1], yrange=[z0,z1], /nodata, xstyle=1, ystyle=1
  plots, xx[*,1], xx[*,2], psym=sym

endif

if (ps) then begin
  if (not quiet) then print, 'pc_plot_par: writing '+imgdir+'/'+filename
  device, /close
  set_plot, 'x'
endif

!p.charthick=pcharthick & !p.thick=pthick
!x.thick=xthick & !y.thick=ythick & !z.thick=zthick

end
