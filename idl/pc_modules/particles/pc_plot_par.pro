;
;  $Id: pc_plot_par.pro,v 1.2 2005-06-20 13:58:27 ajohan Exp $
;
pro pc_plot_par, xx, pos=pos, ps=ps, color=color, filename=filename

default, ps, 0
default, color, 0
default, filename, 'particles.eps'
default, pos, [0.1,0.1,0.9,0.9]

if (ps) then begin
  set_plot, 'ps'
  device, /encapsulated, color=color, xsize=14.0, ysize=11.0, $
      font_size=11, filename=filename
  !p.font=1
  !p.charsize=2.0
endif

if (color) then loadct, 12
frame_color=100
par_color=200

pc_read_param, obj=par, /quiet
pc_read_dim, obj=dim, /quiet

x0=par.xyz0[0]
x1=x0+par.Lxyz[0]
y0=par.xyz0[1]
y1=y0+par.Lxyz[1]
z0=par.xyz0[2]
z1=z0+par.Lxyz[2]

if ( (dim.nxgrid ne 1) and (dim.nygrid ne 1) and (dim.nzgrid ne 1) ) then begin

  surface,  [[0.0,0.0,0.0],[0.0,0.0,0.0]], col=frame_color, $
      xrange=[x0,x1], yrange=[y0,y1], zrange=[z0,z1], $
      xstyle=1, ystyle=1, zstyle=1, /save, /nodata, pos=pos
  axis, xaxis=1, 0.0, y1, z1, /t3d, xtickformat='noticknames_aj',col=frame_color
  axis, xaxis=1, 0.0, y1, z0, /t3d, xtickformat='noticknames_aj',col=frame_color
  axis, yaxis=1, x1, 0.0, z0, /t3d, ytickformat='noticknames_aj',col=frame_color
  axis, yaxis=1, x1, 0.0, z1, /t3d, ytickformat='noticknames_aj',col=frame_color
  axis, zaxis=1, x1, y1, 0.0, /t3d, ztickformat='noticknames_aj',col=frame_color
  axis, zaxis=1, x1, y0, 0.0, /t3d, ztickformat='noticknames_aj',col=frame_color

  plots, xx[*,0], xx[*,1], xx[*,2], psym=3, col=par_color, /t3d

  axis, zaxis=1, x0, y0, 0.0, /t3d, ztickformat='noticknames_aj',col=frame_color
  axis, yaxis=1, x0, 0.0, z1, /t3d, ytickformat='noticknames_aj',col=frame_color
  axis, xaxis=1, 0.0, y0, z1, /t3d, xtickformat='noticknames_aj',col=frame_color
 
endif else if ( (dim.nxgrid ne 1) and (dim.nygrid eq 1) and (dim.nzgrid ne 1) ) then begin
  
  plot, xx[*,0], xrange=[x0,x1], yrange=[z0,z1], /nodata
  plots, xx[*,0], xx[*,2], psym=3

endif else if ( (dim.nxgrid eq 1) and (dim.nygrid ne 1) and (dim.nzgrid ne 1) ) then begin
  
  plot, xx[*,0], xrange=[y0,y1], yrange=[z0,z1], /nodata
  plots, xx[*,1], xx[*,2], psym=3

endif

if (ps) then begin
  device, /close
  print, 'Written file '+filename
  set_plot, 'x'
endif

end
