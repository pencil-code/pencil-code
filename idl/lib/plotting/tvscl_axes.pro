;;
;;  $Id: tvscl_axes.pro,v 1.1 2006-07-10 13:45:18 ajohan Exp $
;;
;;  tvscl_axes - tvscl contour plot with axes.
;;
;;  Author : Anders Johansen
;;  Date   : 10.07.2006
;;
pro tvscl_axes, ff, xax=xax, yax=yax, $
    min=min, max=max, $
    title=title, xtitle=xtitle, ytitle=ytitle, $
    zoom=zoom, xstyle=xstyle, ystyle=ystyle, noaxis=noaxis, position=position, $
    nonewwindow=nonewwindow, ps=ps, filename=filename
;;
;;  Set default values
;;
default, ps, 0
default, zoom, 1
default, xstyle, 0
default, ystyle, 0
default, noaxis, 0
default, filename, 'tvscl_axes.eps'
default, nonewwindow, 0
default, position, [0.1,0.1,0.9,0.9]
;;
;;  Quick way to get no axes.
;;
if (noaxis) then begin
  xstyle=4 & ystyle=4
endif
;;
;;  Resolution.
;;
nx=n_elements(ff[*,0])
ny=n_elements(reform(ff[0,*]))
;;
;;  Either axes are provided, otherwise they default to index arrays.
;;
if (n_elements(xax) eq 0) then xax=findgen(nx)
if (n_elements(yax) eq 0) then yax=findgen(ny)
x0=xax[0] & x1=xax[nx-1] & y0=yax[0] & y1=yax[ny-1]
;;
;;  Need 2-D axes for contour plot.
;;
xx=spread(xax,1,ny) & yy=spread(yax,0,nx)
;;
;;  The window is larger than the plot. Use the position keyword to specify
;;  where to place the plot.
;;
xsize=nx*1.0/(position[2]-position[0])
ysize=ny*1.0/(position[3]-position[1])
;;
;;  Treat min and max of the plot by replacing all values larger/smaller than
;;  min/max by the value min/max.
;;
if (defined(min)) then begin
  ii=where(ff gt min)
  if (ii[0] ne -1) then ff[ii]=min
endif
;;
if (defined(max)) then begin
  ii=where(ff gt max)
  if (ii[0] ne -1) then ff[ii]=max
endif
;;
;;  Plot either to postscript or to X. If one is going to see a movie series,
;;  it is only necessary to open a new window for the first frame, so one can
;;  set the keyword nonewwindow for all other frames than the first.
;;
if (ps) then begin
  set_plot, 'ps'
  device, filename=filename, /encapsulated, xsize=xsize, ysize=ysize, /pixels
endif else begin
  if (not nonewwindow) then window, retain=2, xsize=zoom*xsize, ysize=zoom*ysize
endelse
;;
;;  Plot with tvscl, then put axes on with contour.
;;
tvscl, rebin(ff,zoom*[nx,ny]), position[0]*zoom*xsize, position[1]*zoom*xsize
;;
contour, ff, xx, yy, /nodata, /noerase, position=position, $
    xstyle=xstyle, ystyle=ystyle, $
    title=title, xtitle=xtitle, ytitle=ytitle
;;
end
