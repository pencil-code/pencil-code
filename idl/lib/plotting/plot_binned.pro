;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   plot_binned.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  $Date: 2004-09-15 02:37:47 $
;;;  $Revision: 1.2 $
;;;  Description:
;;;    Scatter-plot data, but bin them first in order to reduce the
;;;    size of thusly created PostScript files. For 60^3 data points,
;;;    this can easily result in a 5 times smaller file.
;;;      Alternatively, the density of the binned data can be displayed
;;;    using PLOTIMAGE (recommended) or CONTOUR 
;;;
;;;  Usage: like `plot', but with the additional key words
;;;      NX=nx, NY=ny       --- set vertical and horizontal number of points
;;;      /CONTOUR           --- do contour plot (by default filled)
;;;      FILL=fill          --- /CONTOUR, FILL=0 does contour lines
;;;      NLEVELS=nlevels    --- set NLEVELS for /CONTOUR (default=60)
;;;      /PLOTIMAGE         --- use PLOTIMAGE to display density
;;;                             (faster than CONTOUR and PS files are
;;;                             much smaller)
;;;      /INVERT            --- invert data for /CONTOUR (default is
;;;                             white for no data, dark for dense data
;;;                             with standard colour maps)
;;;      OVERPLOT=overplot  --- do overplotting (kind of `oplot_binned')
;;;      DENSITY=density    --- return density map (only values 0 an 1
;;;                             unless /CONTOUR was sepcified
;;;      XCOORDS=xcoords    --- corresponding x coordinates
;;;      YCOORDS=ycoords    --- corresponding x coordinates
;;;
;;;    Some key words like PSYM, SYMSIZE and COLOR are intecepted, such as
;;;    to work in overplot mode as well as in normal plot mode. For the
;;;    rest, draw the frame first explicitly and then call plot_binned
;;;    with /OVERPLOT.
;;;
;;;  Examples:
;;;    N=8000 & xx=1+5.*randomu(seed,N) & ff=2.5+cos(4*xx)+0.2*randomn(seed,N)
;;;    plot_binned, xx, ff
;;;    plot_binned, xx, ff, NX=60, NY=40, /PLOTIMAGE, /XLOG, /YLOG
;;;    plot_binned, xx, ff, NX=60, NY=40, /CONTOUR
;;;    plot_binned, xx, ff, NX=60, NY=40, /CONTOUR, /INVERT
;;;    
;;;  To do:
;;;    Reduce memory usage: xx,yy,good,xind,yind all use the same
;;;    amount of space as the original data xvar, yvar
;;;    [Partly accomplished]

pro plot_binned, xvar, yvar, $
                 NX=nmx, NY=nmy, $
                 CONTOUR=contour, $
                 FILL=fill, NLEVELS=nlevels, $
                 PLOTIMAGE=plotimage, $
                 INVERT=invert, $
                 OVERPLOT=overplot, $
                 PSYM=psym, $
                 SYMSIZE=symsize, $
                 COLOR=color, $
                 DENSITY=map, XCOORDS=mapx, YCOORDS=mapy, $
                 _EXTRA=extra

if (n_elements(nmx)       eq 0) then nmx = 500/(!p.multi[1] > 1)
if (n_elements(nmy)       eq 0) then nmy = 350/(!p.multi[2] > 1)

if (n_elements(psym)      eq 0) then psym      = 3
if (n_elements(contour)   eq 0) then contour   = 0
if (n_elements(plotimage) eq 0) then plotimage = 0
if (n_elements(invert)    eq 0) then invert    = 0
if (n_elements(fill)      eq 0) then fill      = 1
if (n_elements(nlevels)   eq 0) then nlevels   = 60
if (n_elements(color)     eq 0) then color     = !p.color

;; Sanitize keywords
if (plotimage ne 0) then contour=0
if ((plotimage ne 0) or (contour ne 0)) then scatter=0 else scatter=1 

;; Save so we can reset later to avoid plotting to a new frame
pmulti1 = !p.multi

;; Plot the frame (necessary to get x- and y-range)
if (n_elements(overplot) eq 0) then begin
  if (scatter eq 0) then begin  ; don't plot, as this will be overplotted later
    xstyle=4 & ystyle=4
  endif else begin              ; really plot frame; this will remain
    xstyle=!x.style & ystyle=!y.style
  endelse
  plot, xvar, yvar, /NODATA, $
      XSTYLE=xstyle, YSTYLE=ystyle, $
      COLOR=color, _EXTRA=extra
  pmulti2 = !p.multi       ; save so we can exit with correct !p.multi
endif

xr = !x.crange                  ; >B: for /xlog, this is alog10(data_range)
yr = !y.crange

dmx = (xr[1]-xr[0])*1./Nmx
dmy = (yr[1]-yr[0])*1./Nmy

;; Handle /xlog and /ylog
if (!x.type ne 1) then begin
  xx=xvar
  mapx = linspace(xr,Nmx)
  xbin = xr[0] + (findgen(Nmx)+0.5)*dmx
  imgxrange = xr
endif else begin
  xx=alog10(xvar)
  mapx = logspace(xr,Nmx)
  xbin = 10^(xr[0] + (findgen(Nmx)+0.5)*dmx)
  imgxrange = 10^xr
endelse
;
if (!y.type ne 1) then begin
  yy=yvar
  mapy = linspace(yr,Nmy)
  ybin = yr[0] + (findgen(Nmy)+0.5)*dmy
  imgyrange = xr
endif else begin
  yy=alog10(yvar)
  mapy = logspace(yr,Nmy)
  ybin = 10^(yr[0] + (findgen(Nmy)+0.5)*dmy)
  imgyrange = 10^yr
endelse

map = hist_2d(xx, yy, $
              MIN1=xr[0], MAX1=xr[1], BIN1=dmx, $
              MIN2=yr[0], MAX2=yr[1], BIN2=dmy)
;; hist and hist_2d do really bizarre things to the last bin. Here
;; this just means we should remove the last entry in each
;; direction:
map = map[0:NMx-1,0:NMy-1]

if (scatter eq 0) then !p.multi = pmulti1

;; Plot
if (scatter ne 0) then begin    ; scatter plot
  ;; Re-use xx, yy here to save memory
  ;; 2-d x- and y-coordinates for the map
  xx = spread(indgen(NMx),[1],[NMy])
  yy = spread(indgen(NMy),[0],[NMx])
  goodm = where(map ne 0)
  ;
  oplot, xbin[xx[goodm]], ybin[yy[goodm]], $
      PSYM=psym, SYMSIZE=symsize, COLOR=color, $
      _EXTRA=extra

endif else begin
  if (invert ne 0) then begin
    levels = linspace(0.,max(map),NLEVELS,GHOST=0.5)
    fact = 1.
  endif else begin
    levels = linspace(-max(map),0,NLEVELS,GHOST=0.5)
    fact = -1.
  endelse
  if (contour ne 0) then begin  ; contour plot
    contour, fact*map, mapx, mapy, $
        /OVERPLOT, $
        FILL=fill, LEVELS=levels, $
        _EXTRA=extra
    ;; Replot the frame as CONTOUR will have erased most of it
    if (n_elements(overplot) eq 0) then begin
      plot, xvar, yvar, /NODATA, COLOR=color, /NOERASE, _EXTRA=extra
      !p.multi = pmulti2
    endif
  endif else begin              ; plotimage
    plotimage, fact*map, RANGE=minmax(fact*map), $
        IMGXRANGE=imgxrange, IMGYRANGE=imgyrange, _EXTRA=extra
  endelse
endelse

end
; End of file plot_binned.pro
