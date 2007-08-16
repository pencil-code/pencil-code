;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   plot_binned.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  $Date: 2007-08-16 23:22:43 $
;;;  $Revision: 1.9 $
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
;;;      /FILL              --- /CONTOUR, FILL=0 does contour lines
;;;      NLEVELS=nlevels    --- set NLEVELS for /CONTOUR (default=60)
;;;      /PLOTIMAGE         --- use PLOTIMAGE to display density
;;;                             (faster than CONTOUR and PS files are
;;;                             much smaller)
;;;      /INVERT            --- invert data for /CONTOUR (default is
;;;                             white for no data, dark for dense data
;;;                             with standard colour maps)
;;;      /OVERPLOT          --- do overplotting (kind of `oplot_binned')
;;;      /EQUX, /EQUY       --- equalize in x/y, i.e. normalize
;;;                             vertical/horizontal sum of counts to
;;;                             unity (except where there is no
;;;                             point). Useful with /PLOTIMAGE or
;;;                             /CONTOUR when plotting over radius
;;;                             data in cubic box
;;;      ENHANCE=enhance    --- multiply density map with ENHANCE,
;;;                             which must be a non-negative 2d array
;;;                             of size NX x NY
;;;      ASINHSCALE=scale   --- plot asinh(scale*density) instead of density
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
;;;    plot_binned, xx, ff, NX=60, NY=40, /PLOTIMAGE, /EQUX
;;;    plot_binned, xx, ff, /PLOTIMAGE, /EQUX, ASINHSCALE=5.
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
                 PSYM=psym, SYMSIZE=symsize, $
                 COLOR=color, $
                 XRANGE=xrange, YRANGE=yrange, $
                 XSTYLE=xstyle, YSTYLE=ystyle, $
                 EQUX=equx, EQUY=equy, $
                 ENHANCE=enhance, $
                 ASINHSCALE=asinhscale, $
                 DENSITY=map, XCOORDS=mapx, YCOORDS=mapy, $
                 HELP=help, $
                 _EXTRA=extra

if (keyword_set(help)) then extract_help, 'plot_binned'

if (n_elements(nmx)        eq 0) then nmx = 500/(!p.multi[1] > 1)
if (n_elements(nmy)        eq 0) then nmy = 350/(!p.multi[2] > 1)

if (n_elements(psym)       eq 0) then psym       = 3
if (n_elements(contour)    eq 0) then contour    = 0
if (n_elements(plotimage)  eq 0) then plotimage  = 0
if (n_elements(invert)     eq 0) then invert     = 0
if (n_elements(fill)       eq 0) then fill       = 1
if (n_elements(nlevels)    eq 0) then nlevels    = 60
if (n_elements(equx)       eq 0) then equx       = 0
if (n_elements(equy)       eq 0) then equy       = 0
if (n_elements(enhance)    eq 0) then enhance    = -1
if (n_elements(asinhscale) eq 0) then asinhscale = 0
if (n_elements(xrange)     eq 0) then xrange     = [0,0]
if (n_elements(yrange)     eq 0) then yrange     = [0,0]
if (n_elements(xstyle)     eq 0) then xstyle     = !x.style
if (n_elements(ystyle)     eq 0) then ystyle     = !y.style

;; Sanitize keywords
if (plotimage ne 0) then contour=0
if ((plotimage ne 0) or (contour ne 0)) then scatter=0 else scatter=1 
if ((equx ne 0) and (equy ne 0)) then $
    message, '/EQUX and /EQUY keywords are mutualy exclusive'

;; Save so we can reset later to avoid plotting to a new frame
pmulti1 = !p.multi

;; Plot the frame (necessary to get x- and y-range)
if (n_elements(overplot) eq 0) then begin
  if (scatter eq 0) then begin  ; don't plot and avoid extra PS frame
    xstyle1 = (xstyle or 4)
    ystyle1 = (ystyle or 4)
  endif else begin
    xstyle1 = xstyle
    ystyle1 = ystyle
  endelse
  plot, xvar, yvar, /NODATA, $
      XRANGE=xrange, YRANGE=yrange, $
      XSTYLE=xstyle1, YSTYLE=ystyle1, $
      COLOR=color, _EXTRA=extra
  pmulti2 = !p.multi       ; save so we can exit with correct !p.multi
endif

xr = !x.crange                  ; NB: for /xlog, this is alog10(data_range)
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
  imgyrange = yr
endif else begin
  yy=alog10(yvar)
  mapy = logspace(yr,Nmy)
  ybin = 10^(yr[0] + (findgen(Nmy)+0.5)*dmy)
  imgyrange = 10^yr
endelse

map = hist_2d(xx, yy, $
              MIN1=xr[0], MAX1=xr[1], BIN1=dmx, $
              MIN2=yr[0], MAX2=yr[1], BIN2=dmy) * 1.D0
;; hist and hist_2d do really bizarre things to the last bin. Here
;; this just means we should remove the last entry in each
;; direction:
map = map[0:NMx-1,0:NMy-1]

;; Equalize if requested
if (equx ne 0) then begin
  for ix=0,NMx-1 do begin
    sum = total(1.D0*abs(map[ix,*]))
    if (sum gt 0) then map[ix,*] = map[ix,*]/sum
  endfor
endif
;
if (equy ne 0) then begin
  for iy=0,NMy-1 do begin
    sum = total(1.D0*abs(map[*,iy]))
    if (sum gt 0) then map[*,iy] = map[*,iy]/sum
  endfor
endif

;; Multiply with ENHANCE, if requested
if (enhance[0] ge 0) then begin
  s = size(enhance)
  if ((s[0] ne 2) or (s[1] ne NMx) or (s[2] ne NMy)) then $
      message, 'ENHANCE must be of size '+strtrim(NMx,2)+'x'+strtrim(NMy,2)
  map = map*enhance
endif

;; Apply asinh if requested
if (asinhscale ne 0) then begin
  map = map / max(abs(map))
  map = asinh(asinhscale*map)
endif


if (scatter eq 0) then !p.multi = pmulti1

;; Plot
if (scatter ne 0) then begin    ; scatter plot
  ;; Re-use xx, yy here to save memory
  ;; 2-d x- and y-coordinates for the map
  xx = spread(indgen(NMx),[1],[NMy])
  yy = spread(indgen(NMy),[0],[NMx])
  goodm = where(map ne 0)
  ;
  if (goodm[0] ge 0) then begin
    oplot, xbin[xx[goodm]], ybin[yy[goodm]], $
           PSYM=psym, SYMSIZE=symsize, COLOR=color, $
           _EXTRA=extra
  endif else begin
    message, /INFO, 'No valid points in domain' $
             + '[' + wdtrim(xr[0]) + ', '+ wdtrim(xr[1]) +']'$
             + ' x ' $
             + '[' + wdtrim(yr[0]) + ', '+ wdtrim(yr[1]) +']'
  endelse

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
      plot, xvar, yvar, /NODATA, COLOR=color, /NOERASE, $
          XRANGE=xr, YRANGE=yr, $
          XSTYLE=((xstyle or 1) and 13), $ ; set 1, unset 2
          YSTYLE=((ystyle or 1) and 13), $
          _EXTRA=extra
      !p.multi = pmulti2
    endif
  endif else begin              ; plotimage
    plotimage, fact*map, RANGE=minmax(fact*map), $
        IMGXRANGE=imgxrange, IMGYRANGE=imgyrange, $
          XSTYLE=((xstyle or 1) and 13), $ ; set 1, unset 2
          YSTYLE=((ystyle or 1) and 13), $
        _EXTRA=extra
  endelse
endelse

end
; End of file plot_binned.pro
