;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   plot_binned.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  $Date: 2004-05-30 12:55:47 $
;;;  $Revision: 1.1 $
;;;  Description:
;;;   Scatter-plot data, but bin them first in order to reduce the
;;;   size of thusly created PostScript files. For 60^3 data points,
;;;   this can easily result in a 5 times smaller file.
;;;
;;;   Usage: like `plot', but with the additional key words
;;;     NX=nx, NY=ny       --- set vertical and horizontal number of points
;;;     OVERPLOT=overplot  --- do overplotting (kind of `oplot_binned')
;;;
;;;   Some key words like PSYM, SYMSIZE and COLOR are trapped, such as
;;;   to work in overplot mode as well as in normal plot mode. For the
;;;   rest, draw the frame first explicitly and then call plot_binned
;;;   with /OVER.
;;;
;;;  To do:
;;;   Reduce memory usage: xx,yy,good,xind,yind all use the same
;;;   amount of space as the original data xvar, yvar

pro plot_binned, xvar, yvar, $
                 NX=nmx, NY=nmy, $
                 OVERPLOT=overplot, $
                 PSYM=psym, $
                 SYMSIZE=symsize, $
                 COLOR=color, $
                 _EXTRA=extra

if (n_elements(nmx) eq 0) then nmx = 500/(!p.multi[1] > 1)
if (n_elements(nmy) eq 0) then nmy = 350/(!p.multi[2] > 1)

if (n_elements(psym) eq 0) then psym = 3

if (n_elements(color) eq 0) then color = !p.color


;; Plot the frame
if (n_elements(overplot) eq 0) then begin
  plot, xvar, yvar, PSYM=3, /NODATA, COLOR=color, _EXTRA=extra
endif

xr = !x.crange
yr = !y.crange

dmx = (xr[1]-xr[0])*1./Nmx
dmy = (yr[1]-yr[0])*1./Nmy
map = bytarr(NMx,NMy)

;; Handle /xlog and /ylog
if (!x.type ne 1) then begin
  xx=xvar
  xbin = xr[0] + (findgen(Nmx)+0.5)*dmx
endif else begin
  xx=alog10(xvar)
  xbin = 10^(xr[0] + (findgen(Nmx)+0.5)*dmx)
endelse
;
if (!y.type ne 1) then begin
  yy=yvar
  ybin = yr[0] + (findgen(Nmy)+0.5)*dmy
endif else begin
  yy=alog10(yvar)
  ybin = 10^(yr[0] + (findgen(Nmy)+0.5)*dmy)
endelse

good = where((xx ge xr[0]) and (xx le xr[1]))

xind = (xx[good]-xr[0]) / dmx
yind = (yy[good]-yr[0]) / dmy

map[xind,yind] = 1

;loadct,21
;tvscl,congridscale(map,5),11
;loadct,5

;; 2-d x- and y-coordinates for the map
mapxx = rebin(reform(findgen(Nmx), Nmx,1), Nmx,Nmy)
mapyy = rebin(reform(findgen(Nmy), 1,Nmy), Nmx,Nmy)
goodm = where(map ne 0)

if (n_elements(overplot) eq 0) then begin
  oplot, xbin[mapxx[goodm]], ybin[mapyy[goodm]], $
      PSYM=psym, SYMSIZE=symsize, COLOR=color
endif else begin
  oplot, xbin[mapxx[goodm]], ybin[mapyy[goodm]], $
      PSYM=psym, SYMSIZE=symsize, COLOR=color, $
      _EXTRA=extra 
endelse

end
; End of file plot_binned.pro
