;;;;;;;;;;;;;;;;;;;
;;;   pbb.pro   ;;;
;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   13-May-2002
;;;
;;;  Description:
;;;   Plot isocontours of magnetic field strength

bb = curl(aa)
b2 = dot2(bb)

;; Add a small sphere to mark the origin
;b2 = b2 + max(b2)*exp(-0.5*(xx^2+yy^2+zz^2)/0.2^2)

default, winsize, [450,450]

save_state

loadct,0

s = size(b2)
xr = [0,s[1]]
yr = [0,s[2]]
zr = [0,s[3]]

;; The following is adapted from D. Fanning's web site
;; (http://www.dfanning.com/tips/volume_axes.html ``Adding Axes to a
;; Rendered 3D Volume'' 

orig_device = !d.name
set_plot, 'Z'
erase
device, SET_RESOLUTION=winsize, SET_COLORS=!d.n_colors

;; Establish 3d coordinate system
;scale3, XRANGE=xr, YRANGE=yr,ZRANGE=zr, AX=90, AZ=0
surface, fltarr(2,2), /NODATA, $
    XRANGE=xr, YRANGE=yr,ZRANGE=zr, $
    XTITLE='!8x!X', YTITLE='!8y!X', ZTITLE='!8z!X', $
    AX=35, AZ=25, $
    /NOERASE, /SAVE

shade_volume, b2, 0.3*max(b2), vertices, polygon, /LOW
image = polyshade(vertices, polygon, /T3D, XSIZE=winsize[0], YSIZE=winsize[1])

;; Take a snapshot of the display window

snapshot = tvrd()

set_plot, orig_device
tv, snapshot

; axis, 0,0,0, /XAXIS, XRANGE=[0,32], XTITLE='!8x!X', /T3D,/NOERASE
; axis, 0,0,0, /YAXIS, YRANGE=[0,32], YTITLE='!8y!X', /T3D,/NOERASE
; axis, 0,0,0, /ZAXIS, ZRANGE=[0,32], ZTITLE='!8z!X', /T3D,/NOERASE

restore_state

end
; End of file pbb.pro
