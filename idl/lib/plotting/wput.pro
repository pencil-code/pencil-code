pro wput

;+
; WGET / WPUT
;
; uses a pixmap-window for "smooth" displaying a series of IDL-plots
; this avoids the "blinking"
;
; use wput to direct output to graphics window
; after plotting use wget to copy the pimap to the actual window
;
; SIDE EFFECTS:
;
; common block window_pixmap_copy
;
; EXAMPLE:
; 
; for i=0,100 DO $
;     BEGIN
;     wput              ; prepares the graphics output to go to a pixmap window
;     plot, data(*,i)   ; plot the data
;     wget              ; copy the contents of the pixmap to actual window
;     END
;
; Written: Hardi Peter, hpeter@linmpi.mpg.de
;
;-

common window_pixmap_copy, pix_win, act_win

;; only makes sense with X11
if (!d.name ne 'X') then return

if (!d.window lt 0) then window,0
;
act_win = !d.window
;
window, /FREE, /PIXMAP, XSIZE=!d.x_size, YSIZE=!d.y_size
pix_win = !d.window

return
end
   
