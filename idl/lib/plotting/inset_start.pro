;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   inset_start.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@ncl.ac.uk)
;;;  Date:   12-Feb-2001
;;;  Version: 0.11
;;;
;;;  Description:
;;;   Set up plotting for inset into the previously plotted figure.
;;;  Usage
;;;      inset_start, xll, yll, xur, yur
;;;      inset_start, [xll, xur], [yll, yur]
;;;    where [xll,yll] is the lower left point, [xur,yur] the upper
;;;    right point of the inset data window -- expressed in normalised
;;;    coordinates of the original plot.
;;;  Example
;;;     plot, [0,1]
;;;     inset_start, [0.6,0.9], [0.1,0.4]
;;;     plot, [1,0]
;;;     inset_end
;;;   where [xll,yll] is the lower left point, [xur,yur] the upper
;;;   right point of the inset data window -- expressed in normalised
;;;   coordinates of the original plot.
;;;  To do:
;;;   Turn _inset into a stack, so insets within insets are possible
;;;   (does anybody really need this?).
;;;

pro inset_start, xll, yll, xur, yur

  common _inset, excl_p, excl_x, excl_y

  ;; Map the two possible input conventions onto a unique state
  if (n_elements(xll) gt 1) then begin ; Vector form
    xur = xll[1] & yur = yll[1]
    xll = xll[0] & yll = yll[0]
  endif

  xr = !x.window
  yr = !y.window
  dxr = xr[1]-xr[0]
  dyr = yr[1]-yr[0]
  xr1 = xr[0]+dxr*[xll,xur]
  yr1 = yr[0]+dyr*[yll,yur]

  excl_p = !p                   ; Save !p
  excl_x = !x                   ; Save !x
  excl_y = !y                   ; Save !y

  !p.multi[0] = !p.multi[0]+1   ; Keep any previous plots
  !p.position = [xr1[0],yr1[0], xr1[1],yr1[1]] ; Position window

  ;; 20% smaller characters
  if (!p.charsize eq 0 ) then !p.charsize=1
  !p.charsize = 0.8*!p.charsize

end
; End of file inset_start.pro
