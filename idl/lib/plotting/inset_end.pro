;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   inset_end.pro   ;;;
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
;;;    plot, [0,1]
;;;    inset_start, [0.6,0.9], [0.1,0.4]
;;;    plot, [1,0]
;;;    inset_end
;;;  where [xll,yll] is the lower left point, [xur,yur] the upper
;;;  right point of the inset data window -- expressed in normalised
;;;  coordinates of the original plot.

pro inset_end

  common _inset, excl_p, excl_x, excl_y

  !p = excl_p                   ; Restore !p
  !x = excl_x                   ; Restore !x
  !y = excl_y                   ; Restore !y

end

; End of file inset_end.pro
