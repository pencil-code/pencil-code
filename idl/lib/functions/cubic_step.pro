;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   cubic_step.pro  ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: bing (mail@bingert.com)
;;;  Date:   06. mai 2013
;;;
;;;  Description:
;;;   Resamples the cubic_step function implemented
;;;   in the fortran code. The calling sequence is the same
;;;
;;;  Parameters:
;;;   x     : coordinates
;;;   x0    : middle coordinate of the transition
;;;   width : width is the half width of the transition
;;;
;;;  Result:
;;;   provides a smooth transition from 0 to 1 between
;;;   the coordinates x0-width/2 (=0) to x0+width/2 (=1)

function cubic_step,x,x0,width,shift

xi = (x-x0)/width - shift
xi = xi > (-1.)
xi = xi < (+1.)

cubic_step = 0.5 + xi*(0.75-xi^2.*0.25)

return,cubic_step

end
