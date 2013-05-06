;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   cubic_der_step.pro  ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: bing (mail@bingert.com)
;;;  Date:   06. mai 2013
;;;
;;;  Description:
;;;   Resamples the cubic_der_step function implemented
;;;   in the fortran code. The calling sequence is the same.
;;;   This routine returns the first derivative of the
;;;   cubic_step function

function cubic_der_step,x,x0,width

xi = (x-x0)/width
xi = xi > (-1.)
xi = xi < (+1.)

cubic_der_step = 0.75*(1. - xi^2)/width

return,cubic_der_step

end
