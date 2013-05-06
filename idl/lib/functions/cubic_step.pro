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

function cubic_step,x,x0,width

xi = (x-x0)/width
xi = xi > (-1.)
xi = xi < (+1.)

cubic_step = 0.5 + xi*(0.75-xi^2.*0.25)

return,cubic_step

end
