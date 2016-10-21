;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;   quintic_step.pro  ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: PABourdin
;;;  Date:   21. Oct 2016
;;;
;;;  Description:
;;;   Resamples the quintic_step function implemented in Pencil Code.

function quintic_step, x, x0, width

xi = (x-x0)/width
xi = xi > (-1.)
xi = xi < (+1.)
xi2 = xi^2
quintic_step = 0.5 + xi*(0.9375 + xi2*(-0.625 + xi2*0.1875))

return, quintic_step

end
