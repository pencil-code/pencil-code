;;;;;;;;;;;;;;;;;;;;;;;;
;;;   fllevels.pro   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;

;;;
;;;  Author: wd (Wolfgang.Dobler@kis.uni-freiburg.de)
;;;  Date:   08-Apr-2002
;;;
;;;  Description:
;;;   Construct levels for plotting field lines in axisymmetry such
;;;   that a uniform vertical field will get equidistant lines.
;;;  Usage:
;;;   contour, xx*Aphi, x, z, LEVELS=fllevels(xx*Aphi)
;;;

function fllevels, var, nlev, $
                  GHOST=ghost

  default, ghost, 0.5           ; pretend data fill slighlty larger interval
  default, nlev, 20             ; number of field lines

  vmin = min(var) & vmax = max(var)
  smin = sign(vmin) & smax = sign(vmax)

  Nmin = nlev*smin*sqrt(abs(vmin))/(smax*sqrt(abs(vmax))-smin*sqrt(abs(vmin)))
; Nmax = nlev*smax*sqrt(abs(vmax))/(smax*sqrt(abs(vmax))-smin*sqrt(abs(vmin)))
;;  ii = Nmin + indgen(nlev+1)
  ii = Nmin + indgen(nlev+1-2*ghost) + ghost

  lev = ii*abs(ii)/nlev^2 $
        * (abs(vmax) + abs(vmin) - smin*smax*2.*sqrt(abs(vmin*vmax)))

  return, lev

end
; End of file fllevels.pro
