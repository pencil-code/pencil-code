pro func, t, C, f, der 

  om=C[1]
  expt=exp(2*C[0]*t)
  f=expt*(C[2]*cos(om*t)^2 + (C[3]*cos(om*t) + C[4]*sin(om*t))*sin(om*t))

  IF N_PARAMS() GE 4 THEN $
    der=[[2.*t*f],[t*expt*(sin(2*om*t)*(C[3]-C[2]) + cos(2*om*t)*C[4])],[expt*cos(om*t)^2],[expt*cos(om*t)*sin(om*t)],[expt*sin(om*t)^2]]
end

function lmfunc, t, C

  om=C[1]
  expt=exp(2*C[0]*t)
  f=expt*(C[2]*cos(om*t)^2 + (C[3]*cos(om*t) + C[4]*sin(om*t))*sin(om*t))
  return, [[f], [2.*t*f],[t*expt*(sin(2*om*t)*(C[3]-C[2]) + cos(2*om*t)*C[4])],[expt*cos(om*t)^2],[expt*cos(om*t)*sin(om*t)],[expt*sin(om*t)^2]]

end

pro fit_rms, t, rms, lam=lam, om=om, res=result
;
; Provides an estimate for growth rate and oscillation frequency of a field based on its rms. 
; Requires a rather good inital guess for the frequency (om), if the growth rate is not small compared to the frequency,
; also for the growth rate (lam).
;
  default, lam, 0.
  default, om, 0.

  C=[lam,om,1.,1.,1.]
  result = sqrt(CURVEFIT( t, rms^2, weights, C, /DOUBLE, FUNCTION_NAME='func', itmax=2500, tol=1e-9))
  ;C = LMFIT(t, rms^2, C, /DOUBLE, FUNCTION_NAME = 'lmfunc', tol=1e-9)

  print, 'lambda=', C[0]
  print, 'omega=', C[1]

end
