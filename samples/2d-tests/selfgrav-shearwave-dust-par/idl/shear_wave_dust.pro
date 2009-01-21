;;
;;  $Id$
;;
;;  Semi-analytical solution of the linear growth of a non-axisymmetric
;;  self-gravitating dust wave in a Keplerian shear flow.
;;
;;  Usage: 
;;    Probably to plot rrho versus tt and compare with the result produced
;;    by the Pencil Code.
;;
ii=complex(0,1)
;;  Physical parameters.
G=1.0d
rhod0=1.0d
kx0=-1.0d & kx=kx0
ky=1.0d
Omega=1.0d
qshear=1.5d
tauf=1.0d
;;  Parameters for the time-solver.
nt=10000
dt=1.0d-3
tt=dblarr(nt)
;;  Runge-Kutta coefficients
alpha= [  0.0d, -5/ 9.0d, -153/128.0d]
beta = [1/3.0d, 15/16.0d,    8/ 15.0d]
dt_beta=dt*beta
;;  The result is stored in these arrays.
wwx  =dblarr(nt)
wwy  =dblarr(nt)
rrhod=dblarr(nt)
;;  Initial condition
t  =0.0d
rhod=complex(0.001d,0.0d)
wx  =complex(0.0d  ,0.0d)
wy  =complex(0.0d  ,0.0d)
;;
;;  Time integration of the linearized equations.
;;
for it=0L,nt-1 do begin

  tt[it]=t
  wwx[it]  =abs(wx)
  wwy[it]  =abs(wy)
  rrhod[it]=abs(rhod)

  for itsub=0,n_elements(beta)-1 do begin
    if (itsub eq 0) then begin
      dwxdt  =complex(0.0d,0.0d)
      dwydt  =complex(0.0d,0.0d)
      drhoddt=complex(0.0d,0.0d)
      ds=0
    endif else begin
      dwxdt  =alpha[itsub]*dwxdt
      dwydt  =alpha[itsub]*dwydt
      drhoddt=alpha[itsub]*drhoddt
      ds=alpha[itsub]*ds
    endelse

    kx=kx0+qshear*Omega*t*ky
    phi = -4*!dpi*G*rhod/(kx^2+ky^2)
;;  Add Coriolis force.
    dwxdt =  dwxdt +            2.0d*Omega*wy
    dwydt =  dwydt - (2.0d - qshear)*Omega*wx
;;  Add gravitational acceleration.
    dwxdt =  dwxdt - ii*kx*phi
    dwydt =  dwydt - ii*ky*phi
;;  Add drag force.
    dwxdt =  dwxdt - 1/tauf*wx
    dwydt =  dwydt - 1/tauf*wy
;;  Continuity equation.
    drhoddt = drhoddt - rhod0*ii*(kx*wx+ky*wy)

    ds=ds+1

    wx   = wx   + dwxdt  *dt_beta[itsub]
    wy   = wy   + dwydt  *dt_beta[itsub]
    rhod = rhod + drhoddt*dt_beta[itsub]
    t = t + ds  *dt_beta[itsub]
  endfor

endfor


end
