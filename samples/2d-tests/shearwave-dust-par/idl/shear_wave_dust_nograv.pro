;;
;;  $Id$
;;
;;  Semi-analytical solution of the linear growth of a non-axisymmetric
;;  dust particle wave in a Keplerian shear flow.
;;
;;  Usage: 
;;    Probably to plot rrhod versus tt and compare with the result produced
;;    by the Pencil Code.
;;
ii=complex(0,1)
;;  Physical parameters.
rho0 =1.0d
rhod0=1.0d
cs2  =1.0
eps0=rhod0/rho0
kx0=0.0d & kx=kx0 & ky =1.0d
Omega =1.0d
qshear=1.5d
tauf  =1.0d
;;  Parameters for the time-solver.
nt=1000
dt=1.0d-3
tt=dblarr(nt)
;;  Runge-Kutta coefficients
alpha= [  0.0d, -5/ 9.0d, -153/128.0d]
beta = [1/3.0d, 15/16.0d,    8/ 15.0d]
dt_beta=dt*beta
;;  The result is stored in these arrays.
uux  =dcomplexarr(nt)
uuy  =dcomplexarr(nt)
rrho =dcomplexarr(nt)
wwx  =dcomplexarr(nt)
wwy  =dcomplexarr(nt)
rrhod=dcomplexarr(nt)
;;  Initial condition
t   =0.0d
ux  =complex(0.0d  ,0.0d)
uy  =complex(0.0d  ,0.0d)
rho =complex(0.0d  ,0.0d)
wx  =complex(1.0d-3,0.0d)
wy  =complex(0.0d  ,0.0d)
rhod=complex(0.0d  ,0.0d)
;;
;;  Time integration of the linearized equations.
;;
for it=0L,nt-1 do begin

  tt[it]=t
  uux[it]  =abs(ux)
  uuy[it]  =abs(uy)
  rrho[it] =abs(rho)
  wwx[it]  =abs(wx)
  wwy[it]  =abs(wy)
  rrhod[it]=abs(rhod)

  for itsub=0,n_elements(beta)-1 do begin
    if (itsub eq 0) then begin
      duxdt  =complex(0.0d,0.0d)
      duydt  =complex(0.0d,0.0d)
      drhodt =complex(0.0d,0.0d)
      dwxdt  =complex(0.0d,0.0d)
      dwydt  =complex(0.0d,0.0d)
      drhoddt=complex(0.0d,0.0d)
      ds=0
    endif else begin
      duxdt  =alpha[itsub]*duxdt
      duydt  =alpha[itsub]*duydt
      drhodt =alpha[itsub]*drhodt
      dwxdt  =alpha[itsub]*dwxdt
      dwydt  =alpha[itsub]*dwydt
      drhoddt=alpha[itsub]*drhoddt
      ds=alpha[itsub]*ds
    endelse

    kx=kx0+qshear*Omega*t*ky
;;  Add Coriolis force.
    duxdt =  duxdt +            2.0d*Omega*uy
    duydt =  duydt - (2.0d - qshear)*Omega*ux
    dwxdt =  dwxdt +            2.0d*Omega*wy
    dwydt =  dwydt - (2.0d - qshear)*Omega*wx
;;  Add pressure force.
    duxdt =  duxdt - 1/rho0*cs2*ii*kx*rho
    duydt =  duydt - 1/rho0*cs2*ii*ky*rho
;;  Add drag force.
    duxdt =  duxdt - eps0/tauf*(ux-wx)
    duydt =  duydt - eps0/tauf*(uy-wy)
    dwxdt =  dwxdt - 1/tauf*(wx-ux)
    dwydt =  dwydt - 1/tauf*(wy-uy)
;;  Continuity equation.
    drhodt  = drhodt  - rho0 *ii*(kx*ux+ky*uy)
    drhoddt = drhoddt - rhod0*ii*(kx*wx+ky*wy)

    ds=ds+1

    ux   = ux   + duxdt  *dt_beta[itsub]
    uy   = uy   + duydt  *dt_beta[itsub]
    rho  = rho  + drhodt *dt_beta[itsub]
    wx   = wx   + dwxdt  *dt_beta[itsub]
    wy   = wy   + dwydt  *dt_beta[itsub]
    rhod = rhod + drhoddt*dt_beta[itsub]
    t = t + ds  *dt_beta[itsub]
  endfor

endfor


end
