;;
;;  $Id$
;;
;;  Semi-analytical solution of the linear growth of a non-axisymmetric
;;  self-gravitating wave in a Keplerian shear flow.
;;
;;  Usage: 
;;    Probably to plot rrho versus tt and compare with the result produced
;;    by the Pencil Code.
;;
ii=complex(0,1)
;;  Physical parameters. Viscosity is added to make comparison with the
;;  code easier.
G=0.1d
rho0=1.0d
cs=1.0d
kx0=-1.0d & kx=kx0
ky=1.0d
Omega=1.0d
qshear=1.5d
nu=1.0d-3
;;  Parameters for the time-solver.
nt=5000
dt=1.0d-3
tt=dblarr(nt)
;;  Runge-Kutta coefficients
alpha= [  0.0d, -5/ 9.0d, -153/128.0d]
beta = [1/3.0d, 15/16.0d,    8/ 15.0d]
dt_beta=dt*beta
;;  The result is stored in these arrays.
uux =dcomplexarr(nt)
uuy =dcomplexarr(nt)
rrho=dcomplexarr(nt)
;;  Initial condition
t  =0.0d
rho=complex(0.001d,0.0d)
ux =complex(0.0d  ,0.0d)
uy =complex(0.001d,0.0d)
;;
;;  Time integration of the linearized equations.
;;
for it=0L,nt-1 do begin

  tt[it]=t
  uux[it]=abs(ux)
  uuy[it]=abs(uy)
  rrho[it]=abs(rho)

  for itsub=0,n_elements(beta)-1 do begin
    if (itsub eq 0) then begin
      duxdt =complex(0.0d,0.0d)
      duydt =complex(0.0d,0.0d)
      drhodt=complex(0.0d,0.0d)
      ds=0
    endif else begin
      duxdt =alpha[itsub]*duxdt
      duydt =alpha[itsub]*duydt
      drhodt=alpha[itsub]*drhodt
      ds=alpha[itsub]*ds
    endelse

    kx=kx0+qshear*Omega*t*ky
    phi = -4*!dpi*G*rho/(kx^2+ky^2)
;;  Add Coriolis force and pressure gradient force.
    duxdt =  duxdt +            2.0d*Omega*uy-cs^2*1/rho0*ii*kx*rho
    duydt =  duydt - (2.0d - qshear)*Omega*ux-cs^2*1/rho0*ii*ky*rho
;;  Add gravitational acceleration.
    duxdt =  duxdt - ii*kx*phi
    duydt =  duydt - ii*ky*phi
;;  Add viscosity.
    duxdt  = duxdt  - nu*(kx^2+ky^2)*ux
    duydt  = duydt  - nu*(kx^2+ky^2)*uy
;;  Continuity equation.
    drhodt = drhodt - rho0*ii*(kx*ux+ky*uy)
;;  Add mass diffusion.
    drhodt = drhodt - nu*(kx^2+ky^2)*rho

    ds=ds+1

    ux  = ux  + duxdt *dt_beta[itsub]
    uy  = uy  + duydt *dt_beta[itsub]
    rho = rho + drhodt*dt_beta[itsub]
    t = t + ds  *dt_beta[itsub]
  endfor

endfor


end
