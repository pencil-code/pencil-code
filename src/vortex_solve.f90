!
! Functions and subroutines to use with KEPVOR subroutine in INITCOND.
! Coded 26-03-2003 (AJ)
!
module PARAMS
!
! Parameters to be saved. All outside calls to functions and subroutines
! must transfer necessary parameters.
! Functions called by NR library programs can not receive parameters
! in other ways
!
  implicit none
  real :: x,y,z,a_ell,b_ell,c_ell,q_ell,kappa,qshear
  save :: x,y,z,a_ell,b_ell,c_ell,q_ell,kappa,qshear
  
end module PARAMS
!***********************************************************************
real function KEPVORFCN (x,y,c_ell)
!
! Calculates xi(x,y,c) as solution to 
! x**2/(c**2*cosh(xi)**2)+y**2/(c**2*sinh(xi)**2)=1.
!
  use PARAMS
  implicit none
  real :: x_du,y_du,c_ell_du
  real :: xi, xi1, xi2, fn, df, RTSAFE
  logical :: succes
  real, parameter :: tol=1.e-5
  external ZBRAC, RTSAFE

  x=x_du
  y=y_du
  c_ell=c_ell_du
!
! Initial guess for bracketing root (need not be correct)
!
  xi1=1.
  xi2=10.
!
! ZBRAC brackets the root
!
  call ZBRAC(XIFUNC,xi1,xi2,succes)
!
! RTSAFE solves for xi
!
  KEPVORFCN = RTSAFE(XIFUNC,xi1,xi2,tol)

end function KEPVORFCN
!***********************************************************************
subroutine XIFUNC (xi,fn,df)
!
! Calculates x**2/(c**2*cosh(xi)**2)+y**2/(c**2*sinh(xi)**2)-1
! and derivative for use in RTSAFE
!
  use PARAMS
  real :: xi, fn, df

  fn = x**2/(c**2*cosh(xi)**2)+y**2/(c**2*sinh(xi)**2)-1.
  df = x**2/c**2*(-2/cosh(xi)**3)*sinh(xi)+y**2/c**2*(-2/sinh(xi)**3)*cosh(xi)

end subroutine XIFUNC
!***********************************************************************
real function PSIDER (dertype,x_du,y_du,z_du,xi0_du,a_ell_du,b_ell_du,c_ell_du,q_ell_du,kappa_du)
!
! Calculates derivative of Psi. dertype=PSIX gives dPsi/dx,
! dertype=PSIX gives dPsi/dy
!
  implicit none
  real :: x_du,y_du,z_du,a_ell_du,b_ell_du,c_ell_du,q_ell_du,kappa_du
  external DFRIDR, dertype

  x=x_du
  y=y_du
  z=z_du
  xi0=xi0_du
  a_ell=a_ell_du
  b_ell=b_ell_du
  c_ell=c_ell_du
  q_ell=q_ell_du
  kappa=kappa_du
!
! DRIDR calculates derivative of known function
!
  call DFRIDR(dertype,x,h,err)

end function PSIDER
!***********************************************************************
function PSI (x_du,y_du,z_du,xi0_du,a_ell_du,b_ell_du,c_ell_du,q_ell_du,kappa_du)!
! Calculate Psi(x,y,z)
!
  implicit none
  use PARAMS
  real :: x_du,y_du,z_du,a_ell_du,b_ell_du,c_ell_du,q_ell_du,kappa_du
  real :: xi, eta

  x=x_du
  y=y_du
  z=z_du
  xi0=xi0_du
  a_ell=a_ell_du
  b_ell=b_ell_du
  c_ell=c_ell_du
  q_ell=q_ell_du
  kappa=kappa_du
!
! Find xi(x,y,z)
!
  xi = KEPVORFCN ()
  if (xi .lt. 0) then xi=-xi
!
! eta is solution to x=c cosh(xi) cos(eta)
!
  eta = acos(x/(c*cosh(xi)))
!
! acos gives result from 0 to pi. If y is negative, result must be
! subtracted from 2 pi to give proper eta.
!
  if (y .lt. 0) then eta = 2*pi-eta
!
! Stream function inside and outside vortex
!
  if (xi .le. xi0) then
    PSI = kappa*b_ell**2/(2*q_ell)*(q_ell+1.)*(cosh(xi)**2*cos(eta)**2 + &
          q_ell**2*sinh(xi)**2*sin(eta)**2
  else
    PSI = 0.25*kappa*b_ell**2*(q_ell**2-1.)*sinh(xi)**2*(1.-cos(2*eta)) + &
          0.25*b**22*kappa*(q_ell+1.)/(q_ell-1.) + &
          0.5*b_ell**2*kappa*(q_ell+1.)/(q_ell-1.)*(xi-xi0) + &
          0.25*kappa*b_ell**2*exp(2*(xi0-xi))*cos(2*eta)
  endif

end function PSI
!***********************************************************************
subroutine PSIX (x_dfridr)
!
! Used for dPsi/dx
!
  implicit none
  real :: x_dfridr

  x=x_dridr

  PSIX = PSI(x,y,z,xi0,a_ell,b_ell,c_ell,q_ell,kappa)

end subroutine PSIX
!***********************************************************************
subroutine PSIY (y_dfridr)
!
! Used for dPsi/dy
!
  implicit none
  real :: y_dfridr

  y=y_dridr

  PSIY = PSI(x,y,z,xi0,a_ell,b_ell,c_ell,q_ell,kappa)

end subroutine PSIY
