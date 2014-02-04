! $Id: cig_dynamo_benchmark.f90 21304 2013-11-19 22:10:55Z joern $
!
!  Set the initial conditions for the CIG Community Accuracy & 
!  Performance Benchmark: 
!  http://geodynamics.org/cig/community/workinggroups/geodyn/benchmark/
!  following Christensen et al. Physics of the Earth and Planetary 
!  Interiors, 128, 25-34 (2001).
!
!  19-nov-13/joern: coded with adaptations from spherical_convection.f90
!
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use EquationOfState
  use Sub, only: step, der_step
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: rhoin=1., A=0.1, DeltaT0=0.1, Tout
! non-dimensional paramethers
  real :: Ekman=1e-3, Rayleigh=100., Prandtl=1.0, mag_Prandtl=5.0
!
  namelist /initial_condition_pars/ &
      Ekman,Rayleigh,Prandtl, mag_Prandtl,DeltaT0, rhoin, A, Tout
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      if (lroot) call svn_id( &
           "$Id: cig_dynamo_benchmark.f90 21304 2013-11-15 22:10:55Z pkapyla $")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f)
!
!  Initializes all the f arrays in one call. This subroutine is called last.
!
!  19-nov-13/joern: coded + adapted from spherical_convection.f90
!
      use SharedVariables, only: get_shared_variable
      use EquationOfState, only: gamma, gamma_m1, rho0, cs20
      use General, only: safe_character_assign
      use Mpicomm, only: stop_it
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (nx) :: TT, dlnTdr, lnrho, dlnrhodr,xx, ss_prof, TT_prof
      real :: rin, rout, chi, eta, nu, Bnorm, DeltaT
      real, pointer :: gravx, cp, cv
      integer :: ierr, unit=1, i
!
      character (len=120) :: wfile
!
!     Retrieve cp, cv, and gravx
!
      call get_shared_variable('cp', cp, ierr)
      if (ierr/=0) call stop_it(" initialize_initial_condition: "//&
           "there was a problem when getting cp")
      call get_shared_variable('cv', cv, ierr)
      if (ierr/=0) call stop_it(" initialize_initial_condition: "//&
           "there was a problem when getting cv")
      call get_shared_variable('gravx', gravx, ierr)
      if (ierr/=0) call stop_it(" initialize_initial_condition: "//&
           "there was a problem when getting gravx")
!
!  radial extent
!
      rin=x0
      rout=x0+Lxyz(1)
      xx=2*x(l1:l2)-rin-rout
!
!  setting T0
!  T0 is related to gravx to be sufficient in an anelastic regime.
!     d*gravx*rout/cs = 0.01
!     cs2= around 100*d^2*gravx^2*rout^2. ; cs2top=cs20*(Tout)*cv*gamma*(gamma-1.)
!
!    cs2top = 1.e3*Lxyz(1)**2.*gravx**2.*rout**2.
!
!    Tout = cs2top/cs20/(cv*gamma*(gamma-1.))
    DeltaT = DeltaT0*Tout
!
!  radial temperature profile
!
      TT=Tout+DeltaT*(rout*rin/x(l1:l2)-rin)
!
!  derivative
!     
      dlnTdr=-rout*rin/x(l1:l2)**2*DeltaT/TT
!
!  calculating then density profile
!
      lnrho(1)=log(rhoin/rho0)
      dlnrhodr=-dlnTdr-gravx*x(l1:l2)/(cv*(gamma-1)*TT)
      do i=2, nx
        lnrho(i)=lnrho(i-1)+dlnrhodr(i-1)/dx_1(i-1)
      enddo
!
!  Calculate the viscosity and chi
!
      nu=sqrt(DeltaT/(Tout+DeltaT)*gravx*rout*Ekman/Rayleigh*Lxyz(1)**3)
      chi=nu/Prandtl
      eta=nu/mag_Prandtl
!
!  normalisation of the magentic field
!
      Bnorm=rho0*rhoin*eta

!  loop over the pencils
!
      do m=m1,m2
      do n=n1,n2
!
!  using the full profile to calculate ss and put it togehter with lnrho in the f-array
!
        TT_prof = TT+DeltaT*(201.*A/sqrt(17920.*pi)*(1.-3.*xx**2.+3.*xx**4.-xx**6.) &
                       *sin(y(m))**4*cos(4.*z(n)))
!
        ss_prof = log(TT_prof*cv*gamma*(gamma-1.))/gamma - & 
              (gamma-1.)/(gamma)*(lnrho-log(rho0))
!
        f(l1:l2,m,n,ilnrho) = lnrho
        f(l1:l2,m,n,iss) = ss_prof
!
!  initial magnetic field for the insulation boundary case
!  Br =  5./8.*(8.*rout-6.*r-2*rin**4/r**3)*cos(y(m))
!  Btheta = -5./8.*(8.*rout-9.*r+rin**4/r**3)*sin(y(m))
!  Bphi = 5*sin(pi*(r-rin))*sin(2*y(m))
!
!        f(l1:l2,m,n,iax) = 0.
!        f(l1:l2,m,n,iay) = (5./pi**2/x(l1:l2)*sin(pi*(x(l1:l2)-rin)) &
!                             -5./pi*cos(pi*(x(l1:l2)-rin)))*sin(2*y(m))*0.003
!        f(l1:l2,m,n,iaz)=5./4.*(8.*rout-6.*x(l1:l2)-2*rin**4/x(l1:l2)**3)*sin(y(m))*0.003
!
!  initial magnetic field for the pseudo-vacum boundary case
!  The magnetic field is directly expressed in terms of the Vectorpotential
!  following Jackson (2013) http://jupiter.ethz.ch/~ajackson/pseudo.pdf
!  taking f1=f2=K=0
!
        f(l1:l2,m,n,iax) = 1./sqrt(2.)*15./16.*x(l1:l2)*sin(pi*(x(l1:l2)-rin))*cos(2*y(m))*bnorm
        f(l1:l2,m,n,iay) = 0.
        f(l1:l2,m,n,iaz) = 1./sqrt(2.)*5./16.*(-48*rin*rout+(4*rout+rin*(4+3*rout))*6*x(l1:l2) &
                           -4*(4+3*(rin+rout))*x(l1:l2)**2.+9*x(l1:l2)**3.)*sin(y(m))*bnorm
!
      enddo
      enddo

!
!  print outs to put in run.in
!
       if (iproc .eq. root) then
         print*,''
         print*,'cs2top = ', cs20*(Tout)*cv*gamma*(gamma-1.)
         print*,'cs2bot = ', cs20*(Tout+DeltaT)*cv*gamma*(gamma-1.)
         print*,'rotation rate =', nu/Ekman/Lxyz(1)**2
         print*,'viscosity =', nu
         print*,'mag. diffusivity =', eta
         print*,'thermal diffusivity =', chi
         print*,'density stratification =',exp(lnrho(1)-lnrho(nx))
         print*,'B normalization =', Bnorm
         print*,''
       endif
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: gamma,gamma_m1,gamma1,cs20,rho0,lnrho0
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!********************************************************************
    subroutine read_initial_condition_pars(unit,iostat)
!
!  07-may-09/wlad: coded
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=initial_condition_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=initial_condition_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
