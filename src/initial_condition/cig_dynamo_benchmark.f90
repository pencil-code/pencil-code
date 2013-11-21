! $Id: cig_dynamo_benchmark.f90 21304 2013-11-19 22:10:55Z joern $
!
!  Set the initial conditions for the CIG Community Accuracy & 
!  Performance Benchmark: 
!  http://geodynamics.org/cig/community/workinggroups/geodyn/benchmark/
!  following Christensen et al. Physics of the Earth and Planetary 
!  Interiors, 128, 25-34 (2001).
!
!  19-nov-13/joern: coded with adaptations from spherical_convection.f90za
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!za
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
  real :: Tin=5000., rhoin=1., A=0.1
! non-dimensional paramethers
  real :: Ekman=1e-3, Rayleigh=100., Prandtl=1.0, mag_Prandtl=5.0
!
  namelist /initial_condition_pars/ &
      Ekman,Rayleigh,Prandtl, mag_Prandtl,Tin,rhoin, A
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
           "$Id: spherical_convection.f90 21304 2013-11-15 22:10:55Z pkapyla $")
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
      real, dimension (nx) :: TT, dlnTdr, lnrho, dlnrhodr,xx, ss_prof
      real :: rin, rout, DeltaT, chi, eta, nu
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
!  setting DeltaT
!
      DeltaT=Tin/400.
!
!  radial temperature profile
!
      TT=Tin-DeltaT*(rout*rin/x(l1:l2)-rin)
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
!  loop over the pencils
!
      do m=m1,m2
      do n=n1,n2
!
!  using the full profile to calculate ss and put it and lnrho
!
          TT=TT+DeltaT*(201.*A/sqrt(17920*pi)*(1.-3*xx**2+3*xx**4-xx**6) &
               *sin(y(m))**4*cos(4*z(n)))
          ss_prof=log(TT*cv*gamma*(gamma-1.))/gamma - & 
              (gamma-1.)/(gamma)*(lnrho-log(rho0))
          f(l1:l2,m,n,ilnrho) = lnrho
          f(l1:l2,m,n,iss) = ss_prof
!
!  initial magnetic field
!  Br =  5./8.*(8.*rout-6.*r-2*rin**4/r**3)*cos(y(n))
!  Btheta = -5./8.*(8.*rout-9.*r+rin**4/r**3)*sin(y(n))
!  Bphi = 5*sin(pi*(r-rin))*sin(2*y(n))
!
          f(l1:l2,m,n,iax) = 0.
          f(l1:l2,m,n,iay) = (5./pi**2/x(l1:l2)*sin(pi*(x(l1:l2)-rin)) &
                             -5./pi*cos(pi*(x(l1:l2)-rin)))*sin(2*y(n))
          f(l1:l2,m,n,iaz)=5./4.*(8.*rout-6.*x(l1:l2)-2*rin**4/x(l1:l2)**3)*sin(y(n))
      enddo
      enddo
!
!  Calculate the viscosity and chi
!
      nu=sqrt(DeltaT/Tin*gravx*rout*Ekman/Rayleigh*Lxyz(1)**3)
      chi=nu/Prandtl
      eta=nu/mag_Prandtl

!
!  print outs to put in run.in
!
       if (iproc .eq. root) then
         print*,''
         print*,'viscosity =', nu
         print*,'rotation rate =', nu/Ekman/Lxyz(1)**2
         print*,'mag. diffusivity =', eta
         print*,'thermal diffusivity =', chi
         print*,'density stratification =',exp(lnrho(1)-lnrho(nx))
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
