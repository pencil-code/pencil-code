! $Id: ionization.f90,v 1.1 2003-02-02 15:12:52 brandenb Exp $

!  This modules contains the routines for simulation with
!  simple hydrogen ionization.

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  real, parameter :: eV=1.602177e-12, chiH=13.6*eV
  real, parameter :: k_B=1.380658e-16, m_p=1.672623e-24, kB_over_mp=k_B/m_p
  real, parameter :: nqe = 2.414703e15, nqp = 1.899883e20

  !  lionization initialized to .true.
  !  it can be reset to .false. in namelist
  logical :: lionization=.true.

  ! input parameters
  integer :: dummy 
  namelist /ionization_init_pars/ lionization

  ! run parameters
  namelist /ionization_run_pars/ lionization

  contains

!***********************************************************************
    subroutine register_ionization()
!
!   2-feb-03/axel: adapted from Interstellar module
!
      use Cdata
      use Mpicomm, only: stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_ionization called twice')
      first = .false.
!
      if ((ip<=8) .and. lroot) then
        print*, 'Register_ionization'
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: ionization.f90,v 1.1 2003-02-02 15:12:52 brandenb Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call stop_it('Register_ionization: nvar > mvar')
      endif
!
    endsubroutine register_ionization
!***********************************************************************
    subroutine initialize_ionization()
!
!  Perform any post-parameter-read initialization, e.g. set derived 
!  parameters.
!
!   2-feb-03/axel: adapted from Interstellar module
!
      use Cdata
      use General
!
      logical, save :: first=.true.
      logical :: exist
!
    endsubroutine initialize_ionization
!***********************************************************************
    subroutine thermodynamics(lnrho,ss,cs2,TT1,cp1tilde)
!
!  Calculate thermodynamical quantities, cs2, 1/T, and cp1tilde
!  cs2=(dp/drho)_s is the adiabatic sound speed
!  TT1=1/T is the inverse temperature
!  neutral gas: cp1tilde=kB_over_mp/cp=0.4 ("nabla_ad" maybe better name)
!  in general: cp1tilde=dlnPdS/dlnPdlnrho
!
!   2-feb-03/axel: simple example coded
!
      use Cdata
      use General
!
      real, dimension (nx) :: lnrho,ss,rho1,cs2,TT1,cp1tilde
      real, dimension (nx) :: TT,dlnPdlnrho,dlnPdS
      real :: ss0=-5.5542
!
!  calculate cs2, 1/T, and cp1tilde
!
      if(lionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1'
        dlnPdlnrho=gamma
        dlnPdS=gamma1
        TT=exp(gamma1*(lnrho+ss-ss0))
        TT1=1./TT
        cs2=kB_over_mp*TT*dlnPdlnrho
        cp1tilde=dlnPdS/dlnPdlnrho
      else
!
!  if ionization turned off, continue assuming cp=1
!
        if(headtt) print*,'thermodynamics: assume cp=1'
        cs2=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)
        TT1=gamma1/cs2            ! 1/(c_p T) = (gamma-1)/cs^2
        cp1tilde=1.
      endif
!
    endsubroutine thermodynamics
!***********************************************************************
endmodule ionization
