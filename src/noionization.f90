! $Id: noionization.f90,v 1.7 2003-03-29 11:14:32 brandenb Exp $

!  Dummy routine for noionization

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  !  These parameters are used if lionization were .true.
  real :: lnTT_ion,lnrho_ion,ss_ion

  !  lionization initialized to .false.
  !  cannot currently be reset to .true. in namelist
  !  because the namelist is now not even read
  logical :: lionization=.false.

  ! input parameters
  integer :: dummy 
  namelist /ionization_init_pars/ dummy

  ! run parameters
  namelist /ionization_run_pars/  dummy

  contains

!***********************************************************************
    subroutine register_ionization()
    endsubroutine register_ionization
!***********************************************************************
    subroutine initialize_ionization()
    endsubroutine initialize_ionization
!***********************************************************************
    subroutine thermodynamics(lnrho,ss,cs2,TT1,cp1tilde,Temperature)
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
      use Sub
!
      real, dimension (nx),optional :: Temperature
      real, dimension (nx) :: lnrho,ss,cs2,TT1,cp1tilde
      real, dimension (nx) :: TT,dlnPdlnrho,dlnPdS,rho,ee
      real :: ss0=-5.5542
!
      intent(out) :: Temperature
!
!  calculate cs2, 1/T, and cp1tilde
!  leave this in, in case we may activate it
!
      if(lionization) then
        if(headtt) print*,'thermodynamics: assume cp is not 1'
        dlnPdlnrho=gamma
        dlnPdS=gamma1
        TT=exp(gamma1*(lnrho+ss-ss0)); TT1=1./TT
        cs2=ss_ion*TT*dlnPdlnrho
        cp1tilde=dlnPdS/dlnPdlnrho
        if(present(Temperature)) Temperature=TT
      else
!
!  if ionization turned off, continue assuming cp=1
!  with IONIZATION=noionization this is the only option
!
        !if(headtt) print*,'thermodynamics: assume cp=1'
        cs2=cs20*exp(gamma1*(lnrho-lnrho0)+gamma*ss)
        TT1=gamma1/cs2            ! 1/(c_p T) = (gamma-1)/cs^2
        if(present(Temperature)) Temperature=cs2/gamma1
        cp1tilde=1.
      endif
!
!  calculation of internal energy
!
      if (ldiagnos.and.i_eth/=0) then
        rho=exp(lnrho)
        ee=cs2/(gamma*gamma1) !(not correct with ionization)
        call sum_mn_name(rho*ee,i_eth)
      endif
!
    endsubroutine thermodynamics
!***********************************************************************
endmodule ionization
