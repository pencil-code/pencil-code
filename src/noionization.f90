! $Id: noionization.f90,v 1.14 2003-05-05 18:48:52 brandenb Exp $

!  Dummy routine for noionization

module Ionization

  use Cparam
  use Cdata
  use Density

  implicit none

  ! global ionization parameter for yH (here a scalar set to 0)
  real :: yyH=0.

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
    subroutine ionfrac(f)
      real, dimension (mx,my,mz,mvar), intent(in) :: f
      if(ip==0) print*,f(1,1,1,1)  !(keep compiler quiet)
    endsubroutine ionfrac
!***********************************************************************
    subroutine output_ionization(lun)
      integer, intent(in) :: lun
      if(ip==0) print*,lun  !(keep compiler quiet)
    endsubroutine output_ionization
!***********************************************************************
    subroutine thermodynamics(lnrho,ss,cs2,TT1,cp1tilde, &
      Temperature,InternalEnergy,IonizationFrac)
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
      real, dimension (nx), intent(out), optional :: Temperature
      real, dimension (nx), intent(out), optional :: InternalEnergy
      real, dimension (nx), intent(out), optional :: IonizationFrac
      real, dimension (nx) :: lnrho,ss,cs2,TT1,cp1tilde
      real, dimension (nx) :: TT,dlnPdlnrho,dlnPdS,rho,ee
      real :: ss0=-5.5542
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
        if(ldiagnos) ee=cs2/(gamma*gamma1)
      endif
!
!  optional output
!
      if(present(Temperature)) Temperature=1./TT1
      if(present(InternalEnergy)) InternalEnergy=ee
      if(present(IonizationFrac)) IonizationFrac=0.
!
    endsubroutine thermodynamics
!***********************************************************************
    subroutine ioncalc(lnrho,ss,yH,dlnPdlnrho,dlnPdss,TT,kappa)
!
!   calculates thermodynamic quantities under partial ionization
!
!   28-mar-03/tobi: added kappa
!
      real, dimension(nx),intent(in)   :: lnrho,ss,yH
      real, dimension(nx), optional    :: dlnPdlnrho,dlnPdss,TT,kappa
                           intent(out) :: dlnPdlnrho,dlnPdss,TT,kappa
      real                             :: ss0=-5.5542
!
!  if kappa is needed, use electron scattering value, kappa_es
!
      if (present(dlnPdlnrho)) dlnPdlnrho=gamma
      if (present(dlnPdss)) dlnPdss=gamma1
      if (present(TT)) TT=exp(gamma1*(lnrho+ss-ss0))
      if (present(kappa)) kappa=kappa_es
      if (ip==0) print*,yH  !(to keep compiler quiet)
    endsubroutine ioncalc
!***********************************************************************
endmodule ionization
