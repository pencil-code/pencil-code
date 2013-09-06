! $Id: $
!
!  Isentropic initial condition (density, entropy) for convection in
!  spherical coordinates. Produces an isentropic stratification with
!  given surface gravity and surface temperature, and a heat
!  conduction profile proportional to r^-15. The setup was originally
!  introduced in Astron. Nachr., 332, 883.
!
!  02-sep-13/pete: coded
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
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: star_luminosity=1.0, Rstar=1.0
  real :: xi0=1.0, npoly1=1.5, npoly_jump=1.0, nad=1.5
  real :: Fbottom
!
  namelist /initial_condition_pars/ &
      star_luminosity, Rstar, nad, npoly1, npoly_jump, xi0
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
           "$Id: baroclinic_init.f90 19193 2012-06-30 12:55:46Z wdobler $")
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
!  02-sep-13/pete: coded
!
      use SharedVariables, only: get_shared_variable
      use EquationOfState, only: gamma, gamma_m1, rho0, cs20
      use General, only: safe_character_assign
      use Mpicomm, only: stop_it
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (mx) :: TT, dTdr, rho_prof, ss_prof, cs2_prof 
      real, dimension (nxgrid) :: kappa, gkappa, npoly2, gnpoly2
      real :: T00, rho00, Rsurf, Tsurf, coef1, L00, sigma, cs2_surf
      real :: Lsun=3.84e26, Rsun=7e8, Omsun=2.6e-6, Msun=2e30
      real :: GG=6.67348e-11, rhosun=200., fluxratio, Omsim, gratio, rratio
      real, pointer :: gravx, cp, cv
!
      integer :: i, n, m, q, ix, ierr, unit=1
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
!  Surface of the star
!
      Rsurf=x0+Lxyz(1)
!
!  Temperature using polytropic index npoly1
!
      TT(l1:l2)=gravx/(cv*(gamma-1.))*(xi0/Rstar + 1./(npoly1+1.)*(1./x(l1:l2) - 1./Rsurf))
      T00=gravx/(cv*(gamma-1.))*(xi0/Rstar + 1./(npoly1+1.)*(1./x0 - 1./Rsurf))
      Tsurf=gravx/(cv*(gamma-1.))*xi0/Rstar
!
!  Density stratification assuming an isentropic atmosphere with ss=0. 
!
      rho_prof=(TT/T00)**(1./(gamma-1.))
      rho00=(T00/T00)**(1./(gamma-1.))
!
!  Renormalize entropy with rho0 and cs20
!
      cs2_prof=TT*cv*gamma*(gamma-1.)
      ss_prof=log(cs2_prof/cs20)/gamma - (gamma-1.)/(gamma)*log(rho_prof/rho0)
!
!  Put lnrho and ss into the f-array
!
      do m=m1,m2
      do n=n1,n2
          f(l1:l2,m,n,ilnrho) = log(rho_prof(l1:l2))
          f(l1:l2,m,n,iss) = ss_prof(l1:l2)
      enddo
      enddo
!
!  Compute hcond and glhcond using global x-array
!
      coef1=star_luminosity*rho0*sqrt(gravx*Rstar)*cv*(gamma-1.)/(4.*pi)
!
      do n=1,nxgrid
         npoly2(n)=npoly_jump*(xglobal(nghost+n)/x0)**(-15.)+nad-npoly_jump
         gnpoly2(n)=15./xglobal(nghost+n)*(nad-npoly_jump-npoly2(n))
      enddo
!
      kappa=coef1*(npoly2+1.)
      gkappa=coef1*gnpoly2
!
!  Write kappa and gkappa to file to be read by run.in
!
      if (iproc .eq. root) then 
         call safe_character_assign(wfile,'hcond_glhc.dat')
         open(unit,file=wfile,status='unknown')
         do ix=1,nxgrid
            write(unit,'(2(2x,1pe12.5))') kappa(ix),gkappa(ix)
         enddo
         close(unit)
      endif
!
!  Compute flux at the bottom and a modified Stefan-Boltzmann constant
!  assuming the the total flux at the outer boundary radiates through
!  the surface for the Fgs boundary condition.
!
      L00=star_luminosity*rho0*gravx**1.5*sqrt(Rstar)
      Fbottom=L00/(4.*pi*x0**2)
      sigma=(L00/(4.*pi*Rsurf**2))/Tsurf**4
      cs2_surf=Tsurf*cv*gamma*(gamma-1.)
!
!  Compute the ratio of the dimensionless luminosity in the simulation
!  versus the Sun in order to determine a rotation rate for the
!  simulation which reproduces the rotational effect in the Sun.
!
      fluxratio=star_luminosity*rhosun*(GG*Msun)**1.5*sqrt(Rsun)/Lsun
      gratio=gravx/(GG*Msun)
      rratio=Rsun/Rstar
      Omsim=fluxratio**(1./3.)*sqrt(gratio)*rratio**(1.5)*Omsun
!
      if (iproc .eq. root) then 
         print*,'initial_condition: Fbottom   =',Fbottom
         print*,'initial_condition: SigmaSBt  =',sigma
         print*,'initial_condition: cs2top    =',cs2_surf
         print*,'initial_condition: fluxratio =',fluxratio
         print*,'initial_condition: Omsim     =',Omsim
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
