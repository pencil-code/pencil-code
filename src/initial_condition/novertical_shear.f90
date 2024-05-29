!  $Id$
!
!  Initial condition for run without vertical shear using the 
!  (unphysical) modified gravity profile of Klahr et al. (2023)
!
!  22-may-24/wlad: adapted from mhs_equilibrium.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
module InitialCondition
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: g0=1
  real :: density_power_law=1.5,temperature_power_law=1.0
!
  namelist /initial_condition_pars/ &
      g0,density_power_law,temperature_power_law
!
  real :: gamma, gamma_m1, cp1
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  07-oct-09/wlad: coded
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
   subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: get_gamma_etc

      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: cp
      call keep_compiler_quiet(f)
!
      call get_gamma_etc(gamma,cp)
      gamma_m1=gamma-1.
      cp1=1./cp

    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      use Gravity,        only: acceleration
      use Sub,            only: get_radial_distance, power_law
      use FArrayManager,  only: farray_use_global

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl,g_r,g_rad,g_zed
      real, dimension (mx) :: OOK2,OO2,tmp,H2
      real :: p,q,OOcorot
      integer, pointer :: iglobal_cs2
      integer :: ics2
!
      if (.not.lspherical_coords) call fatal_error("initial_condition_uu",&
           "This method is only for spherical coordinates. Use centrifugal_balance for cylindrical.")
!
!  Set the sound speed
!
      call set_sound_speed(f)
!
!  Get the sound speed
!
      if (llocal_iso) then 
        nullify(iglobal_cs2)
        call farray_use_global('cs2',iglobal_cs2)
        ics2=iglobal_cs2 
      elseif (lentropy) then
        ics2=iss
      elseif (ltemperature) then 
        ics2=iTT
      endif
!
!  Analytical expression that leads to an equilibrium configuration. 
!
      if (lcorotational_frame) then 
        OOcorot=rcorot**(-1.5)
      else
        OOcorot=0.
      endif
!
      p=-density_power_law
      q=-temperature_power_law
!
      do m=1,my;do n=1,mz
        call get_radial_distance(rr_sph,rr_cyl)
!
        OOK2=g0/rr_cyl**3
        H2=f(:,m,n,ics2)/(gamma*OOK2)
!
!  pressure correction
!
        OO2=OOK2*(1 + H2/rr_cyl**2*(p+q))
!
        f(:,m,n,iuz) = f(:,m,n,iuz) + rr_cyl*sqrt(OO2)
!
      enddo;enddo
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: rho0
      use FArrayManager,   only: farray_use_global
      use Gravity,         only: potential
      use Sub,             only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: rr_sph,rr_cyl
      real, dimension (mx) :: lnrhomid,strat
      real, dimension (mx) :: cs2,tmp1,tmp2
      real    :: p,q
      integer :: ics2
      integer, pointer :: iglobal_cs2
!
      if (ldensity_linearstart) call fatal_error("initial_condition_lnrho",&
           "Switch off ldensity_linearstart. This routine assumes lnrho in initial condition")
!
      p=-density_power_law 
      q=-temperature_power_law
!
      if (llocal_iso) then 
        if (lroot) print*,&
             'initial_condition_lnrho: locally isothermal approximation'
      else
        if (lroot) print*,&
             'initial_condition_lnrho: gamma=',gamma
      endif
      if (lroot) print*,'Radial density stratification with power law=',p
      if (lroot) print*,'Radial temperature stratification with power law=',q
!
!  Get the sound speed globals
!
      if (llocal_iso) then 
        nullify(iglobal_cs2)
        call farray_use_global('cs2',iglobal_cs2)
        ics2=iglobal_cs2 
      elseif (lentropy) then 
        ics2=iss
      elseif (ltemperature) then 
        ics2=iTT
      endif
!
!  Pencilize the density allocation.
!
      do n=1,mz
        do m=1,my
!
!  Midplane density
!
          call get_radial_distance(rr_sph,rr_cyl)
          lnrhomid=log(rho0)+p*log(rr_cyl/r_ref) 
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho)+lnrhomid
!
!  Vertical stratification
!
          cs2=f(:,m,n,ics2)
          strat = -.5*gamma*cotth(m)**2/(cs2*rr_cyl)
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho)+strat
!
        enddo
      enddo
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: cs20,lnrho0
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: cs2,lnrho

      do m=1,my; do n=1,mz
        if (ldensity_nolog) then 
          lnrho=log(f(:,m,n,irho))
        else
          lnrho=f(:,m,n,ilnrho)
        endif     
!
!  The sound speed is stored in the energy slot
!
        if (lentropy) then 
          cs2=f(:,m,n,iss)
          if (pretend_lnTT) then
            f(:,m,n,iss)=log(cs2*cp1/gamma_m1)
          else
            f(:,m,n,iss)=1./(gamma*cp1)*(log(cs2/cs20)-gamma_m1*(lnrho-lnrho0))
          endif
        elseif (ltemperature) then 
          cs2=f(:,m,n,iTT)
          f(:,m,n,iTT)=cs2*cp1/gamma_m1
        endif
!
      enddo;enddo
!
    endsubroutine initial_condition_ss
!***********************************************************************     
    subroutine set_sound_speed(f)
!
!  Set the thermo-related quantities. Illustrates that 
!  the user can define as many internal routines as wanted.
!
!  10-may-09/wlad : moved from initial_condition_lnrho
!
      use FArrayManager,   only: farray_use_global
      use EquationOfState, only: cs20
      use Sub,             only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray) :: f      
      real, dimension (mx) :: rr_sph,rr_cyl,cs2
      integer, pointer :: iglobal_cs2,iglobal_glnTT
      real    :: q
      integer :: ics2,iglnTT
!
      q=-temperature_power_law
!
!  Get the globals needed to store sound speed and temperature gradient
!
      if (llocal_iso) then 
        nullify(iglobal_cs2)
        call farray_use_global('cs2',iglobal_cs2);ics2=iglobal_cs2
        nullify(iglobal_glnTT)
        call farray_use_global('glnTT',iglobal_glnTT);iglnTT=iglobal_glnTT
      elseif (lentropy) then 
          ics2=iss
      elseif (ltemperature) then 
          ics2=iTT
      endif
!
!  Set the sound speed - a power law in cylindrical radius.
!
      do m=1,my
        do n=1,mz
          call get_radial_distance(rr_sph,rr_cyl)
          cs2=cs20*(rr_cyl/r_ref)**q
!
!  Store cs2 in one of the free slots of the f-array
!
          f(:,m,n,ics2)=cs2
!
!  Same for the temperature gradient
!
          if (llocal_iso) then             
            f(:,m,n,iglnTT  )=q/rr_sph
            f(:,m,n,iglnTT+1)=q/rr_sph*cotth(m)
            f(:,m,n,iglnTT+2)=0.
          endif
!
        enddo
      enddo
!
    endsubroutine set_sound_speed
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
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
