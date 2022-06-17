! $Id$
!
!  This module provide a way for users to specify custom initial
!  conditions.
!
!  The module provides a set of standard hooks into the Pencil Code
!  and currently allows the following customizations:
!
!   Description                               | Relevant function call
!  ------------------------------------------------------------------------
!   Initial condition registration            | register_special
!     (pre parameter read)                    |
!   Initial condition initialization          | initialize_special
!     (post parameter read)                   |
!                                             |
!   Initial condition for momentum            | initial_condition_uu
!   Initial condition for density             | initial_condition_lnrho
!   Initial condition for entropy             | initial_condition_ss
!   Initial condition for magnetic potential  | initial_condition_aa
!
!   And a similar subroutine for each module with an "init_XXX" call.
!   The subroutines are organized IN THE SAME ORDER THAT THEY ARE CALLED.
!   First uu, then lnrho, then ss, then aa, and so on.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!    linitial_condition = .true.
! to enable use of custom initial conditions.
!
! The rest of this file may be used as a template for your own initial
! conditions. Simply fill out the prototypes for the features you want
! to use.
!
! Save the file with a meaningful name, e.g. mhs_equilibrium.f90, and
! place it in the $PENCIL_HOME/src/initial_condition directory. This
! path has been created to allow users to optionally check their
! contributions in to the Pencil Code SVN repository. This may be
! useful if you are working on/using an initial condition with
! somebody else or may require some assistance from one from the main
! Pencil Code team. HOWEVER, less general initial conditions should
! not go here (see below).
!
! You can also place initial condition files directly in the run
! directory. Simply create the folder 'initial_condition' at the same
! level as the *.in files and place an initial condition file there.
! With pc_setupsrc this file is linked automatically into the local
! src directory. This is the preferred method for initial conditions
! that are not very general.
!
! To use your additional initial condition code, edit the
! Makefile.local in the src directory under the run directory in which
! you wish to use your initial condition. Add a line that says e.g.
!
!    INITIAL_CONDITION =   initial_condition/binarycomparison
!
! Here mhs_equilibrium is replaced by the filename of your new file,
! not including the .f90 extension.
!
! This module is based on Tony's special module.
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
  real :: g0=1.,v0=1d-4, npot=4.0, Omegab=1.0
  real :: Rcav = 2.5, density_floor = 1d-5, Rout=10.0
  real :: temperature_power_law=1.0
!
  namelist /initial_condition_pars/ g0, Rcav, density_floor, &
       Rout, npot, Omegab, v0, temperature_power_law
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
      use SharedVariables, only: put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  This subroutine is a general routine that takes
!  the gravity acceleration and adds the centrifugal force
!  that numerically balances it.
!
!  Pressure corrections to ensure centrifugal equilibrium are
!  added in the respective modules
!
!  24-feb-05/wlad: coded
!  04-jul-07/wlad: generalized for any shear
!  08-sep-07/wlad: moved here from initcond
!
      use Sub,     only: get_radial_distance
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rr_cyl,rr_sph,OO,Omega0,vphi,vr
      real, dimension (nx) :: cosp,sinp
!
      if (lroot) &
           print*,'binarycomparison: initializing velocity field'
!
      do m=m1,m2
        do n=n1,n2
!
          call  get_radial_distance(rr_cyl,rr_sph)

          cosp = x(l1:l2)/rr_cyl
          sinp = y(  m  )/rr_cyl
!        
          vr = v0*sinp*rr_cyl * exp(-(rr_cyl/3.5)**6)
!        
          Omega0 = sqrt(g0/rr_cyl**3 * (1-cs0**2))
        
          OO = (Omega0**(-npot) + Omegab**(-npot))**(-1.0/npot)
!        
          vphi = rr_cyl*OO
!
          f(l1:l2,m,n,iux) = f(l1:l2,m,n,iux) + vr*cosp - vphi*sinp
          f(l1:l2,m,n,iuy) = f(l1:l2,m,n,iuy) + vr*sinp + vphi*cosp
          f(l1:l2,m,n,iuz) = f(l1:l2,m,n,iuz) + 0.
!
        enddo
      enddo
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
      use FArrayManager
      use Gravity,         only: potential,acceleration
      use Sub,             only: get_radial_distance,grad,power_law
      use EquationOfState, only: rho0
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx)   :: strat,tmp1,tmp2,cs2
      real, dimension (mx)   :: rr_sph,rr,rr_cyl,lgSigma,fr,cavity
      integer, pointer       :: iglobal_cs2,iglobal_glnTT
      integer                :: i
      logical                :: lheader
!
      if (lroot) print*,&
           'initial_condition_lnrho: locally isothermal approximation'
!
!  Pencilize the density allocation.
!
      do n=1,mz
        do m=1,my
!
          lheader=lroot.and.(m==1).and.(n==1)
!
!  Midplane density
!
          call get_radial_distance(rr_sph,rr_cyl)
          if (lcylindrical_gravity.or.lcylinder_in_a_box.or.lcylindrical_coords) then
            rr=rr_cyl
          elseif (lsphere_in_a_box.or.lspherical_coords) then
            rr=rr_sph
          endif
!
          fr = 1 - 1/(1+exp(-2*(rr_cyl-Rout)))
          cavity = (Rcav/rr_cyl)**12
          lgSigma = log(rho0) - cavity + log(fr)
!        
          f(:,m,n,ilnrho) = f(:,m,n,ilnrho) + max(lgSigma,log(density_floor))
        enddo
      enddo
!
!  Set the thermodynamical variable
!
      call set_thermodynamical_quantities&
           (f,temperature_power_law,iglobal_cs2,iglobal_glnTT)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine set_thermodynamical_quantities&
         (f,temperature_power_law,iglobal_cs2,iglobal_glnTT)
!
!  Subroutine that sets the thermodynamical quantities
!   - static sound speed, temperature or entropy -
!  based on a sound speed which is given as input.
!  This routine is not general. For llocal_iso (locally
!  isothermal approximation, the temperature gradient is
!  stored as a static array, as the (analytical) derivative
!  of an assumed power-law profile for the sound speed
!  (hence the parameter temperature_power_law)
!
!  05-jul-07/wlad: coded
!  16-dec-08/wlad: moved pressure gradient correction to
!                  the density module (correct_pressure_gradient)
!                  Now this subroutine really only sets the thermo
!                  variables.
!
      use FArrayManager
      use EquationOfState, only: gamma,gamma_m1,get_cp1,&
           cs20,cs2bot,cs2top,lnrho0
      use Sub,             only: power_law,get_radial_distance
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rr,rr_sph,rr_cyl,cs2,lnrho
      real, dimension (nx) :: gslnTT
      real :: temperature_power_law
      integer, pointer, optional :: iglobal_cs2,iglobal_glnTT
      logical :: lheader
!
      intent(in)    :: temperature_power_law
      intent(inout) :: f
!
!  Set the sound speed
!
      if (lroot) print*,'Temperature gradient with power law=',temperature_power_law
!
      nullify(iglobal_cs2)
      call farray_use_global('cs2',iglobal_cs2)
      nullify(iglobal_glnTT)
      call farray_use_global('glnTT',iglobal_glnTT)
  
      do m=1,my; do n=1,mz
        lheader=((m==1).and.(n==1).and.lroot)
        call get_radial_distance(rr_sph,rr_cyl)
        
        if (lcylindrical_gravity.or.lcylinder_in_a_box.or.lcylindrical_coords) then
          rr=rr_cyl 
        elseif (lsphere_in_a_box.or.lspherical_coords) then
          rr=rr_sph
        else
          call fatal_error("initial_condition_lnrho",&
               "no valid coordinate system")
        endif
!
        call power_law(cs20,rr,temperature_power_law,cs2,r_ref)
        f(l1:l2,m,n,iglobal_cs2) = cs2
!
        gslnTT=-temperature_power_law/((rr/r_ref)**2+rsmooth**2)*rr/r_ref**2
!
        f(l1:l2,m,n,iglobal_glnTT  )=gslnTT*x(l1:l2)/rr_cyl
        f(l1:l2,m,n,iglobal_glnTT+1)=gslnTT*y(m)    /rr_cyl
        f(l1:l2,m,n,iglobal_glnTT+2)=0.
        
      enddo
      enddo
!
!  Word of warning...
!
      if (lroot.and.llocal_iso) then
        if (associated(iglobal_cs2)) then
          print*,"Max global cs2 = ",&
               maxval(f(l1:l2,m1:m2,n1:n2,iglobal_cs2))
          print*,"Sum global cs2 = ",&
               sum(f(l1:l2,m1:m2,n1:n2,iglobal_cs2))
        endif
        if (associated(iglobal_glnTT)) then
          print*,"Max global glnTT(1) = ",&
               maxval(f(l1:l2,m1:m2,n1:n2,iglobal_glnTT))
          print*,"Sum global glnTT(1) = ",&
               sum(f(l1:l2,m1:m2,n1:n2,iglobal_glnTT))
        endif
      endif
!
      cs2bot=cs20
      cs2top=cs20
!
      if (lroot) &
           print*,"thermodynamical quantities successfully set"
!
    endsubroutine set_thermodynamical_quantities
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
