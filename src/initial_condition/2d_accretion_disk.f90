! $Id$
!
!  21-feb-23/hongzhe: adapted from initial_condition/centrifugal_balance.f90
! 
!  This module sets up the initial condition for an axisymmetric accretion
!  disk in spherical coordinates, using the following structure:
!  (denoting the radial distance in spherical coordinates in r)
!  (and that in cylindrical coordinates in R)
!
!  gravity:  potential propto r**(-1)
!  
!  In the R direction:
!  velocity: keplerian R**(-0.5), plus noise propto cs, i.e. also R**(-0.5)
!  sound speed: cs propto R**(-0.5)
!  density: propto R**(-1.5)
!  pressure: propto R**(-2.5)
!
!  In the z direction:
!  All the quantities depend on a single variable l=z/R=1/tan(theta) and
!  have some analytical form.
!  
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
!    INITIAL_CONDITION =   initial_condition/mhs_equilibrium
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
! global variables to this module
!
  real :: R0=1.0
  real :: gamma
  real :: ct0   !  thermal sound speed. ct0=cs0/sqrt(gamma)
  real :: vk0   !  Keplerian speed at R0
  real :: h0    !  scale height at R0, h0=ct0*R0/vk0
  real, dimension (mx,my,mz) :: r_cyl,r_sph,sin_th
!
  character (len=labellen) :: inituu='0',initlnrho='uniform',initaa='0'
  real :: ampluu=0.,ampluu_noise=0.
  real :: rho_exp=0., rho_rmax=0.,rho_rmin=0., amplrho=1.0
!
  namelist /initial_condition_pars/ &
       inituu,ampluu,ampluu_noise,amplrho,&
       initlnrho,rho_exp,rho_rmax,rho_rmin

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
!  21-feb-23/hongzhe: coded
!
      use EquationOfState, only: get_gamma_etc,cs20
      use Gravity, only : g0
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call get_gamma_etc(gamma)
!
      ct0 = sqrt(cs20/gamma)
      vk0 = sqrt(g0/R0)
      h0  = ct0*R0/vk0
      r_cyl = max(R0,spread(spread( x,2,my),3,mz) * spread(spread( sin(y), 1,mx),3,mz))
      r_sph = spread(spread(x,2,my),3,mz)
      sin_th = spread(spread( sin(y), 1,mx),3,mz)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  21-feb-23/hongzhe: adapted
!
      use Gravity, only: g0
      use Initcond, only: gaunoise
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my) :: fact_y
!
      integer :: mm
      real :: lat_crit
!
      lat_crit=acos(1.0/(1.0+0.9*(xyz1(1)*gamma/(gamma-1.)/g0*ct0**2)**(1+rho_exp*(gamma-1.0))))
!
!  add noise and rescale to the local keplerian velocity
!
      if (ampluu_noise>0.0) then
        call gaunoise(ampluu_noise,f,iux,iuz)
        do mm=iux,iuz
          f(:,:,:,mm) = f(:,:,:,mm)*(r_cyl/R0)**(-0.5)
        enddo
      endif
!
      select case (inituu)
      case ('const')
        f(:,:,:,iux:iuz) = f(:,:,:,iux:iuz) + ampluu
      case ('force-balance')
        !  horizontally balancing gravity and centrifugal forces
        f(:,:,:,iuz) = f(:,:,:,iuz) + ampluu*sqrt(g0/r_sph)*sin_th
      case ('force-balance-disk')
        !  horizontally balancing gravity and centrifugal forces in the disk
        where (abs(y-pi/2.0)>=lat_crit)
          fact_y = 0.0
        elsewhere
          fact_y = sin(y)
        endwhere
        f(:,:,:,iuz) = f(:,:,:,iuz) + ampluu*sqrt(g0/r_sph)*spread(spread(fact_y,1,mx),3,mz)
      case default
        call fatal_error('initial_condition_uu', &
            'initial_condition_uu: No profile of inituu="'//trim(inituu)//'"')
      endselect
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  21-feb-23/hongzhe: adapted
!
      use FArrayManager
      use EquationOfState, only: rho0
      use Gravity, only: g0
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: epsilon
      integer :: ll
!
      select case (initlnrho)
      case ('uniform')
        f(:,:,:,ilnrho) = log(amplrho*rho0)
      case ('rpower*exp(-z2/H2)')
        epsilon = ct0/sqrt(g0/1.0)
        f(:,:,:,ilnrho) = log(amplrho*rho0) &
                              + spread(spread(rho_exp*log(x/R0),                      2,my),3,mz) &
                              + spread(spread(rho_exp*log(0.001+sin(y)),              1,mx),3,mz) &
                              - spread(spread(0.5*cos(y)**2/epsilon**2, 1,mx),3,mz) &
                              * spread(spread((x/R0)**(-rho_exp*(gamma-1.)-1),        2,my),3,mz)
        if (rho_rmin>0.0) then
          do ll=1,mx
            if (x(ll)<=rho_rmin) f(ll,:,:,ilnrho) = f(ll,:,:,ilnrho) - min(10.,(x(ll)/rho_rmin-1)**2)
          enddo
        endif
        if (rho_rmax>0.0) then
          do ll=1,mx
            if (x(ll)>=rho_rmax) f(ll,:,:,ilnrho) = f(ll,:,:,ilnrho) - min(10.,(x(ll)/rho_rmax-1)**2)
          enddo
        endif
      case default
        call fatal_error('initial_condition_lnrho', &
            'initial_condition_lnrho: No profile of initlnrho="'//trim(initlnrho)//'"')
      endselect
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize vector potential.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aa
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
