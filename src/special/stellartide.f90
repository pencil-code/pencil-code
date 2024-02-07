! $Id: stellartide.f90,v 1.1 2014/03/13 16:31:16 wlyra Exp $

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
!
!***************************************************************

!-------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the 
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code CVS repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or 
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module 
! selections to say something like:
!
!    SPECIAL=special/nstar
!
! Where nstar it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------

!
!  This file adds a 2D turbulent potential to emulate the effects of
!  turbulence in a simulation, in a more realistic way than simply 
!  using alpha viscosity. The description can be found in the following papers
!  
!    Laughlin, G., Steinacker, A., & Adams, F. C. 2004, ApJ, 608, 489
!    Ogihara, M., Ida S., & Morbidelli, A. 2007, Icarus, 188, 522
!    Baruteau, C., & Lin, D. 2010, ApJ, 709, 759
!    Horn, R. B., Lyra, W., Mac Low, M.-M., & Sandor, Zs. 2012, ApJ, 750, 34
!  
!  03-oct-12/wlad: coded
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
!
  implicit none
!
  include '../special.h'
!
  real :: Omegap=1.0,mstar=1000.
  logical :: lgravity_second_order=.true.
  logical :: lgravity_third_order=.true.
  logical :: lgravity_fourth_order=.true.
!
  namelist /special_init_pars/ Omegap,mstar,&
       lgravity_second_order,lgravity_third_order,lgravity_fourth_order
!
  namelist /special_run_pars/ Omegap,mstar,&
       lgravity_second_order,lgravity_third_order,lgravity_fourth_order
!
  real, dimension(nx,ny,2) :: gravity
  real, dimension(nx,3) :: fcoriolis,fgravity
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
!
!  14-jul-09/wlad: coded
!
      use Cdata
      use FArrayManager, only: farray_register_auxiliary
!
      if (lroot) call svn_id( &
           "$Id: stellartide.f90,v 1.1 2014/03/13 16:31:16 wlyra Exp $")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  14-jul-09/wlad: coded
!
      use Cdata 
      use Mpicomm
      use EquationOfState, only: cs0
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: mu,mu13,mu23,ym
      real, dimension (nx) :: rr_cyl,rr2_cyl,rr3_cyl
      integer :: mm1
!
      call keep_compiler_quiet(f)
!
      rr_cyl=x(l1:l2)
      rr2_cyl=rr_cyl**2
      rr3_cyl=rr_cyl**3
!
      mu=1./mstar
      mu13=mu**(1./3)
      mu23=mu**(2./3)
!
      do m=1,ny
        mm1=m-1+m1
        ym=y(mm1)

        gravity(:,m,1) = -1/rr2_cyl
        gravity(:,m,2) = 0.

        if (lgravity_second_order) then
          gravity(:,m,1) = gravity(:,m,1) + 1.5*rr_cyl*(1+cos(2*ym))
          gravity(:,m,2) = gravity(:,m,2) - 1.5*rr_cyl*   sin(2*ym)
        endif
!
        if (lgravity_third_order) then
          gravity(:,m,1) = gravity(:,m,1) - 3./8 * mu13 * rr2_cyl * (3*cos(ym) -5*cos(3*ym))
          gravity(:,m,2) = gravity(:,m,2) + 3./8 * mu13 * rr2_cyl * (3*sin(ym) -5*sin(3*ym))
        endif
!
        if (lgravity_fourth_order) then
          gravity(:,m,1) = gravity(:,m,1) + 1./16 * mu23 * rr3_cyl * (9 + 20*cos(2*ym) + 35*cos(4*ym))
          gravity(:,m,2) = gravity(:,m,2) - 1./16 * mu23 * rr3_cyl * (    10*sin(2*ym) + 35*sin(4*ym))
        endif
!
      enddo

    endsubroutine initialize_special
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=special_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************   
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine pencil_criteria_special()
!
      use Mpicomm
!
      lpenc_requested(i_uu)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Add tides and Coriolis force to gas.
!
      use Sub, only: grad
!
      real, dimension(mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      fgravity(:,1) = gravity(:,m-m1+1,1)
      fgravity(:,2) = gravity(:,m-m1+1,2)
!
      fcoriolis(:,1) = -2*Omegap*p%uu(:,2)
      fcoriolis(:,2) =  2*Omegap*p%uu(:,1)
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   momentum equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
      use Cdata
      use Diagnostics
!      
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
      integer :: j,ju
!
!  Modified momentum equation
!
      do j=1,3 
        ju=j+iuu-1
        df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) + fgravity(:,j) - fcoriolis(:,j)
      enddo
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
!***********************************************************************
!********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special

