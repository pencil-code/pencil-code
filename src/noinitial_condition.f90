! $Id: nospecial.f90 10661 2009-05-02 14:38:19Z ajohan@strw.leidenuniv.nl $

!  This module provide a way for users to specify custom initial conditions
!
!  The module provides a set of standard hooks into the Pencil-Code and
!  currently allows the following customizations:
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
! CPARAM logical, parameter :: linitial_condition = .false.
!
!***************************************************************
!
!--------------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above 
!    linitialcondition = .true.
! to enable use of special hooks.
!
! The rest of this file may be used as a template for your own
! initial conditions. Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. mhs_equilibrium.f90 
! and place it in the $PENCIL_HOME/src/initial_condition directory. 
! This path has been created to allow users ot optionally check their 
! contributions in to the Pencil-Code CVS repository.  This may be useful 
! if you are working on/using an initial condition with somebody else or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional initial condition code edit the Makefile.local
! in the src directory under the run directory in which you wish to
! use your initial condition.  Add a line with all the module
! selections to say something like:
!
!    INITIAL_CONDITION=initial_condition/mhs_equilibrium
!
! Where mhs_equilibrium is replaced by the filename of your new file
! upto and not including the .f90
!
! This module is based on Tony's special module.
!
!--------------------------------------------------------------------

module InitialCondition

  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'initial_condition.h'

!!  integer :: dummy

!!! input parameters
!!  namelist /initial_condition_init_pars/ dummy

  contains

!***********************************************************************
    subroutine register_initial_condition()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/wlad: coded
!
      use Cdata
!
!  Identify CVS/SVN version information.
!
      if (lroot) call cvs_id( &
           "$Id: nospecial.f90 10661 2009-05-02 14:38:19Z ajohan@strw.leidenuniv.nl $")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub,     only: keep_compiler_quiet
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Initialize any module variables which are parameter dependent
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

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
!  Initialize logarithmic density. init_lnrho 
!  will take care of converting it to linear 
!  density if you use ldensity_nolog
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
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
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic potential.
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine initial_condition_aatest(f)
!
!  Initialize testfield
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aatest
!********************************************************************
    subroutine initial_condition_uutest(f)
!
!  Initialize testflow
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uutest
!********************************************************************
    subroutine initial_condition_lncc(f)
!
!  Initialize passive scalar
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lncc
!********************************************************************
    subroutine initial_condition_chiral(f)
!
!  Initialize chiral
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chiral
!********************************************************************
    subroutine initial_condition_chemistry(f)
!
!  Initialize chemistry
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chemistry
!********************************************************************
    subroutine initial_condition_uud(f)
!
!  Initialize dust fluid velocity
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uud
!********************************************************************
    subroutine initial_condition_nd(f)
!
!  Initialize dust fluid density
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_nd
!********************************************************************
    subroutine initial_condition_uun(f)
!
!  Initialize neutral fluid velocity
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uun
!********************************************************************
    subroutine initial_condition_lnrhon(f)
!
!  Initialize neutral fluid density
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrhon
!********************************************************************
    subroutine initial_condition_ecr(f)
!
!  Initialize cosmic rays
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ecr
!********************************************************************
    subroutine initial_condition_fcr(f)
!
!  Initialize cosmic ray flux
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_fcr
!********************************************************************
    subroutine initial_condition_solid_cells(f)
!
!  Initialize solid cells
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_solid_cells
!********************************************************************
    subroutine initial_condition_cctest(f)
!
!  Initialize testscalar
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_cctest
!********************************************************************
    subroutine initial_condition_xxp(f,fp)
!
!  Initialize particles' positions
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_xxp
!***********************************************************************
    subroutine initial_condition_vvp(f,fp)
!
!  Initialize particles' velocities
!
!  07-may-09/wlad: coded
!
      use Cdata
      use Sub, only: keep_compiler_quiet
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_vvp
!***********************************************************************
    subroutine read_initial_condition_init_pars(unit,iostat)
!
      use Sub, only: keep_compiler_quiet
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) call keep_compiler_quiet(iostat)
      call keep_compiler_quiet(unit)

    endsubroutine read_initial_condition_init_pars
!***********************************************************************
    subroutine write_initial_condition_init_pars(unit)
!
      use Sub, only: keep_compiler_quiet
!
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_initial_condition_init_pars
!***********************************************************************

!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include 'initial_condition_dummies.inc'
!********************************************************************

endmodule InitialCondition

