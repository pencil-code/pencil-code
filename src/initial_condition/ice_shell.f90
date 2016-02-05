! $Id: centrifugal_balance.f90 19193 2012-06-30 12:55:46Z wdobler $
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
!!  integer :: dummy
 
  real :: Tupp=100.  ! Surface temperature in K
  real :: Tbot=270.  ! Ocean temperature in K 
  real :: ampltt=1.  ! Random initial temperature amplitude
  real :: kx_TT=2*pi ! 
  real :: kz_TT=pi   !
  character (len=labellen) :: initTT='nothing'
!
  namelist /initial_condition_pars/ Tupp,Tbot, ampltt, kx_TT, kz_TT, initTT
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
         "$Id: centrifugal_balance.f90 19193 2012-06-30 12:55:46Z wdobler $")
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
    subroutine initial_condition_ss(f)
!
!  Initial condition for temperature
!
      use Initcond, only: gaunoise
      use Boundcond, only: update_ghosts
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: deltaT
!
!  Amplitude of fluctuation
!
      select case (initTT) 
        case ('single-mode')
        do n=n1,n2
          do m=m1,m2
            f(l1:l2,m,n,iTT) = f(l1:l2,m,n,iTT) + &
                 ampltt*cos(kx_TT*x(l1:l2)/Lxyz(1))*sin(kz_TT*z(n)/Lxyz(3))
          enddo
        enddo
!          
        case ('gaussian-noise')  
          call gaunoise(ampltt,f,iTT)
        case ('nothing')
!
        case default
!
!  Catch unknown values.
!
          write(unit=errormsg,fmt=*) 'No such value initTT = ',trim(initTT)
          call fatal_error('initial_condition_ss',errormsg)
!
      endselect
!
!  Add to that the temperature gradient
!
      deltaT=(Tupp-Tbot)/(n2-n1)
      do n=n1,n2
        do m=m1,m2
          f(l1:l2,m,n,iTT) = f(l1:l2,m,n,iTT) + Tbot + (n-n1)*deltaT
        enddo
      enddo
!
!  Is it needed to apply boundaries here? The code will do it already,
!  unless it is needed for the streamfunction iteration.       
!
      call update_ghosts(f,iTT)
!
      print*,'initial condition for temperature set'
      print*,'sum,max,min TT',sum(   f(l1:l2,m1:m2,n1:n2,iTT)),&
                              maxval(f(l1:l2,m1:m2,n1:n2,iTT)),&
                              minval(f(l1:l2,m1:m2,n1:n2,iTT))      
!
    endsubroutine initial_condition_ss
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
