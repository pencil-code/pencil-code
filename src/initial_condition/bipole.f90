! $Id: bipole.f90
! 
!
!  26-apr-12/piyali: adapted from fieldloop.f90
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
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../initial_condition.h'
!
!!  integer :: dummy
!
  real :: amplaa=1e-3,r0=0.2,beta=0.0,rho0
  character (len=8) :: extfieldfile
!
  namelist /initial_condition_pars/ amplaa, r0, beta, extfieldfile,rho0
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
         "$Id: fieldloop.f90 17229 2011-07-19 20:14:00Z sven.bingert $")
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
    subroutine initial_condition_lnrho(f)
!
!  Initialize the magnetic vector potential.
!
!  02-may-12/piyali: coded
!
      use Sub, only: curl
      use IO,  only: input_snap,input_snap_finalize
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx,3) :: bb
      real :: xi
      integer :: l
!
      call initial_condition_aa(f) 
      do m=m1,m2
        do n=n1,n2
          call curl(f,iaa,bb)
          f(l1:l2,m,n,ilnrho)=log(rho0+bb(:,1)**2+bb(:,2)**2+bb(:,3)**2)
        end do
      end do
    endsubroutine initial_condition_lnrho  
!
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  02-may-12/piyali: coded
!
      use Sub, only: curl
      use IO,  only: input_snap,input_snap_finalize
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,3) :: bbext
      real :: xi
      integer :: l
!
!  IMPLEMENTATION OF INSERTION OF BIPOLES (Non-Potential part) 
!  (Yeates, Mackay and van Ballegooijen 2008, Sol. Phys, 247, 103)
!
      do n=1,mz
        do m=1,my
          do l=1,mx
!
            if (lcartesian_coords) then 
              xi=((x(l)**2+z(n)**2)/2.+y(m)**2)/r0**2
              f(l,m,n,iax) = amplaa*beta*z(n)*exp(-2*xi)
              f(l,m,n,iay) = amplaa*r0*exp(-xi)
              f(l,m,n,iaz) = -amplaa*beta*x(l)*exp(-2*xi)
            else if (lcylindrical_coords) then 
              call fatal_error('initial_condition_aa','Bipoles not coded  &
              for cylindrical coordinates')
            else if (lspherical_coords) then 
!
!  First set up the external potential field
!
              f(l,m,n,iax) = 0.0
              f(l,m,n,iay) = 0.0
              f(l,m,n,iaz) = 0.0
            endif
          enddo
!          call input_snap(extfieldfile,bbext,1,0)
!          call input_persistent(extfieldfile)
!          call input_snap_finalize()
!          f(l,m,n,iglobal_bx_ext) = bbext(l,m,n,1)
!          f(l,m,n,iglobal_by_ext) = bbext(l,m,n,2)
!          f(l,m,n,iglobal_bz_ext) = bbext(l,m,n,3)

!
!  Then set up the helical field
!
          if (lspherical_coords) then 
            do l=1, mx
              xi=(((x(l)-xyz0(1))**2+z(n)**2)/2.+log(tan(y(m)/2.))**2)/r0**2
              f(l,m,n,iax) = -amplaa*beta*z(n)*exp(-2*xi)
              f(l,m,n,iay) = amplaa*r0*exp(-xi)/x(l)
              f(l,m,n,iaz) = amplaa*beta*(x(l)-xyz0(1))*exp(-2*xi)/x(l)
            enddo
          endif
        enddo
      enddo
!
    endsubroutine initial_condition_aa
!***********************************************************************
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
