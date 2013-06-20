! $Id: khi_colin.f90,v 1.1 2011-02-02 20:47:30 wlyra Exp $
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
  use General, only: keep_compiler_quiet
  use EquationOfState
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: lsmooth=0.025
!
! 
  namelist /initial_condition_pars/ lsmooth
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
         "$Id: khi_colin.f90,v 1.1 2011-02-02 20:47:30 wlyra Exp $")
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
      real, dimension (mx,my,mz,mfarray) :: f
      real :: u1,u2
!
      u1= 0.5
      u2=-0.5
!
      do n=n1,n2
        do m=m1,m2 
!
          if (y(m)>0 .and. y(m)<=0.25) then 
            f(l1:l2,m,n,iux) = u1 - .5*(u1-u2)*exp( ( y(m)-0.25)/lsmooth)  
          else if (y(m)>0.25 .and. y(m)<=0.50) then 
            f(l1:l2,m,n,iux) = u2 + .5*(u1-u2)*exp( (-y(m)+0.25)/lsmooth)
          else if (y(m)>0.50 .and. y(m)<=0.75) then 
            f(l1:l2,m,n,iux) = u2 + .5*(u1-u2)*exp(-(0.75 -y(m))/lsmooth)
          else if (y(m)>0.75 .and. y(m)<=1.00) then 
            f(l1:l2,m,n,iux) = u1 - .5*(u1-u2)*exp(-(y(m) -0.75)/lsmooth)
          endif
!
          f(l1:l2,m,n,iuy) = 0.01*sin(4*pi*x(l1:l2))
!
        enddo
      enddo
!
      
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
      use Mpicomm,         only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: rho1,rho2
!
      rho1=1.
      rho2=2.
!
      do n=n1,n2
        do m=m1,m2 
!
          if (y(m)>0 .and. y(m)<=0.25) then 
            f(l1:l2,m,n,ilnrho) = rho1 - .5*(rho1-rho2)*exp( ( y(m)-0.25)/lsmooth)  
          else if (y(m)>0.25 .and. y(m)<=0.50) then 
            f(l1:l2,m,n,ilnrho) = rho2 + .5*(rho1-rho2)*exp( (-y(m)+0.25)/lsmooth)
          else if (y(m)>0.50 .and. y(m)<=0.75) then 
            f(l1:l2,m,n,ilnrho) = rho2 + .5*(rho1-rho2)*exp(-(0.75 -y(m))/lsmooth)
          else if (y(m)>0.75 .and. y(m)<=1.00) then 
            f(l1:l2,m,n,ilnrho) = rho1 - .5*(rho1-rho2)*exp(-(y(m) -0.75)/lsmooth)
          endif
!
        enddo
      enddo
!
      f(l1:l2,m1:m2,n1:n2,ilnrho) = alog(f(l1:l2,m1:m2,n1:n2,ilnrho))
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

      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: lnrho,lnTT,TT,rho
      real :: cp,cv,cp1,lnTT0,pp0,TT0
      integer :: irho
!
!  SAMPLE IMPLEMENTATION
!
      cp=1.
      cp1=1/cp
      cv=gamma1*cp
!
      TT0 = cs20*cp1/gamma_m1 ; lnTT0=log(TT0)
      pp0=(cp-cv)*TT0*rho0
!
      irho=ilnrho
!
      do m=m1,m2;do n=n1,n2
!
        rho = f(l1:l2,m,n,ilnrho)
        TT=(pp0/((cp-cv)*rho))/TT0
!
        lnTT = log(TT)
        lnrho=log(rho)
!
! cp=1
!
        f(l1:l2,m,n,iss)=f(l1:l2,m,n,iss) + &
             cv*(lnTT-gamma_m1*lnrho)
!
      enddo;enddo
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
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
